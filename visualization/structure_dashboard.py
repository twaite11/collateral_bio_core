"""
Generate Structure Dashboard HTML: 3D viewer + domain coloring + homology table.
Reads: input FASTA, InterPro domains.json, homology_scores.json, OmegaFold PDB paths.
Optional: synthesis_metadata.csv or family_manifest for SRA + CRISPR repeats.
Predicts NUC/REC lobes for structure display.
Serve the HTML via HTTP (e.g. python -m http.server) to load PDBs.
"""
import argparse
import csv
import json
import re
from pathlib import Path
from typing import List, Tuple

from Bio import SeqIO

HEPN_REGEX = re.compile(r"R.{4,6}H")

# Pfam: HEPN = nuclease lobe; HEL/WYL = recognition lobe
NUC_SIGNATURES = {"PF05168"}  # HEPN
REC_SIGNATURES = {"PF01228", "PF18456"}  # HEL, WYL


def count_hepn(seq_str: str) -> int:
    return len(HEPN_REGEX.findall(str(seq_str)))


def get_hepn_positions(seq_str: str) -> List[Tuple[int, int]]:
    """Return (start, end) 0-based for each HEPN motif."""
    return [(m.start(), m.end()) for m in HEPN_REGEX.finditer(str(seq_str))]


def predict_nuc_rec_spans(
    sequence: str,
    domains: list,
    hepn_padding: int = 50,
) -> dict:
    """
    Predict NUC (nuclease) and REC (recognition) lobe residue spans.
    NUC = HEPN domains (PF05168) or expanded HEPN motifs if no InterPro.
    REC = HEL/WYL domains or remaining regions.
    Returns {"nuc_spans": [(start,end),...], "rec_spans": [(start,end),...], "nuc_summary": "100-250, 400-550"}.
    """
    nuc_spans = []
    rec_spans = []
    seq_len = len(sequence)

    if domains:
        for d in domains:
            sig = (d.get("signature") or "").split(".")[0]
            start = int(d.get("start", 0))
            end = int(d.get("end", 0))
            if sig in NUC_SIGNATURES:
                nuc_spans.append((start, end))
            elif sig in REC_SIGNATURES:
                rec_spans.append((start, end))

    if not nuc_spans:
        # Fallback: expand HEPN motifs ±hepn_padding aa
        for start, end in get_hepn_positions(sequence):
            s = max(0, start - hepn_padding)
            e = min(seq_len, end + hepn_padding)
            nuc_spans.append((s, e))

    # Sort and merge: keep each HEPN domain separate (do not merge nuc_spans so count matches visual)
    def merge_spans(spans: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        if not spans:
            return []
        spans = sorted(spans, key=lambda x: x[0])
        merged = [spans[0]]
        for s, e in spans[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        return merged

    nuc_spans = sorted(nuc_spans, key=lambda x: x[0])  # do not merge; one span per HEPN domain
    rec_spans = merge_spans(rec_spans)

    # REC = regions not in NUC (if we have InterPro REC domains, use them; else infer)
    if not rec_spans and nuc_spans:
        nuc_set = set()
        for s, e in nuc_spans:
            nuc_set.update(range(s, e))
        pos = 0
        while pos < seq_len:
            if pos not in nuc_set:
                end = pos
                while end < seq_len and end not in nuc_set:
                    end += 1
                rec_spans.append((pos, end))
                pos = end
            else:
                pos += 1
        rec_spans = merge_spans(rec_spans)

    nuc_summary = ", ".join(f"{s+1}-{e}" for s, e in nuc_spans) if nuc_spans else "-"
    rec_summary = ", ".join(f"{s+1}-{e}" for s, e in rec_spans) if rec_spans else "-"

    return {
        "nuc_spans": nuc_spans,
        "rec_spans": rec_spans,
        "nuc_summary": nuc_summary,
        "rec_summary": rec_summary,
    }


def parse_family_id(seq_id: str) -> str:
    parts = str(seq_id).split("_")
    return parts[0] if parts else seq_id


def functional_score(criteria: dict) -> tuple:
    """Composite (pass_count, tiebreaker) for ranking; higher is better."""
    pass_count = 0.0
    tiebreaker = 0.0
    if criteria.get("hepn1_tm") is not None:
        if criteria["hepn1_tm"] >= 0.7:
            pass_count += 1
        tiebreaker += criteria["hepn1_tm"]
    if criteria.get("hepn2_tm") is not None:
        if criteria["hepn2_tm"] >= 0.7:
            pass_count += 1
        tiebreaker += criteria["hepn2_tm"]
    d = criteria.get("catalytic_distance_angstrom")
    if d is not None:
        if d < 25:
            pass_count += 1
        tiebreaker += max(0, 25 - d) / 25.0
    if criteria.get("plddt_dip_ok") is True:
        pass_count += 1
        tiebreaker += 1
    if (criteria.get("linker_net_charge") or 0) > 0:
        tiebreaker += 0.5
    return (pass_count, tiebreaker)


def load_synthesis_metadata(proj: Path, metadata_path: str = None) -> dict:
    """Load new_id -> {sra_accession, repeat_domains} from synthesis_metadata or family_manifest."""
    lookup = {}
    paths = []
    if metadata_path:
        p = proj / metadata_path if not Path(metadata_path).is_absolute() else Path(metadata_path)
        if p.exists():
            paths.append(p)
    if not paths:
        mined = proj / "data" / "mined_sequences"
        if mined.exists():
            paths.extend(sorted(mined.glob("synthesis_metadata.csv")))
            paths.extend(sorted(mined.glob("family_manifest_*.csv"))[-1:])  # latest manifest
    for p in paths:
        try:
            with open(p, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    nid = row.get("new_id", "").strip()
                    if nid:
                        lookup[nid] = {
                            "sra_accession": row.get("sra_accession", "").strip(),
                            "repeat_domains": row.get("repeat_domains", "").strip(),
                        }
        except Exception as e:
            print(f"[!] Warning: could not read metadata {p}: {e}")
    return lookup


def find_pdb_paths(structures_dir: Path, proj: Path) -> dict:
    """Scan OmegaFold (flat *.pdb) or ColabFold (subdirs) output for {seq_id: relative_path}."""
    result = {}
    # OmegaFold: flat *.pdb in output dir
    for pdb in structures_dir.glob("*.pdb"):
        rel = (structures_dir / pdb.name).relative_to(proj)
        result[pdb.stem] = str(rel).replace("\\", "/")
    # ColabFold fallback: subdirs with *_rank_001*.pdb
    if not result:
        for sub in structures_dir.iterdir():
            if not sub.is_dir():
                continue
            for pdb in sub.glob("*_rank_001*.pdb"):
                name = pdb.stem.replace("_unrelaxed", "").replace("_relaxed", "").split("_rank")[0]
                rel = pdb.relative_to(proj)
                result[name] = str(rel).replace("\\", "/")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate structure dashboard with 3D viewer and homology table."
    )
    parser.add_argument(
        "--input",
        default="data/structure_pipeline/passed_structures.fasta",
        help="Input FASTA (pipeline: passed_structures.fasta; or input_2-3_hepn.fasta)",
    )
    parser.add_argument(
        "--domains",
        default="data/structure_pipeline/interpro/domains.json",
        help="InterPro domains JSON",
    )
    parser.add_argument(
        "--homology",
        default="data/structure_pipeline/structures/homology_scores.json",
        help="Homology scores JSON (pipeline writes to structures/ next to PDBs)",
    )
    parser.add_argument(
        "--filter-results",
        default="data/structure_pipeline/structures/structure_filter_results.json",
        help="Structure filter results JSON (pipeline writes to structures/ next to PDBs)",
    )
    parser.add_argument(
        "--structures-dir",
        default="data/structure_pipeline/structures/omegafold",
        help="OmegaFold output directory",
    )
    parser.add_argument(
        "--references-dir",
        default="data/structure_pipeline/references",
        help="Reference PDB directory (5W1H, 6DTD, 6IV9) for side-by-side viewing",
    )
    parser.add_argument(
        "--metadata",
        default=None,
        help="synthesis_metadata.csv or family_manifest CSV (new_id, sra_accession, repeat_domains); searches mined_sequences/ if not set",
    )
    parser.add_argument(
        "--output",
        default="visualization/structure_dashboard.html",
        help="Output HTML path",
    )
    args = parser.parse_args()

    proj = Path(__file__).resolve().parent.parent
    input_path = proj / args.input
    domains_path = proj / args.domains
    homology_path = proj / args.homology
    filter_results_path = proj / args.filter_results
    references_dir = (proj / args.references_dir).resolve()
    struct_dir = proj / args.structures_dir
    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = proj / out_path

    # Load sequences
    seqs = {}
    if input_path.exists():
        for rec in SeqIO.parse(str(input_path), "fasta"):
            seq_str = str(rec.seq)
            seqs[rec.id] = {
                "sequence": seq_str,
                "length": len(seq_str),
                "hepn_count": count_hepn(seq_str),
                "family": parse_family_id(rec.id),
            }
    else:
        seqs = {}

    # Load synthesis metadata (SRA, CRISPR repeats)
    synth_meta = load_synthesis_metadata(proj, args.metadata)

    # Load domains
    domains = {}
    if domains_path.exists():
        with open(domains_path) as f:
            domains = json.load(f)

    # Homology: pipeline writes to structures/homology_scores.json; fallback to structure_pipeline/ or from filter_results
    homology = {}
    if homology_path.exists():
        with open(homology_path) as f:
            homology = json.load(f)
    if not homology and (proj / "data/structure_pipeline/homology_scores.json").exists():
        with open(proj / "data/structure_pipeline/homology_scores.json") as f:
            homology = json.load(f)
    if not homology and filter_results_path.exists():
        with open(filter_results_path) as f:
            raw = json.load(f)
        for sid, r in raw.items():
            h = {k: r[k] for k in ("cas13a", "cas13b", "cas13d") if k in r and r[k] is not None}
            if h:
                homology[sid] = h

    # Filter results: pipeline writes to structures/structure_filter_results.json; fallback to structure_pipeline/
    filter_results = {}
    if filter_results_path.exists():
        with open(filter_results_path) as f:
            filter_results = json.load(f)
    if not filter_results and (proj / "data/structure_pipeline/structure_filter_results.json").exists():
        with open(proj / "data/structure_pipeline/structure_filter_results.json") as f:
            filter_results = json.load(f)

    # PDB paths (relative for fetch from project root)
    pdb_paths = find_pdb_paths(struct_dir, proj) if struct_dir.exists() else {}

    # Build rows for table
    all_ids = sorted(set(seqs.keys()) | set(domains.keys()) | set(homology.keys()) | set(pdb_paths.keys()) | set(filter_results.keys()))
    rows = []
    for seq_id in all_ids:
        s = seqs.get(seq_id, {})
        d = domains.get(seq_id, [])
        h = homology.get(seq_id, {})
        fr = filter_results.get(seq_id, {})
        cas13a = h.get("cas13a")
        cas13b = h.get("cas13b")
        cas13d = h.get("cas13d")
        best = None
        best_name = None
        if cas13a is not None and (best is None or (cas13a or 0) > best):
            best = cas13a
            best_name = "Cas13a"
        if cas13b is not None and (best is None or (cas13b or 0) > best):
            best = cas13b
            best_name = "Cas13b"
        if cas13d is not None and (best is None or (cas13d or 0) > best):
            best = cas13d
            best_name = "RfxCas13d"

        # NUC/REC lobe prediction
        seq_str = s.get("sequence", "")
        nuc_rec = predict_nuc_rec_spans(seq_str, d) if seq_str else {}

        # Functional criteria (soft metrics for ranking)
        hepn1_tm = fr.get("hepn1_tm")
        hepn2_tm = fr.get("hepn2_tm")
        catalytic_dist = fr.get("catalytic_distance_angstrom")
        linker_charge = fr.get("linker_net_charge")
        linker_plddt = fr.get("linker_plddt_mean")
        plddt_dip = fr.get("plddt_dip_ok")
        func_pass, func_tie = functional_score(fr)

        meta = synth_meta.get(seq_id, {})
        nuc_spans_list = nuc_rec.get("nuc_spans", [])
        hepn_count = len(nuc_spans_list)
        if hepn_count == 0:
            continue
        rows.append({
            "id": seq_id,
            "length": s.get("length", 0),
            "hepn_count": hepn_count,
            "family": s.get("family", parse_family_id(seq_id)),
            "domains": d,
            "cas13a": round(cas13a * 100, 1) if cas13a is not None else None,
            "cas13b": round(cas13b * 100, 1) if cas13b is not None else None,
            "cas13d": round(cas13d * 100, 1) if cas13d is not None else None,
            "best_match": best_name,
            "best_pct": round((best or 0) * 100, 1) if best is not None else None,
            "pdb_path": pdb_paths.get(seq_id),
            "sra_accession": meta.get("sra_accession", ""),
            "repeat_domains": meta.get("repeat_domains", ""),
            "nuc_summary": nuc_rec.get("nuc_summary", "-"),
            "rec_summary": nuc_rec.get("rec_summary", "-"),
            "nuc_spans": nuc_rec.get("nuc_spans", []),
            "rec_spans": nuc_rec.get("rec_spans", []),
            "hepn1_tm": round(hepn1_tm, 3) if hepn1_tm is not None else None,
            "hepn2_tm": round(hepn2_tm, 3) if hepn2_tm is not None else None,
            "catalytic_distance_angstrom": round(catalytic_dist, 1) if catalytic_dist is not None else None,
            "linker_net_charge": linker_charge,
            "plddt_dip_ok": plddt_dip,
            "linker_plddt_mean": linker_plddt,
            "func_pass": func_pass,
            "func_tie": func_tie,
        })

    # Sort by functional criteria first (best pass count, then tiebreaker), then homology, then id
    rows.sort(key=lambda r: (
        -(r["func_pass"] or 0),
        -(r["func_tie"] or 0),
        -(r["best_pct"] or 0),
        r["id"],
    ))

    # Reference PDBs for side-by-side viewing (path relative to project root); NUC/REC from sequence
    # Dynamically discover ALL .pdb files in the references directory (supports full Cas13 family)
    refs = []
    try:
        from modules.structure_filter.functional_criteria import get_sequence_from_pdb_single_chain
    except Exception:
        get_sequence_from_pdb_single_chain = None

    # Load catalog if available (from download_cas13_references.py) for subtype labels
    _catalog_subtype: dict = {}
    _catalog_path = references_dir / "cas13_reference_catalog.json"
    if _catalog_path.exists():
        try:
            with open(_catalog_path) as _cf:
                for _entry in json.load(_cf):
                    pid = _entry.get("pdb_id", "").upper()
                    sub = _entry.get("subtype", "Cas13")
                    title = _entry.get("title", "")
                    if pid:
                        _catalog_subtype[pid] = f"{sub} ({pid})" if sub else pid
        except (json.JSONDecodeError, KeyError):
            pass

    # Fallback labels for the original 3 if no catalog
    _legacy_labels = {"5W1H": "Cas13a (5W1H)", "6DTD": "Cas13b (6DTD)", "6IV9": "RfxCas13d (6IV9)"}

    ref_dirs_to_try = [references_dir]
    alt_ref = (proj / "data/references").resolve()
    if alt_ref.exists() and alt_ref != references_dir:
        ref_dirs_to_try.append(alt_ref)

    # Discover all PDB files across reference directories
    seen_ids: set = set()
    for d in ref_dirs_to_try:
        if not d.is_dir():
            continue
        for pdb_file in sorted(d.glob("*.pdb")):
            pdb_id = pdb_file.stem.upper()
            if pdb_id in seen_ids:
                continue
            seen_ids.add(pdb_id)
            label = _catalog_subtype.get(pdb_id) or _legacy_labels.get(pdb_id) or pdb_id
            try:
                rel = pdb_file.resolve().relative_to(proj)
            except ValueError:
                rel = Path(args.references_dir) / pdb_file.name
            ref_entry = {"id": pdb_id, "label": label, "path": str(rel).replace("\\", "/")}
            if get_sequence_from_pdb_single_chain:
                try:
                    seq, _ = get_sequence_from_pdb_single_chain(str(pdb_file.resolve()))
                    if seq:
                        nuc_rec = predict_nuc_rec_spans(seq, [])
                        ref_entry["nuc_spans"] = nuc_rec.get("nuc_spans", [])
                        ref_entry["rec_spans"] = nuc_rec.get("rec_spans", [])
                except Exception:
                    pass
            refs.append(ref_entry)

    if not refs:
        print(f"   [!] No reference PDBs found. Looked in: {references_dir}")
    else:
        print(f"   [*] Reference PDBs from: {references_dir} ({len(refs)} files)")

    # Domain color map for 3Dmol; NUC/REC for lobe coloring
    domain_colors = {
        "PF05168": "#22d3ee",   # HEPN (nuclease)
        "PF01228": "#f97316",   # HEL
        "PF18456": "#94a3b8",   # WYL
    }
    lobe_colors = {"nuc": "#22d3ee", "rec": "#f97316"}  # NUC=cyan, REC=orange
    hepn_colors = ["#22d3ee", "#a78bfa", "#22c55e"]  # HEPN1=cyan, HEPN2=purple, HEPN3=green

    data_json = json.dumps({
        "rows": rows,
        "refs": refs,
        "domain_colors": domain_colors,
        "lobe_colors": lobe_colors,
        "hepn_colors": hepn_colors,
    })

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Structure Dashboard - 2-3 HEPN Cas13</title>
  <link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&display=swap" rel="stylesheet">
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
  <style>
    :root {{ --bg:#0f172a; --card:#1e3a5f; --accent:#38bdf8; --hepn:#22d3ee; --text:#f8fafc; --muted:#94a3b8; }}
    * {{ box-sizing: border-box; }}
    html, body {{ margin: 0; padding: 0; font-family: 'DM Sans', system-ui; background: var(--bg); color: var(--text); min-height: 100vh; }}
    .dashboard-header {{ padding: 0.75rem 1.25rem; border-bottom: 1px solid #334155; flex-shrink: 0; }}
    .dashboard-header h1 {{ font-size: 1.25rem; margin: 0 0 0.2rem; font-weight: 600; }}
    .dashboard-header .sub {{ color: var(--muted); font-size: 0.8rem; margin: 0; }}
    .layout {{ display: grid; grid-template-columns: 1fr 380px; grid-template-rows: auto 1fr; gap: 0; min-height: calc(100vh - 4rem); }}
    @media (min-width: 1100px) {{ .layout {{ grid-template-columns: minmax(0, 1fr) 420px; }} }}
    @media (max-width: 900px) {{ .layout {{ grid-template-columns: 1fr; grid-template-rows: auto auto auto; }} }}
    .table-col {{ display: flex; flex-direction: column; min-height: 0; padding: 0.75rem 1rem; border-right: 1px solid #334155; background: var(--bg); }}
    @media (max-width: 900px) {{ .table-col {{ border-right: none; border-bottom: 1px solid #334155; }} }}
    .table-col h2 {{ font-size: 0.95rem; color: var(--muted); margin: 0 0 0.5rem; font-weight: 500; }}
    .table-wrap {{ flex: 1; min-height: 0; overflow: auto; }}
    .viewer-col {{ display: flex; flex-direction: column; padding: 0.75rem 1rem; min-height: 0; overflow: auto; }}
    .viewer-card {{ flex-shrink: 0; border-radius: 8px; overflow: hidden; background: #1e293b; border: 1px solid #334155; width: 100%; }}
    .viewer-card-header {{ display: flex; align-items: center; justify-content: space-between; padding: 0.4rem 0.5rem; background: #1e293b; border-bottom: 1px solid #334155; }}
    .viewer-card-header h2 {{ font-size: 0.85rem; color: var(--muted); margin: 0; font-weight: 500; }}
    .viewer-wrap {{ position: relative; width: 100%; height: 400px; background: #1e293b; }}
    .viewer {{ position: absolute; top: 0; left: 0; right: 0; bottom: 0; width: 100%; height: 100%; background: #1e293b; }}
    .popout-btn {{ padding: 0.25rem 0.5rem; font-size: 0.75rem; border-radius: 4px; border: 1px solid #334155; background: var(--card); color: var(--muted); cursor: pointer; font-family: inherit; }}
    .popout-btn:hover {{ border-color: var(--accent); color: var(--text); }}
    .modal-overlay {{ display: none; position: fixed; inset: 0; background: rgba(0,0,0,0.7); z-index: 1000; align-items: center; justify-content: center; padding: 1rem; }}
    .modal-overlay.open {{ display: flex; }}
    .modal-content {{ background: var(--bg); border: 1px solid #334155; border-radius: 8px; overflow: hidden; width: 85vw; max-width: 1200px; max-height: 90vh; }}
    .modal-content .viewer {{ width: 100%; height: 75vh; min-height: 400px; position: relative; }}
    .modal-close {{ position: absolute; top: 0.5rem; right: 0.5rem; z-index: 10; padding: 0.35rem 0.6rem; font-size: 0.8rem; border: none; border-radius: 4px; background: var(--card); color: var(--text); cursor: pointer; }}
    .refs-section {{ display: block; min-height: 2.5rem; }}
    .ref-btn {{ appearance: none; padding: 0.35rem 0.6rem; font-size: 0.78rem; border-radius: 6px; border: 1px solid #334155; background: var(--card); color: var(--text); cursor: pointer; font-family: inherit; display: inline-flex; align-items: center; gap: 0.35rem; min-height: 28px; box-shadow: 0 1px 2px rgba(0,0,0,0.2); }}
    .ref-btn .lobe-dots {{ display: inline-flex; align-items: center; gap: 2px; }}
    .ref-btn .lobe-dot {{ width: 6px; height: 6px; border-radius: 50%; flex-shrink: 0; }}
    .ref-btns-empty {{ font-size: 0.75rem; color: var(--muted); padding: 0.4rem 0; }}
    .viewer-actions {{ margin-top: 0.6rem; }}
    .refs-section h3 {{ font-size: 0.8rem; color: var(--muted); margin: 0 0 0.35rem; font-weight: 500; }}
    .ref-btns {{ display: flex; flex-wrap: wrap; gap: 0.4rem; }}
    .ref-btn:hover {{ border-color: var(--accent); background: #1e3a5f; }}
    .ref-btn.selected {{ border-color: var(--accent); background: #1e3a5f; outline: 1px solid var(--accent); }}
    .load-hint {{ color: var(--muted); font-size: 0.75rem; margin-top: 0.5rem; }}
    #domain-legend {{ font-size: 0.8rem; color: var(--muted); margin-top: 0.4rem; }}
    .domain-chart-wrap {{ margin-top: 0.6rem; display: none; }}
    .domain-chart-wrap.visible {{ display: block; }}
    .domain-chart-wrap h4 {{ font-size: 0.78rem; color: var(--muted); margin: 0 0 0.35rem; font-weight: 500; }}
    .domain-chart-bar {{ position: relative; height: 28px; background: #334155; border-radius: 4px; overflow: hidden; width: 100%; }}
    .domain-chart-bar .segment {{ position: absolute; top: 0; height: 100%; border-radius: 2px; }}
    .domain-chart-bar .segment.rec {{ background: #f97316; }}
    .domain-chart-bar .segment.repeat {{ background: #94a3b8; }}
    .domain-chart-labels {{ display: flex; justify-content: space-between; font-size: 0.7rem; color: var(--muted); margin-top: 0.2rem; }}
    .domain-chart-legend {{ display: flex; flex-wrap: wrap; gap: 0.5rem; margin-top: 0.35rem; font-size: 0.7rem; color: var(--muted); }}
    .domain-chart-legend span {{ display: inline-flex; align-items: center; gap: 0.25rem; }}
    .domain-chart-legend .dot {{ width: 8px; height: 8px; border-radius: 50%; }}
    .domain-chart-bar .segment.hepn {{ border-radius: 2px; }}
    #detail-panel {{ margin-top: 0.5rem; padding: 0.45rem 0.5rem; background: #1e293b; border-radius: 6px; font-size: 0.78rem; color: var(--muted); border: 1px solid #334155; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.8rem; }}
    th, td {{ padding: 0.35rem 0.5rem; text-align: left; border-bottom: 1px solid #334155; }}
    th {{ color: var(--muted); font-weight: 500; position: sticky; top: 0; background: var(--bg); z-index: 1; }}
    tr:hover {{ background: #1e3a5f; }}
    tr.selected {{ background: #1e3a5f; outline: 1px solid var(--accent); }}
    .pct {{ font-variant-numeric: tabular-nums; }}
  </style>
</head>
<body>
  <header class="dashboard-header">
    <h1>Structure Dashboard – 2-3 HEPN Cas13</h1>
    <p class="sub">Select a row to load structure · References below · Serve from project root: python -m http.server 8000</p>
  </header>

  <div class="layout">
    <div class="table-col">
      <h2>Discovered</h2>
      <div class="table-wrap">
      <table>
        <thead>
          <tr>
            <th>ID</th>
            <th>Family</th>
            <th>Len</th>
            <th>HEPN</th>
            <th>NUC</th>
            <th>REC</th>
            <th>HEPN1 TM</th>
            <th>HEPN2 TM</th>
            <th>Cat. dist &#197;</th>
            <th>Linker chg</th>
            <th>pLDDT dip</th>
            <th>Func</th>
            <th>Cas13a %</th>
            <th>Cas13b %</th>
            <th>RfxCas13d %</th>
            <th>Best</th>
            <th>SRA</th>
          </tr>
        </thead>
        <tbody id="tbody"></tbody>
      </table>
      </div>
    </div>
    <div class="viewer-col">
      <div class="viewer-card">
        <div class="viewer-card-header">
          <h2>3D structure</h2>
          <button type="button" class="popout-btn" id="popout-btn" title="Open in larger view">Pop out</button>
        </div>
        <div class="viewer-wrap">
          <div id="molviewer" class="viewer"></div>
        </div>
      </div>
      <div id="viewer-modal" class="modal-overlay" aria-hidden="true">
        <div class="modal-content" style="position:relative;">
          <button type="button" class="modal-close" id="modal-close">Close</button>
          <div id="molviewer-modal" class="viewer"></div>
        </div>
      </div>
      <div class="viewer-actions">
        <div class="refs-section">
          <h3>Reference structures</h3>
          <div id="ref-buttons" class="ref-btns"></div>
        </div>
        <p class="load-hint">Click a row or reference to load</p>
        <div id="domain-legend"></div>
        <div id="domain-chart-wrap" class="domain-chart-wrap">
          <h4>Domain map</h4>
          <div id="domain-chart-bar" class="domain-chart-bar"></div>
          <div class="domain-chart-labels"><span id="domain-chart-n">N</span><span id="domain-chart-c">C</span></div>
          <div class="domain-chart-legend">
            <span><span class="dot" style="background:#22d3ee"></span> HEPN1</span>
            <span><span class="dot" style="background:#a78bfa"></span> HEPN2</span>
            <span><span class="dot" style="background:#22c55e"></span> HEPN3</span>
            <span><span class="dot" style="background:#f97316"></span> REC</span>
            <span><span class="dot" style="background:#94a3b8"></span> CRISPR repeats</span>
          </div>
        </div>
        <div id="detail-panel"></div>
      </div>
    </div>
  </div>

  <script>
    const DATA = {data_json};

    function deselectRefButtons() {{
      document.querySelectorAll(".ref-btn").forEach(b => b.classList.remove("selected"));
    }}

    function selectRef(ref, btn) {{
      document.querySelectorAll("#tbody tr").forEach(tr => tr.classList.remove("selected"));
      deselectRefButtons();
      if (btn) btn.classList.add("selected");
      currentLoad = {{ type: "ref", ref: ref }};
      loadReference(ref);
      document.getElementById("detail-panel").innerHTML = "<b>Reference:</b> " + ref.label;
      document.getElementById("domain-chart-wrap").classList.remove("visible");
    }}

    function loadReference(ref) {{
      const wrap = document.querySelector(".viewer-wrap");
      let viewer = document.getElementById("molviewer");
      const div = document.createElement("div");
      div.className = "viewer";
      div.id = "molviewer";
      if (viewer && viewer.parentNode) viewer.parentNode.replaceChild(div, viewer);
      else if (wrap) wrap.appendChild(div);
      viewer = document.getElementById("molviewer");
      const view = $3Dmol.createViewer(div, {{ backgroundColor: "#1e293b" }});
      const base = window.location.pathname.includes("/visualization/") ? "../" : "./";
      const url = base + (ref.path || "");
      fetch(url).then(res => res.text()).then(pdb => {{
        view.addModel(pdb, "pdb");
        view.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});
        const nucSpans = ref.nuc_spans || [];
        const recSpans = ref.rec_spans || [];
        const hepnColors = DATA.hepn_colors || ["#22d3ee", "#a78bfa", "#22c55e"];
        if (nucSpans.length > 0 || recSpans.length > 0) {{
          nucSpans.forEach(([s,e], i) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: hepnColors[i % hepnColors.length] }} }}));
          recSpans.forEach(([s,e]) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: DATA.lobe_colors.rec }} }}));
          const legend = document.getElementById("domain-legend");
          legend.innerHTML = nucSpans.map((_, i) => "<span style=\\"color:" + hepnColors[i % hepnColors.length] + "\\">&#9679;</span> HEPN" + (i+1)).join(" ") + "  <span style=\\"color:" + DATA.lobe_colors.rec + "\\">&#9679;</span> REC";
        }} else {{
          document.getElementById("domain-legend").innerHTML = "";
        }}
        view.zoomTo();
        view.render();
      }}).catch(() => {{
        view.addLabel("Failed to load reference PDB", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }});
        view.zoomTo();
        view.render();
        document.getElementById("domain-legend").innerHTML = "";
      }});
    }}

    function renderRefButtons() {{
      const section = document.querySelector(".refs-section");
      const container = document.getElementById("ref-buttons");
      container.innerHTML = "";
      if (!(DATA.refs && DATA.refs.length)) {{
        container.classList.add("ref-btns-empty");
        container.innerHTML = "No reference PDBs found. Add 5W1H.pdb, 6DTD.pdb, 6IV9.pdb to <code>data/structure_pipeline/references</code> (or <code>data/references</code>) and regenerate the dashboard.";
        return;
      }}
      container.classList.remove("ref-btns-empty");
      const nuc = (DATA.lobe_colors && DATA.lobe_colors.nuc) ? DATA.lobe_colors.nuc : "#22d3ee";
      const rec = (DATA.lobe_colors && DATA.lobe_colors.rec) ? DATA.lobe_colors.rec : "#f97316";
      DATA.refs.forEach(ref => {{
        const btn = document.createElement("button");
        btn.className = "ref-btn";
        btn.type = "button";
        btn.innerHTML = '<span class="lobe-dots"><span class="lobe-dot" style="background:' + nuc + '" title="NUC"></span><span class="lobe-dot" style="background:' + rec + '" title="REC"></span></span> ' + (ref.label || ref.id || "");
        btn.addEventListener("click", () => selectRef(ref, btn));
        container.appendChild(btn);
      }});
    }}

    function renderTable() {{
      const tbody = document.getElementById("tbody");
      tbody.innerHTML = "";
      DATA.rows.forEach((r, i) => {{
        const tr = document.createElement("tr");
        tr.dataset.id = r.id;
        tr.dataset.index = i;
        const sra = (r.sra_accession || "").slice(0, 12);
        const fmt = (v, d) => (v != null && v !== '') ? v : (d || '-');
        const hepn1 = r.hepn1_tm != null ? r.hepn1_tm.toFixed(2) : '-';
        const hepn2 = r.hepn2_tm != null ? r.hepn2_tm.toFixed(2) : '-';
        const catDist = r.catalytic_distance_angstrom != null ? r.catalytic_distance_angstrom + ' &#197;' : '-';
        const dip = r.plddt_dip_ok === true ? '&#10003;' : (r.plddt_dip_ok === false ? '&#10007;' : '-');
        tr.innerHTML = `<td>${{r.id}}</td><td>${{r.family}}</td><td>${{r.length}}</td><td>${{r.hepn_count}}</td>
          <td title="${{r.nuc_summary || ''}}">${{(r.nuc_summary || '-').slice(0,20)}}${{(r.nuc_summary || '').length > 20 ? '…' : ''}}</td>
          <td title="${{r.rec_summary || ''}}">${{(r.rec_summary || '-').slice(0,20)}}${{(r.rec_summary || '').length > 20 ? '…' : ''}}</td>
          <td class="pct" title="HEPN1 local TM &gt;= 0.7">${{hepn1}}</td>
          <td class="pct" title="HEPN2 local TM &gt;= 0.7">${{hepn2}}</td>
          <td class="pct" title="Distance between RxxxxH His &lt; 25 &#197;">${{catDist}}</td>
          <td title="Linker net charge (R+K-D-E)">${{fmt(r.linker_net_charge)}}</td>
          <td title="Linker pLDDT dip (flexible)">${{dip}}</td>
          <td title="Functional criteria pass count (max 4)">${{r.func_pass != null ? r.func_pass : '-'}}</td>
          <td class="pct">${{r.cas13a != null ? r.cas13a + '%' : '-'}}</td>
          <td class="pct">${{r.cas13b != null ? r.cas13b + '%' : '-'}}</td>
          <td class="pct">${{r.cas13d != null ? r.cas13d + '%' : '-'}}</td>
          <td>${{r.best_match || '-'}}</td>
          <td title="${{r.sra_accession || ''}}">${{sra || '-'}}</td>`;
        tr.addEventListener("click", () => selectRow(i));
        tbody.appendChild(tr);
      }});
    }}

    function selectRow(index) {{
      document.querySelectorAll("#tbody tr").forEach(tr => tr.classList.remove("selected"));
      deselectRefButtons();
      const tr = document.querySelector(`#tbody tr[data-index="${{index}}"]`);
      if (tr) tr.classList.add("selected");
      const row = DATA.rows[index];
      currentLoad = {{ type: "row", index: index }};
      loadStructure(row);
      const dp = document.getElementById("detail-panel");
      const parts = [];
      if (row.sra_accession) parts.push(`<b>SRA:</b> ${{row.sra_accession}}`);
      if (row.repeat_domains) parts.push(`<b>CRISPR repeats:</b> ${{row.repeat_domains}}`);
      const fc = [];
      if (row.hepn1_tm != null) fc.push(`HEPN1 TM: ${{row.hepn1_tm.toFixed(2)}}`);
      if (row.hepn2_tm != null) fc.push(`HEPN2 TM: ${{row.hepn2_tm.toFixed(2)}}`);
      if (row.catalytic_distance_angstrom != null) fc.push(`Cat. dist: ${{row.catalytic_distance_angstrom}} &#197;`);
      if (row.linker_net_charge != null) fc.push(`Linker chg: ${{row.linker_net_charge}}`);
      if (row.plddt_dip_ok != null) fc.push(`pLDDT dip: ${{row.plddt_dip_ok ? 'yes' : 'no'}}`);
      if (fc.length) parts.push(`<b>Functional:</b> ${{fc.join(' | ')}}`);
      if (parts.length) dp.innerHTML = parts.join(" &nbsp;|&nbsp; ");
      else dp.innerHTML = "";
      renderDomainChart(row);
    }}

    function parseRepeatDomains(s) {{
      if (!s || typeof s !== 'string') return [];
      const spans = [];
      const re = /(\d+)\s*-\s*(\d+)/g;
      let m;
      while ((m = re.exec(s)) !== null) {{
        const start = parseInt(m[1], 10) - 1;
        const end = parseInt(m[2], 10);
        if (start >= 0 && end > start) spans.push([start, end]);
      }}
      return spans;
    }}

    function renderDomainChart(row) {{
      const wrap = document.getElementById("domain-chart-wrap");
      const barEl = document.getElementById("domain-chart-bar");
      const len = row.length || 1;
      wrap.classList.add("visible");
      barEl.innerHTML = "";
      const nucSpans = row.nuc_spans || [];
      const recSpans = row.rec_spans || [];
      const repeatSpans = parseRepeatDomains(row.repeat_domains);
      function pct(start, end) {{ return (100 * start / len) + "%"; }}
      function width(start, end) {{ return (100 * (end - start) / len) + "%"; }}
      const hepnColors = DATA.hepn_colors || ["#22d3ee", "#a78bfa", "#22c55e"];
      function addSegment(className, start, end, bgColor) {{
        const seg = document.createElement("div");
        seg.className = "segment " + className;
        seg.style.left = pct(start, end);
        seg.style.width = width(start, end);
        if (bgColor) seg.style.background = bgColor;
        seg.title = (start + 1) + "-" + end;
        barEl.appendChild(seg);
      }}
      nucSpans.forEach(([s, e], i) => addSegment("hepn", s, e, hepnColors[i % hepnColors.length]));
      recSpans.forEach(([s, e]) => addSegment("rec", s, e));
      repeatSpans.forEach(([s, e]) => addSegment("repeat", s, e));
      document.getElementById("domain-chart-n").textContent = "N (1)";
      document.getElementById("domain-chart-c").textContent = "C (" + len + ")";
    }}

    function loadStructure(row) {{
      const wrap = document.querySelector(".viewer-wrap");
      let viewer = document.getElementById("molviewer");
      const div = document.createElement("div");
      div.className = "viewer";
      div.id = "molviewer";
      if (viewer && viewer.parentNode) viewer.parentNode.replaceChild(div, viewer);
      else if (wrap) wrap.appendChild(div);
      viewer = document.getElementById("molviewer");

      const view = $3Dmol.createViewer(div, {{ backgroundColor: "#1e293b" }});

      if (!row.pdb_path) {{
        view.addLabel("No PDB available", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }});
        view.zoomTo();
        view.render();
        return;
      }}

      // Fetch PDB - path relative to project root when served via http.server
      const base = window.location.pathname.includes("/visualization/") ? "../" : "./";
      const url = base + (row.pdb_path || "");
      fetch(url).then(res => res.text()).then(pdb => {{
        view.addModel(pdb, "pdb");
        view.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});

        let legendParts = [];
        const hepnColors = DATA.hepn_colors || ["#22d3ee", "#a78bfa", "#22c55e"];
        if ((row.domains || []).length > 0) {{
          const hepnDomains = (row.domains || []).filter(d => (d.signature || "").split(".")[0] === "PF05168").sort((a,b) => (a.start || 0) - (b.start || 0));
          const otherDomains = (row.domains || []).filter(d => (d.signature || "").split(".")[0] !== "PF05168");
          hepnDomains.forEach((d, i) => view.setStyle({{ resi: d.start + "-" + d.end }}, {{ cartoon: {{ color: hepnColors[i % hepnColors.length] }} }}));
          otherDomains.forEach(d => {{
            const sig = (d.signature || "").split(".")[0];
            view.setStyle({{ resi: d.start + "-" + d.end }}, {{ cartoon: {{ color: DATA.domain_colors[sig] || "#64748b" }} }});
          }});
          const legMap = {{ "PF05168": "HEPN", "PF01228": "HEL (REC)", "PF18456": "WYL (REC)" }};
          if (hepnDomains.length) legendParts.push(hepnDomains.map((_, i) => `<span style="color:${{hepnColors[i % hepnColors.length]}}">&#9679;</span> HEPN${{i+1}}`).join(" "));
          const used = [...new Set(otherDomains.map(d => (d.signature || "").split(".")[0]))];
          legendParts = legendParts.concat(used.map(s => `<span style="color:${{DATA.domain_colors[s] || '#64748b'}}">&#9679;</span> ${{legMap[s] || s}}`));
        }} else if ((row.nuc_spans || []).length > 0 || (row.rec_spans || []).length > 0) {{
          (row.nuc_spans || []).forEach(([s,e], i) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: hepnColors[i % hepnColors.length] }} }}));
          (row.rec_spans || []).forEach(([s,e]) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: DATA.lobe_colors.rec }} }}));
          legendParts = (row.nuc_spans || []).map((_, i) => `<span style="color:${{hepnColors[i % hepnColors.length]}}">&#9679;</span> HEPN${{i+1}}`).concat([`<span style="color:${{DATA.lobe_colors.rec}}">&#9679;</span> REC`]);
        }}

        view.zoomTo();
        view.render();

        const legend = document.getElementById("domain-legend");
        legend.innerHTML = legendParts.join(" ") || "";
      }}).catch(() => {{
        view.addLabel("Failed to load PDB", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }});
        view.zoomTo();
        view.render();
      }});
    }}

    function openPopout() {{
      if (!currentLoad) return;
      const modal = document.getElementById("viewer-modal");
      const container = document.getElementById("molviewer-modal");
      container.innerHTML = "";
      const div = document.createElement("div");
      div.className = "viewer";
      div.style.width = "100%";
      div.style.height = "100%";
      container.appendChild(div);
      const view = $3Dmol.createViewer(div, {{ backgroundColor: "#1e293b" }});
      const base = window.location.pathname.includes("/visualization/") ? "../" : "./";

      function done() {{ view.zoomTo(); view.render(); modal.classList.add("open"); modal.setAttribute("aria-hidden", "false"); }}

      const hepnColors = DATA.hepn_colors || ["#22d3ee", "#a78bfa", "#22c55e"];
      if (currentLoad.type === "ref") {{
        const ref = currentLoad.ref;
        const url = base + (ref.path || "");
        fetch(url).then(res => res.text()).then(pdb => {{
          view.addModel(pdb, "pdb");
          view.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});
          const nucSpans = ref.nuc_spans || [];
          const recSpans = ref.rec_spans || [];
          if (nucSpans.length > 0 || recSpans.length > 0) {{
            nucSpans.forEach(([s,e], i) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: hepnColors[i % hepnColors.length] }} }}));
            recSpans.forEach(([s,e]) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: DATA.lobe_colors.rec }} }}));
          }}
          done();
        }}).catch(() => {{ view.addLabel("Failed to load", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }}); done(); }});
      }} else {{
        const row = DATA.rows[currentLoad.index];
        if (!row.pdb_path) {{ view.addLabel("No PDB available", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }}); done(); return; }}
        const url = base + (row.pdb_path || "");
        fetch(url).then(res => res.text()).then(pdb => {{
          view.addModel(pdb, "pdb");
          view.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});
          if ((row.domains || []).length > 0) {{
            const hepnDomains = (row.domains || []).filter(d => (d.signature || "").split(".")[0] === "PF05168").sort((a,b) => (a.start || 0) - (b.start || 0));
            const otherDomains = (row.domains || []).filter(d => (d.signature || "").split(".")[0] !== "PF05168");
            hepnDomains.forEach((d, i) => view.setStyle({{ resi: d.start + "-" + d.end }}, {{ cartoon: {{ color: hepnColors[i % hepnColors.length] }} }}));
            otherDomains.forEach(d => {{ const sig = (d.signature || "").split(".")[0]; view.setStyle({{ resi: d.start + "-" + d.end }}, {{ cartoon: {{ color: DATA.domain_colors[sig] || "#64748b" }} }}); }});
          }} else if ((row.nuc_spans || []).length > 0 || (row.rec_spans || []).length > 0) {{
            (row.nuc_spans || []).forEach(([s,e], i) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: hepnColors[i % hepnColors.length] }} }}));
            (row.rec_spans || []).forEach(([s,e]) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: DATA.lobe_colors.rec }} }}));
          }}
          done();
        }}).catch(() => {{ view.addLabel("Failed to load PDB", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }}); done(); }});
      }}
    }}

    function closePopout() {{
      document.getElementById("viewer-modal").classList.remove("open");
      document.getElementById("viewer-modal").setAttribute("aria-hidden", "true");
    }}

    document.getElementById("popout-btn").addEventListener("click", openPopout);
    document.getElementById("modal-close").addEventListener("click", closePopout);
    document.getElementById("viewer-modal").addEventListener("click", function(e) {{ if (e.target === this) closePopout(); }});

    renderTable();
    renderRefButtons();
    if (DATA.rows.length) selectRow(0);
  </script>
</body>
</html>
'''

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[SUCCESS] Wrote {out_path}")
    print(f"   [*] Loaded {len(rows)} rows, {len(refs)} reference(s)")
    if not rows:
        print(f"   [!] No rows: set --input (e.g. passed_structures.fasta), --filter-results (structures/structure_filter_results.json), --structures-dir")
    print(f"   [*] Serve via: python -m http.server 8000  (from project root)")
    print(f"   [*] Open: http://localhost:8000/visualization/structure_dashboard.html")


if __name__ == "__main__":
    main()
