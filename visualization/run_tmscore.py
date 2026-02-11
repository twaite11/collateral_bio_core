"""
Compute TM-score folding homology vs Cas13a/b reference structures.
Downloads references from RCSB, runs US-align (or tmtools), outputs homology_scores.json.
"""
import argparse
import json
import os
import re
import subprocess
import urllib.request
from pathlib import Path
from typing import Optional

# Legacy baseline – always downloaded if nothing else is present.
_LEGACY_REFS = {
    "cas13a": "5W1H",
    "cas13b": "6DTD",
    "cas13d": "6IV9",  # RfxCas13d
}

RCSB_URL = "https://files.rcsb.org/download/{pdb}.pdb"

# Catalog written by download_cas13_references.py (subtype info)
_CATALOG_NAME = "cas13_reference_catalog.json"


def download_reference(pdb_id: str, out_dir: str) -> str:
    """Download PDB from RCSB; return path."""
    out_path = Path(out_dir) / f"{pdb_id}.pdb"
    if out_path.exists():
        return str(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = RCSB_URL.format(pdb=pdb_id.upper())
    urllib.request.urlretrieve(url, str(out_path))
    return str(out_path)


def discover_all_references(ref_dir: Path) -> dict:
    """
    Build {label: pdb_id} from every .pdb in ref_dir.
    If a catalog JSON exists (from download_cas13_references.py), use its
    subtype field for readable labels.  Otherwise fall back to the PDB stem.
    """
    catalog_path = ref_dir / _CATALOG_NAME
    subtype_map: dict = {}
    if catalog_path.exists():
        try:
            with open(catalog_path) as f:
                for entry in json.load(f):
                    pid = entry.get("pdb_id", "").upper()
                    sub = entry.get("subtype", "Cas13")
                    if pid:
                        subtype_map[pid] = sub
        except (json.JSONDecodeError, KeyError):
            pass

    ref_pdbs: dict = {}
    for pdb_file in sorted(ref_dir.glob("*.pdb")):
        pdb_id = pdb_file.stem.upper()
        subtype = subtype_map.get(pdb_id, pdb_id)
        # Make label unique if multiple PDBs share a subtype
        label = f"{subtype}_{pdb_id}".lower()
        ref_pdbs[label] = pdb_id
    return ref_pdbs


def run_usalign(pdb1: str, pdb2: str) -> Optional[float]:
    """
    Run US-align and parse TM-score. Returns TM-score (0-1) or None.
    US-align output contains: "TM-score= 0.xxxx (normalized by ...)"
    """
    try:
        result = subprocess.run(
            ["USalign", pdb1, pdb2],
            capture_output=True,
            text=True,
            timeout=60,
        )
        text = result.stdout or result.stderr or ""
        # Match "TM-score= 0.72" or "TM-score=0.72"
        m = re.search(r"TM-score\s*=\s*([\d.]+)", text)
        if m:
            return float(m.group(1))
    except (subprocess.TimeoutExpired, FileNotFoundError, ValueError):
        pass
    return None


def run_tmtools(pdb1: str, pdb2: str) -> Optional[float]:
    """Try tmtools Python package for TM-score."""
    try:
        import tmtools
        res = tmtools.tm_score(pdb1, pdb2)
        return res.tm_norm_chain1 if hasattr(res, "tm_norm_chain1") else res.tm_score
    except Exception:
        return None


def compute_tm_score(pdb_query: str, pdb_ref: str) -> Optional[float]:
    """Compute TM-score; try US-align first, then tmtools."""
    score = run_usalign(pdb_query, pdb_ref)
    if score is not None:
        return score
    return run_tmtools(pdb_query, pdb_ref)


def main():
    parser = argparse.ArgumentParser(
        description="Compute TM-score homology vs Cas13 references."
    )
    parser.add_argument(
        "--structures-dir",
        default="data/structure_pipeline/structures/omegafold",
        help="Directory with predicted PDB files (OmegaFold: flat *.pdb)",
    )
    parser.add_argument(
        "--references-dir",
        default="data/structure_pipeline/references",
        help="Directory for reference PDBs (all Cas13 family; run download_cas13_references.py first)",
    )
    parser.add_argument(
        "--output",
        default="data/structure_pipeline/homology_scores.json",
        help="Output JSON path",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip reference download (use existing)",
    )
    args = parser.parse_args()

    ref_dir = Path(args.references_dir)
    ref_dir.mkdir(parents=True, exist_ok=True)

    # Ensure at least the legacy 3 are present (backwards-compat)
    if not args.skip_download:
        print("[*] Ensuring baseline reference structures are present...")
        for name, pdb_id in _LEGACY_REFS.items():
            download_reference(pdb_id, str(ref_dir))

    # Discover ALL .pdb files in the references dir (curated + downloaded)
    REF_PDBS = discover_all_references(ref_dir)
    if not REF_PDBS:
        print(f"[!] No reference PDBs found in {ref_dir}. Run download_cas13_references.py first.")
        return
    print(f"[*] Using {len(REF_PDBS)} reference structures from {ref_dir}")

    ref_paths = {k: str(ref_dir / f"{v}.pdb") for k, v in REF_PDBS.items()}
    missing = [k for k, p in ref_paths.items() if not os.path.exists(p)]
    if missing:
        for k in missing:
            print(f"[!] Reference not found: {ref_paths[k]}")
        return

    # Collect predicted PDBs (OmegaFold: flat *.pdb in output dir)
    struct_dir = Path(args.structures_dir)
    query_pdbs = {}
    for pdb in struct_dir.glob("*.pdb"):
        query_pdbs[pdb.stem] = str(pdb.resolve())
    # ColabFold fallback: subdirs with *_rank_001*.pdb
    if not query_pdbs:
        for sub in struct_dir.iterdir():
            if not sub.is_dir():
                continue
            for pdb in sub.glob("*_rank_001*.pdb"):
                name = pdb.stem.replace("_unrelaxed", "").replace("_relaxed", "").split("_rank")[0]
                query_pdbs[name] = str(pdb.resolve())

    if not query_pdbs:
        print("[!] No predicted PDB files found. Run OmegaFold first.")
        scores = {}
    else:
        print(f"[*] Computing TM-scores for {len(query_pdbs)} structures...")
        scores = {}
        for i, (seq_id, pdb_path) in enumerate(query_pdbs.items()):
            scores[seq_id] = {}
            for ref_name, ref_path in ref_paths.items():
                tm = compute_tm_score(pdb_path, ref_path)
                scores[seq_id][ref_name] = round(tm, 4) if tm is not None else None
            if (i + 1) % 50 == 0:
                print(f"   [{i+1}/{len(query_pdbs)}]")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(scores, f, indent=2)
    print(f"[SUCCESS] Wrote {out_path}")


if __name__ == "__main__":
    main()
