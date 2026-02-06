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

REF_PDBS = {
    "cas13a": "5W1H",
    "cas13b": "6DTD",
}

RCSB_URL = "https://files.rcsb.org/download/{pdb}.pdb"


def download_reference(pdb_id: str, out_dir: str) -> str:
    """Download PDB from RCSB; return path."""
    out_path = Path(out_dir) / f"{pdb_id}.pdb"
    if out_path.exists():
        return str(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = RCSB_URL.format(pdb=pdb_id.upper())
    urllib.request.urlretrieve(url, str(out_path))
    return str(out_path)


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
        default="data/structure_pipeline/structures/colabfold",
        help="Directory with predicted PDB files (scan subdirs for *_rank_001*.pdb)",
    )
    parser.add_argument(
        "--references-dir",
        default="data/structure_pipeline/references",
        help="Directory for reference PDBs (5W1H, 6DTD)",
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

    if not args.skip_download:
        print("[*] Downloading reference structures...")
        for name, pdb_id in REF_PDBS.items():
            download_reference(pdb_id, str(ref_dir))
            print(f"   [+] {pdb_id}.pdb")

    ref_paths = {k: str(ref_dir / f"{v}.pdb") for k, v in REF_PDBS.items()}
    for k, p in ref_paths.items():
        if not os.path.exists(p):
            print(f"[!] Reference not found: {p}")
            return

    # Collect predicted PDBs (ColabFold structure)
    struct_dir = Path(args.structures_dir)
    query_pdbs = {}
    for sub in struct_dir.iterdir():
        if not sub.is_dir():
            continue
        for pdb in sub.glob("*_rank_001*.pdb"):
            # Extract seq_id from path (e.g. SN01_001_...)
            name = pdb.stem.replace("_unrelaxed", "").replace("_relaxed", "").split("_rank")[0]
            query_pdbs[name] = str(pdb.resolve())
    if not query_pdbs:
        # Fallback: treat structures_dir as flat PDB dir
        for pdb in struct_dir.glob("*.pdb"):
            query_pdbs[pdb.stem] = str(pdb.resolve())

    if not query_pdbs:
        print("[!] No predicted PDB files found. Run ColabFold first.")
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
