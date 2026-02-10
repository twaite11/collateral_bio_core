#!/usr/bin/env python3
"""
Run SRA Cas13 search: Bacteria/Archaea WGS or METAGENOMIC, Library METAGENOMIC.
Uses Magic-BLAST (SRA Toolkit) to search raw reads without downloading.
Saves intact Cas13-like hits (700-1400 aa, 2 HEPN, N/C intact, mandatory CRISPR) as
deep_hits_*.fasta; CRISPR domain sequences are saved in metadata with SRA accession for crRNA binding.

Pagination: The script keeps requesting SRA Run pages from NCBI until (1) it has
collected max_total runs (--max-runs), or (2) a page returns fewer than page_size
results (no more runs). There is no fixed "page limit"—it goes on until NCBI has
no more runs to return or you hit your cap.

Requirements:
  - SRA Toolkit installed (magic-blast on PATH)
  - pip: biopython

Usage:
  python scripts/run_sra_cas13_search.py
  python scripts/run_sra_cas13_search.py --max-runs 1000 --batch-size 25
  SRA_TERM="txid2[ORGN] AND metagenomic" python scripts/run_sra_cas13_search.py

Environment:
  SRA_TERM          - NCBI SRA search query (default: Bacteria|Archaea, WGS|METAGENOMIC, METAGENOMIC lib)
  MAX_SRA_RUNS      - Max Run accessions to process (default 100000)
  RUN_BATCH_SIZE    - Runs per magic-blast batch (default 50)
  REFERENCE_FASTA   - Protein Cas13 reference for back-translation (default data/references/mining_refs.fasta)
  MAGICBLAST_CMD    - magic-blast executable (default magicblast)
  OUTPUT_DIR        - Where to write deep_hits_*.fasta (default data/raw_sequences)
  NUM_THREADS       - threads per magic-blast job (default 4)
  MAX_WORKERS       - parallel batch jobs (default 8; total CPU ~ workers * num_threads, e.g. 32)
"""
from __future__ import annotations

import os
import sys
import argparse

# Project root
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from modules.mining.sra_blast_miner import (
    DEFAULT_SRA_TERM,
    get_all_sra_runs,
    mine_sra_with_magicblast,
    save_discoveries,
)


def main():
    parser = argparse.ArgumentParser(
        description="SRA Cas13 search via Magic-BLAST (no download). Output: deep_hits_*.fasta for structure pipeline."
    )
    parser.add_argument(
        "--sra-term",
        default=os.environ.get("SRA_TERM", DEFAULT_SRA_TERM),
        help="NCBI SRA search query (organism, strategy, library)",
    )
    parser.add_argument(
        "--max-runs",
        type=int,
        default=int(os.environ.get("MAX_SRA_RUNS", "100000")),
        help="Max SRA Run accessions to process",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=int(os.environ.get("RUN_BATCH_SIZE", "50")),
        help="Runs per magic-blast batch",
    )
    parser.add_argument(
        "--reference-fasta",
        default=os.environ.get("REFERENCE_FASTA", os.path.join(ROOT, "data", "references", "mining_refs.fasta")),
        help="Protein Cas13 reference (first seq back-translated for magic-blast)",
    )
    parser.add_argument(
        "--output-dir",
        default=os.environ.get("OUTPUT_DIR", os.path.join(ROOT, "data", "raw_sequences")),
        help="Output directory for deep_hits_*.fasta and metadata",
    )
    parser.add_argument(
        "--magicblast-cmd",
        default=os.environ.get("MAGICBLAST_CMD", "magicblast"),
        help="magic-blast executable name or path",
    )
    parser.add_argument(
        "--num-threads",
        type=int,
        default=int(os.environ.get("NUM_THREADS", "4")),
        help="Threads per magic-blast job",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("MAX_WORKERS", "8")),
        help="Parallel batch jobs (total CPU ~ workers * num-threads; default 8 for 32 cores)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only fetch and print SRA Run list; do not run magic-blast",
    )
    args = parser.parse_args()

    print("[*] SRA Cas13 search (Magic-BLAST, no download)")
    print(f"    SRA term: {args.sra_term[:80]}...")
    print(f"    Max runs: {args.max_runs}, batch size: {args.batch_size}, workers: {args.workers}, threads/job: {args.num_threads}")

    print("[*] Fetching SRA Run accessions...")
    runs = get_all_sra_runs(term=args.sra_term, max_total=args.max_runs, page_size=500)
    print(f"    Got {len(runs)} Run accessions (SRR)")

    if not runs:
        print("[!] No SRA runs found. Check SRA_TERM or NCBI availability.")
        return 1

    if args.dry_run:
        for i, r in enumerate(runs[:20]):
            print(f"      {r}")
        if len(runs) > 20:
            print(f"      ... and {len(runs) - 20} more")
        return 0

    if not os.path.exists(args.reference_fasta):
        print(f"[!] Reference FASTA not found: {args.reference_fasta}")
        return 1

    print("[*] Running Magic-BLAST and filtering (700-1400 aa, 2 HEPN, N/C intact, mandatory CRISPR)...")
    discoveries, metadata = mine_sra_with_magicblast(
        sra_runs=runs,
        reference_fasta=args.reference_fasta,
        output_dir=args.output_dir,
        magicblast_cmd=args.magicblast_cmd,
        run_batch_size=args.batch_size,
        num_threads=args.num_threads,
        max_workers=args.workers,
    )

    if not discoveries:
        print("[*] No Cas13-like hits passed filters.")
        return 0

    fasta_path, csv_path = save_discoveries(discoveries, metadata, args.output_dir)
    print(f"[SUCCESS] Saved {len(discoveries)} hits to {fasta_path}")
    print(f"          Metadata: {csv_path} (includes sra_accession and CRISPR domain sequences for crRNA binding)")
    print("          Use this FASTA in run_full_pipeline.py (design -> structure -> identity -> matchmaker).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
