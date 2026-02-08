"""
Filter FASTA to sequences with 2-3 HEPN motifs (R.{4,6}H).
Output for structure pipeline: InterProScan, ESMFold, ColabFold.
"""
import argparse
import os
import re
from pathlib import Path

from Bio import SeqIO

HEPN_REGEX = re.compile(r"R.{4,6}H")
MIN_HEPN = 2
MAX_HEPN = 3


def count_hepn(seq_str: str) -> int:
    """Count HEPN motifs in sequence."""
    return len(HEPN_REGEX.findall(str(seq_str)))


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA to 2-3 HEPN motif sequences for structure pipeline."
    )
    parser.add_argument(
        "--input",
default="data/mined_sequences/fam_fasta.fasta",
            help="Input FASTA (e.g. fam_fasta.fasta from family_grouper)",
    )
    parser.add_argument(
        "--output",
        default="data/structure_pipeline/input_2-3_hepn.fasta",
        help="Output FASTA path",
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"[!] Input not found: {args.input}")
        return

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    passed = []
    total = 0
    for rec in SeqIO.parse(args.input, "fasta"):
        total += 1
        cnt = count_hepn(str(rec.seq))
        if MIN_HEPN <= cnt <= MAX_HEPN:
            passed.append(rec)

    SeqIO.write(passed, out_path, "fasta")
    print(f"[*] Filtered {total} -> {len(passed)} sequences (2-3 HEPN motifs)")
    print(f"[SUCCESS] Wrote {out_path}")


if __name__ == "__main__":
    main()
