"""
Identity filter: keep only sequences with max pairwise identity to known Cas13
below a threshold (drift goal, e.g. <85% for novel IP).
"""
import os
import sys
import argparse
import csv
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices


def pairwise_identity(seq1: str, seq2: str) -> float:
    """Global alignment identity (fraction of aligned positions that match)."""
    aligner = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"))
    alns = aligner.align(seq1, seq2)
    if not alns:
        return 0.0
    a = alns[0]
    matches = sum(1 for i, j in zip(a[0], a[1]) if i == j and i != "-")
    length = max(len(seq1), len(seq2))
    return matches / length if length else 0.0


def max_identity_to_refs(seq: str, ref_seqs: list) -> tuple:
    """Return (max_identity, reference_name) across ref_seqs (list of (name, seq))."""
    best = 0.0
    best_name = ""
    for name, ref_seq in ref_seqs:
        ident = pairwise_identity(seq, ref_seq)
        if ident > best:
            best = ident
            best_name = name
    return best, best_name


def run_identity_filter(
    input_fasta: str,
    reference_fasta: str,
    output_passed_fasta: str,
    output_metadata_csv: str,
    max_identity: float = 0.85,
) -> int:
    """
    Filter input FASTA: keep only sequences with max identity to any reference < max_identity.
    Writes passed FASTA and metadata CSV (sequence_id, max_identity, reference_name).
    """
    input_path = Path(input_fasta)
    ref_path = Path(reference_fasta)
    if not input_path.exists():
        print(f"[!] Input not found: {input_fasta}")
        return 0
    ref_seqs = []
    if ref_path.exists():
        ref_seqs = [(r.id, str(r.seq)) for r in SeqIO.parse(ref_path, "fasta")]
    else:
        print(f"[!] Reference FASTA not found: {reference_fasta}. All identities will be 0.")

    records = list(SeqIO.parse(input_path, "fasta"))
    passed = []
    metadata = []
    for rec in records:
        seq = str(rec.seq)
        ident, ref_name = max_identity_to_refs(seq, ref_seqs)
        if ident < max_identity:
            passed.append(rec)
        metadata.append((rec.id, f"{ident:.4f}", ref_name or ""))

    Path(output_passed_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(output_metadata_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(output_passed_fasta, "w") as f:
        SeqIO.write(passed, f, "fasta")
    with open(output_metadata_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sequence_id", "max_identity", "closest_reference"])
        for row in metadata:
            w.writerow(row)

    print(f"[+] Identity filter (<{max_identity}): {len(passed)}/{len(records)} passed.")
    return len(passed)


def main():
    parser = argparse.ArgumentParser(description="Filter by max identity to known Cas13 (drift goal).")
    parser.add_argument("--input", default="data/structure_pipeline/passed_structures.fasta",
                        help="Input FASTA (e.g. after structure filter)")
    parser.add_argument("--references", default="data/references/known_cas13.fasta",
                        help="FASTA of known Cas13 (Lwa, Rfx, etc.)")
    parser.add_argument("--output", default="data/identity_filtered/passed.fasta",
                        help="Output FASTA (passed only)")
    parser.add_argument("--metadata", default="data/identity_filtered/identity_metadata.csv",
                        help="Output CSV: sequence_id, max_identity, closest_reference")
    parser.add_argument("--max-identity", type=float, default=0.85,
                        help="Max identity to any reference (default 0.85 for novel IP)")
    args = parser.parse_args()

    n = run_identity_filter(
        args.input,
        args.references,
        args.output,
        args.metadata,
        max_identity=args.max_identity,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
