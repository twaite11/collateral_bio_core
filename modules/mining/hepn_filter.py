"""
Cas13/Type VI HEPN Post-Filter.

Filters prospector hits (deep_hits_*.fasta) using canonical genomics methods:
- Pfam HEPN (PF05168) HMM: require >= 2 domain hits (E-value <= 1e-5)
- Optional motif check: R.{4,6}H with 100-600 aa spacing
"""
import argparse
import glob
import hashlib
import os
import re
from datetime import datetime
from pathlib import Path

from Bio import SeqIO

# HEPN motif for optional secondary filter (Type VI Cas13 topology)
HEPN_REGEX = re.compile(r"R.{4,6}H")
MIN_HEPN_HITS = 2
E_VALUE_THRESHOLD = 1e-5
HEPN_PFAM_ID = "PF05168"


def _has_hepn_motif_topology(seq_str: str) -> bool:
    """Check for two R.{4,6}H motifs with 100-600 aa spacing (Cas13 topology)."""
    matches = [m.start() for m in HEPN_REGEX.finditer(seq_str)]
    if len(matches) < 2:
        return False
    for i in range(len(matches)):
        for j in range(i + 1, len(matches)):
            dist = matches[j] - matches[i]
            if 100 < dist < 600:
                return True
    return False


def _seq_dedup_key(record) -> str:
    """Hash sequence for deduplication."""
    return hashlib.md5(str(record.seq).encode()).hexdigest()


class HEPNFilter:
    """Filter Cas13-like hits by Pfam HEPN domain and optional motif."""

    def __init__(self, hmm_path: str | None = None, require_motif: bool = False):
        self.hmm_path = hmm_path or os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "data", "pfam", "PF05168.hmm"
        )
        self.require_motif = require_motif
        self._pyhmmer_available = False
        try:
            import pyhmmer
            self._pyhmmer_available = True
        except ImportError:
            pass

    def _run_hmmscan_pyhmmer(self, sequences: list) -> dict[str, list[tuple[float, int, int]]]:
        """Use PyHMMER hmmsearch: HMM vs sequences. Returns {seq_id: [(evalue, start, end), ...]}."""
        import pyhmmer.easel
        import pyhmmer.hmmer
        import pyhmmer.plan7

        alphabet = pyhmmer.easel.Alphabet.amino()
        seqs = []
        for rec in sequences:
            ts = pyhmmer.easel.TextSequence(name=rec.id.encode(), sequence=str(rec.seq))
            seqs.append(ts.digitize(alphabet))
        seq_block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)

        results = {rec.id: [] for rec in sequences}
        with pyhmmer.plan7.HMMFile(self.hmm_path) as hmm_file:
            hmm = next(iter(hmm_file))
            for top_hits in pyhmmer.hmmer.hmmsearch(
                [hmm], seq_block, cpus=1, E=E_VALUE_THRESHOLD
            ):
                for hit in top_hits:
                    name = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
                    hits_list = []
                    for dom in hit.domains:
                        ev = float(dom.i_evalue) if hasattr(dom.i_evalue, "__float__") else dom.i_evalue
                        hits_list.append((ev, getattr(dom.alignment, "hmm_from", 0), getattr(dom.alignment, "hmm_to", 0)))
                    if name in results:
                        results[name] = hits_list
        return results

    def _run_hmmscan_subprocess(self, sequences: list, tmp_fasta: str) -> dict[str, list[tuple[float, int, int]]]:
        """Fallback: run hmmscan via subprocess, parse domtblout."""
        import subprocess
        import tempfile
        SeqIO.write(sequences, tmp_fasta, "fasta")
        domtblout = tmp_fasta + ".domtblout"
        try:
            subprocess.run([
                "hmmscan", "--domtblout", domtblout, "-E", str(E_VALUE_THRESHOLD),
                "--noali", self.hmm_path, tmp_fasta  # hmmdb, seqfile
            ], check=True, capture_output=True, timeout=300)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            raise RuntimeError(
                "HMMER hmmscan not found or failed. Install HMMER or use: pip install pyhmmer"
            ) from e
        results = {rec.id: [] for rec in sequences}
        with open(domtblout) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 22:
                    continue
                target = parts[0]
                seq_id = parts[3]
                i_eval = float(parts[12])
                ali_from, ali_to = int(parts[19]), int(parts[20])
                if seq_id in results:
                    results[seq_id].append((i_eval, ali_from, ali_to))
        os.unlink(domtblout)
        return results

    def _run_hmmscan(self, sequences: list) -> dict[str, list[tuple[float, int, int]]]:
        """Scan sequences against Pfam HEPN HMM. Returns {seq_id: [(evalue, start, end), ...]}."""
        if not os.path.exists(self.hmm_path):
            raise FileNotFoundError(
                f"HEPN HMM not found at {self.hmm_path}. "
                "Run: python utils/fetch_pfam_hepn.py"
            )
        if self._pyhmmer_available:
            return self._run_hmmscan_pyhmmer(sequences)
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tf:
            tmp = tf.name
        try:
            return self._run_hmmscan_subprocess(sequences, tmp)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)

    def filter_candidates(
        self,
        fasta_paths: list[str],
        output_dir: str,
        output_basename: str | None = None,
    ) -> tuple[list, list[dict]]:
        """
        Merge FASTA files, dedupe, run HMM + optional motif filter.
        Returns (passed_records, report_rows).
        """
        all_records = []
        seen_hashes = set()
        for path in fasta_paths:
            for rec in SeqIO.parse(path, "fasta"):
                h = _seq_dedup_key(rec)
                if h not in seen_hashes:
                    seen_hashes.add(h)
                    all_records.append(rec)
        if not all_records:
            print("[*] No sequences to filter.")
            return [], []

        print(f"[*] Loaded {len(all_records)} unique sequences from {len(fasta_paths)} file(s).")
        print(f"[*] Running Pfam HEPN (PF05168) HMM scan (E <= {E_VALUE_THRESHOLD})...")
        hmm_hits = self._run_hmmscan(all_records)

        report_rows = []
        passed = []
        for rec in all_records:
            hits = hmm_hits.get(rec.id, [])
            n_hepn = len(hits)
            evals_str = ";".join(f"{e:.2e}" for e, _, _ in hits[:5])
            if n_hepn >= 5:
                evals_str += "..."
            passed_hmm = n_hepn >= MIN_HEPN_HITS
            seq_str = str(rec.seq)
            passed_motif = _has_hepn_motif_topology(seq_str) if self.require_motif else True
            passed_all = passed_hmm and passed_motif
            report_rows.append({
                "seq_id": rec.id,
                "hepn_hits": n_hepn,
                "e_values": evals_str,
                "passed_hmm": passed_hmm,
                "passed_motif": passed_motif if self.require_motif else "N/A",
                "passed": passed_all,
            })
            if passed_all:
                passed.append(rec)

        print(f"[*] HMM: {sum(1 for r in report_rows if r['passed_hmm'])} with >= {MIN_HEPN_HITS} HEPN hits.")
        if self.require_motif:
            print(f"[*] Motif: {sum(1 for r in report_rows if r['passed_motif'])} with valid R.x4-6.H topology.")
        print(f"[*] Final: {len(passed)} sequences passed all filters.")
        return passed, report_rows

    def run(
        self,
        input_dir: str = "data/raw_sequences",
        glob_pattern: str = "deep_hits_*.fasta",
        output_basename: str | None = None,
    ) -> tuple[str, str]:
        """
        Main entry: glob inputs, filter, write FASTA and report TSV.
        Returns (fasta_path, report_path).
        """
        os.makedirs(input_dir, exist_ok=True)
        paths = sorted(glob.glob(os.path.join(input_dir, glob_pattern)))
        if not paths:
            print(f"[!] No files matching {glob_pattern} in {input_dir}")
            return "", ""

        date_str = datetime.now().strftime("%Y%m%d")
        out_base = output_basename or f"cas13_filtered_{date_str}"
        fasta_path = os.path.join(input_dir, f"{out_base}.fasta")
        report_path = os.path.join(input_dir, f"cas13_filter_report_{date_str}.tsv")

        passed, report_rows = self.filter_candidates(paths, input_dir, out_base)
        if passed:
            SeqIO.write(passed, fasta_path, "fasta")
            print(f"[SUCCESS] Wrote {fasta_path} ({len(passed)} sequences).")
        else:
            print("[*] No sequences passed the filter.")

        if report_rows:
            with open(report_path, "w") as f:
                headers = ["seq_id", "hepn_hits", "e_values", "passed_hmm", "passed_motif", "passed"]
                f.write("\t".join(headers) + "\n")
                for r in report_rows:
                    f.write("\t".join(str(r[k]) for k in headers) + "\n")
            print(f"[*] Report: {report_path}")
        return fasta_path, report_path


def main():
    parser = argparse.ArgumentParser(description="Filter Cas13/Type VI hits by Pfam HEPN domain.")
    parser.add_argument("--input-dir", default="data/raw_sequences", help="Directory containing deep_hits FASTA files")
    parser.add_argument("--glob", default="deep_hits_*.fasta", help="Glob pattern for input files")
    parser.add_argument("--output", default=None, help="Output basename (default: cas13_filtered_YYYYMMDD)")
    parser.add_argument("--require-motif", action="store_true", help="Also require R.x4-6.H motif topology")
    parser.add_argument("--hmm", default=None, help="Path to Pfam HEPN (PF05168) HMM file")
    args = parser.parse_args()

    if os.environ.get("HEPN_REQUIRE_MOTIF", "").lower() in ("1", "true", "yes"):
        args.require_motif = True

    filt = HEPNFilter(hmm_path=args.hmm, require_motif=args.require_motif)
    filt.run(
        input_dir=args.input_dir,
        glob_pattern=args.glob,
        output_basename=args.output,
    )


if __name__ == "__main__":
    main()
