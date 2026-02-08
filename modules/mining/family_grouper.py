"""
Family Grouper: Groups mined Cas13d sequences into families by homology (ESM-2)
and HEPN domain count. Outputs FASTA with naming convention SN01_001, SN01_002, etc.
"""
import argparse
import hashlib
import os
import re
import sys
from datetime import datetime
from pathlib import Path

# Ensure project root on path when run as script (same pattern as autonomous_prospector)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from modules.mining.deep_miner_utils import DeepEngine

HEPN_REGEX = re.compile(r"R.{4,6}H")


def count_hepn(seq_str: str) -> int:
    """Count HEPN motifs (R.{4,6}H) in sequence."""
    return len(HEPN_REGEX.findall(str(seq_str)))


def seq_dedup_key(record) -> str:
    """Hash sequence for deduplication."""
    return hashlib.md5(str(record.seq).encode()).hexdigest()


class FamilyGrouper:
    def __init__(
        self,
        input_dir="data/mined_sequences",
        output_dir=None,
        similarity_threshold=0.7,
        prefix="SN",
        glob_pattern="*.fasta",
    ):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir or input_dir)
        self.similarity_threshold = float(
            os.getenv("FAMILY_SIMILARITY_THRESHOLD", str(similarity_threshold))
        )
        self.prefix = os.getenv("FAMILY_PREFIX", prefix)
        self.glob_pattern = glob_pattern
        self.deep_engine = None

    def _load_sequences(self):
        """Load all FASTA files, dedupe by sequence."""
        paths = sorted(self.input_dir.glob(self.glob_pattern))
        if not paths:
            raise FileNotFoundError(
                f"No files matching {self.glob_pattern} in {self.input_dir}"
            )
        records = []
        seen = set()
        for p in paths:
            for rec in SeqIO.parse(p, "fasta"):
                key = seq_dedup_key(rec)
                if key not in seen:
                    seen.add(key)
                    records.append(rec)
        return records

    def _partition_by_hepn(self, records):
        """Partition records by HEPN motif count."""
        partitions = {}
        for rec in records:
            cnt = count_hepn(str(rec.seq))
            if cnt not in partitions:
                partitions[cnt] = []
            partitions[cnt].append(rec)
        return partitions

    def _compute_similarity_matrix(self, records):
        """Compute pairwise cosine similarity via ESM-2 embeddings."""
        if self.deep_engine is None:
            self.deep_engine = DeepEngine()
        sequences = [str(r.seq) for r in records]
        result = self.deep_engine.get_embeddings_batch(sequences)
        if result is None or result[0] is None:
            return None, []
        embeddings, valid_indices = result
        n = len(valid_indices)
        if n == 0:
            return None, []
        # Cosine similarity matrix
        emb = embeddings.to("cpu")
        norm = emb / emb.norm(dim=1, keepdim=True)
        sim = norm @ norm.T
        # Clamp for numerical stability
        sim = sim.clamp(-1.0, 1.0)
        return sim.numpy(), valid_indices

    def _cluster_by_similarity(self, sim_matrix):
        """
        Hierarchical clustering. Distance = 1 - similarity.
        Returns cluster labels (1, 2, 3, ...) for each sequence.
        """
        n = sim_matrix.shape[0]
        if n == 1:
            return [1]
        # Condensed distance matrix (upper triangle)
        dist = 1.0 - sim_matrix
        np.fill_diagonal(dist, 0)
        condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method="average")
        # t = 1 - similarity_threshold; clusters where max distance <= t
        t = 1.0 - self.similarity_threshold
        labels = fcluster(Z, t=t, criterion="distance")
        return labels.tolist()

    def run(self):
        """Run family grouping. Returns (fasta_path, manifest_path)."""
        print(f"[*] Loading sequences from {self.input_dir} ({self.glob_pattern})...")
        records = self._load_sequences()
        if not records:
            raise ValueError("No sequences loaded.")
        print(f"   [+] Loaded {len(records)} unique sequences.")

        partitions = self._partition_by_hepn(records)
        for hepn_cnt, recs in partitions.items():
            print(f"   [+] HEPN count {hepn_cnt}: {len(recs)} sequences")

        all_renamed = []
        manifest_rows = []
        family_counter = 0

        for hepn_cnt in sorted(partitions.keys()):
            recs = partitions[hepn_cnt]
            if len(recs) == 1:
                family_counter += 1
                fid = f"{self.prefix}{family_counter:02d}"
                new_id = f"{fid}_001"
                orig_id = recs[0].id
                recs[0].id = new_id
                recs[0].description = ""
                all_renamed.append(recs[0])
                manifest_rows.append(
                    {
                        "new_id": new_id,
                        "original_id": orig_id,
                        "hepn_count": hepn_cnt,
                        "family_id": fid,
                        "member_index": 1,
                    }
                )
                continue

            # Store original ids before overwriting
            for r in recs:
                r._original_id = r.id

            print(f"[*] Clustering {len(recs)} sequences (HEPN={hepn_cnt})...")
            sim_matrix, valid_indices = self._compute_similarity_matrix(recs)
            skipped_indices = set(range(len(recs))) - set(valid_indices)
            if skipped_indices:
                print(f"   [!] {len(skipped_indices)} sequences skipped (too short or embedding failed). Assigning to own families.")
                for i in sorted(skipped_indices):
                    r = recs[i]
                    family_counter += 1
                    fid = f"{self.prefix}{family_counter:02d}"
                    new_id = f"{fid}_001"
                    orig_id = getattr(r, "_original_id", r.id)
                    r.id = new_id
                    r.description = ""
                    all_renamed.append(r)
                    manifest_rows.append(
                        {
                            "new_id": new_id,
                            "original_id": orig_id,
                            "hepn_count": hepn_cnt,
                            "family_id": fid,
                            "member_index": 1,
                        }
                    )
            if sim_matrix is None or len(valid_indices) == 0:
                continue

            # Map valid_indices back to full list
            sub_recs = [recs[i] for i in valid_indices]
            labels = self._cluster_by_similarity(sim_matrix)

            # Group by cluster
            cluster_to_recs = {}
            for idx, label in enumerate(labels):
                if label not in cluster_to_recs:
                    cluster_to_recs[label] = []
                cluster_to_recs[label].append(sub_recs[idx])

            # Assign IDs: SN01_001, SN01_002, SN02_001, ...
            for _ in sorted(cluster_to_recs.keys()):
                family_counter += 1
                fid = f"{self.prefix}{family_counter:02d}"
                members = cluster_to_recs[_]
                for mi, r in enumerate(members, start=1):
                    new_id = f"{fid}_{mi:03d}"
                    orig = getattr(r, "_original_id", r.id)
                    r.id = new_id
                    r.description = ""
                    all_renamed.append(r)
                    manifest_rows.append(
                        {
                            "new_id": new_id,
                            "original_id": orig,
                            "hepn_count": hepn_cnt,
                            "family_id": fid,
                            "member_index": mi,
                        }
                    )

        # Fix manifest original_id - we overwrote .id, so use the stored _original_id
        # Actually we stored it before overwriting - but the manifest uses getattr(r, "_original_id", r.id)
        # and we've already overwritten r.id. So _original_id should be correct. Good.

        date_str = datetime.now().strftime("%Y%m%d")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Main output: fam_fasta.fasta (used by matchmaker)
        fam_fasta_path = Path("data") / "mined_sequences" / "fam_fasta.fasta"
        fam_fasta_path.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(all_renamed, fam_fasta_path, "fasta")

        # Also write timestamped copies for history
        fasta_path = self.output_dir / f"family_grouped_{date_str}.fasta"
        manifest_path = self.output_dir / f"family_manifest_{date_str}.csv"
        SeqIO.write(all_renamed, fasta_path, "fasta")
        df = pd.DataFrame(manifest_rows)
        df.to_csv(manifest_path, index=False)

        print(f"\n[SUCCESS] Wrote {fam_fasta_path} ({len(all_renamed)} sequences) [matchmaker input].")
        print(f"          Wrote {fasta_path}, {manifest_path}")
        return str(fam_fasta_path), str(manifest_path)


def main():
    parser = argparse.ArgumentParser(description="Group mined sequences into families by homology and HEPN count.")
    parser.add_argument(
        "--input-dir",
        default=os.getenv("FAMILY_INPUT_DIR", "data/mined_sequences"),
        help="Input directory with FASTA files",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: same as input)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=float(os.getenv("FAMILY_SIMILARITY_THRESHOLD", "0.7")),
        help="Cosine similarity threshold for same family",
    )
    parser.add_argument(
        "--prefix",
        default=os.getenv("FAMILY_PREFIX", "SN"),
        help="Family ID prefix (e.g. SN -> SN01, SN02)",
    )
    parser.add_argument(
        "--glob",
        default="*.fasta",
        help="Glob pattern for input files",
    )
    args = parser.parse_args()

    grouper = FamilyGrouper(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        similarity_threshold=args.threshold,
        prefix=args.prefix,
        glob_pattern=args.glob,
    )
    grouper.run()


if __name__ == "__main__":
    main()
