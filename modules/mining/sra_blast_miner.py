"""
SRA Cas13 search using Magic-BLAST (blastn_vdb) against SRA runs without downloading.
Fetches Run accessions (SRR) by ESearch, runs magic-blast -sra, parses hits, applies
full-enzyme filters (700-1400 aa, exactly 2 HEPN, N-term M, C-term tail, mandatory CRISPR).
CRISPR is required: we search 5 kb upstream and 5 kb downstream of the ORF for the nearest
CRISPR repeat and save that sequence in metadata with SRA accession for crRNA binding.
Saves deep_hits FASTA + metadata for the structure pipeline.

Pagination: SRA ESearch is paginated (retstart/retmax). We keep requesting pages until
we have max_total runs or NCBI returns no more. There is no fixed "page limit"—we stop
when the server returns fewer than page_size results or we hit our max_total cap.
"""
from __future__ import annotations

import os
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Optional, Dict

from Bio import Entrez, SeqIO, Seq
from Bio.Seq import Seq

# Project root
_root = Path(__file__).resolve().parents[2]
if str(_root) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(_root))

try:
    from modules.mining.full_orf_checks import passes_n_term, passes_c_term, get_full_orf_config
except ImportError:
    from full_orf_checks import passes_n_term, passes_c_term, get_full_orf_config

Entrez.email = os.environ.get("ENTREZ_EMAIL", "founder@senarybio.com")

HEPN_REGEX = re.compile(r"R.{4,6}H")

# Default SRA search: Bacteria OR Archaea, WGS or METAGENOMIC strategy, METAGENOMIC library source
DEFAULT_SRA_TERM = (
    '(txid2[ORGN] OR txid2157[ORGN]) '
    'AND (wgs[Strategy] OR metagenomic[Strategy]) '
    'AND metagenomic[LibrarySource]'
)

MIN_AA = 700
MAX_AA = 1400
HEPN_COUNT_EXACT = 2


# E. coli / bacterial preferred codons for back-translation (table 11)
_BACK_TABLE = {
    "M": "ATG", "F": "TTT", "L": "CTG", "S": "TCT", "Y": "TAT", "*": "TAA",
    "C": "TGC", "W": "TGG", "P": "CCG", "H": "CAT", "Q": "CAG", "R": "CGT",
    "I": "ATT", "T": "ACT", "N": "AAT", "K": "AAA", "V": "GTG", "A": "GCG",
    "D": "GAT", "E": "GAA", "G": "GGT",
}


def _back_translate(protein: str) -> str:
    """Back-translate protein to nucleotide using bacterial-preferred codons."""
    dna = []
    for aa in protein.upper():
        if aa == "*":
            break
        dna.append(_BACK_TABLE.get(aa, "ATG"))
    return "".join(dna)


def _parse_runs_from_xml(xml_handle) -> List[str]:
    """Parse RUN accessions from SRA EFetch XML."""
    import xml.etree.ElementTree as ET
    runs = []
    try:
        tree = ET.parse(xml_handle)
        root = tree.getroot()
        for elem in root.iter():
            if elem.tag.endswith("RUN") or (elem.tag == "RUN"):
                acc = elem.get("accession") or elem.get("acc")
                if acc and str(acc).startswith("SRR"):
                    runs.append(str(acc))
            if elem.tag.endswith("Run") and elem.text and elem.text.strip().startswith("SRR"):
                runs.append(elem.text.strip())
    except Exception:
        pass
    return runs


def fetch_sra_run_accessions(
    term: str = DEFAULT_SRA_TERM,
    max_records: int = 100_000,
    batch_size: int = 500,
    retstart: int = 0,
) -> List[str]:
    """
    ESearch SRA with term, then EFetch XML to resolve Run accessions (SRR).
    Returns list of SRR accessions. Paginate with retstart for large result sets.
    """
    run_ids = []
    try:
        h = Entrez.esearch(db="sra", term=term, retmax=min(batch_size, max_records), retstart=retstart)
        rec = Entrez.read(h)
        h.close()
        id_list = rec.get("IdList", [])
        if not id_list:
            return run_ids
        # EFetch XML to get RUN accessions (experiments contain one or more runs)
        h = Entrez.efetch(db="sra", id=",".join(id_list), retmode="xml")
        run_ids = _parse_runs_from_xml(h)
        h.close()
        if not run_ids:
            # Fallback: ESummary - some responses include Run list
            h = Entrez.esummary(db="sra", id=",".join(id_list))
            summary = Entrez.read(h)
            h.close()
            for item in (summary if isinstance(summary, list) else [summary]):
                if isinstance(item, dict):
                    for key in ("Runs", "Run", "runs", "run"):
                        if key in item:
                            val = item[key]
                            if isinstance(val, str) and val.startswith("SRR"):
                                run_ids.append(val)
                            elif isinstance(val, list):
                                for r in val:
                                    if isinstance(r, str) and r.startswith("SRR"):
                                        run_ids.append(r)
                                    elif isinstance(r, dict) and "acc" in r:
                                        run_ids.append(r["acc"])
                    acc = item.get("Accession") or item.get("Run") or item.get("acc")
                    if acc and str(acc).startswith("SRR"):
                        run_ids.append(str(acc))
    except Exception as e:
        print(f"[!] SRA fetch error: {e}")
    return list(dict.fromkeys(run_ids))


def get_all_sra_runs(
    term: str = DEFAULT_SRA_TERM,
    max_total: int = 100_000,
    page_size: int = 500,
) -> List[str]:
    """
    Paginate SRA search until we have up to max_total Run accessions.
    Keeps requesting pages (retstart) until NCBI returns no runs (empty batch) or we hit max_total.
    One page of SRA records can yield fewer runs than page_size (e.g. 438 runs from 500 records),
    so we do not stop on len(batch) < page_size—only when batch is empty.
    """
    all_runs = []
    retstart = 0
    while len(all_runs) < max_total:
        batch = fetch_sra_run_accessions(term, max_records=page_size, batch_size=page_size, retstart=retstart)
        if not batch:
            break
        for r in batch:
            if r not in all_runs:
                all_runs.append(r)
        retstart += page_size
        if retstart >= max_total:
            break
    return all_runs[:max_total]


def build_nucleotide_reference(protein_fasta: str, max_aa: int = 600, out_path: Optional[str] = None) -> str:
    """
    Back-translate first max_aa of first protein in FASTA to DNA; write to out_path or temp.
    Returns path to nucleotide FASTA (used as input for makeblastdb).
    """
    path = Path(protein_fasta)
    if not path.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {protein_fasta}")
    rec = next(SeqIO.parse(path, "fasta"), None)
    if not rec:
        raise ValueError(f"No sequence in {protein_fasta}")
    seq = str(rec.seq).strip()
    seq = seq[:max_aa].replace("*", "")
    dna = _back_translate(seq)
    out = out_path or tempfile.mktemp(suffix=".fasta", prefix="cas13_ref_")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write(f">{rec.id}_nt\n{dna}\n")
    return out


def build_blast_db(
    fasta_path: str,
    db_path: str,
    makeblastdb_cmd: str = "makeblastdb",
) -> str:
    """
    Run makeblastdb to create a nucleotide BLAST database for Magic-BLAST -db.
    db_path is the basename/path without extension (e.g. output_dir/_cas13_ref_db).
    Returns db_path.
    """
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        makeblastdb_cmd,
        "-in", str(fasta_path),
        "-dbtype", "nucl",
        "-out", str(db_path),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        if r.stderr:
            raise RuntimeError(f"makeblastdb failed: {r.stderr[:500]}")
        raise RuntimeError("makeblastdb failed")
    return str(db_path)


def build_magicblast_command(
    sra_accessions: List[str],
    db_path: str,
    out_tsv: str,
    magicblast_cmd: str = "magicblast",
    num_threads: int = 4,
) -> Tuple[List[str], str, Optional[str]]:
    """
    Build the Magic-BLAST command and cwd for SRA mode. Returns (cmd, cwd, batch_file_path).
    batch_file_path is None for single accession; otherwise path to write SRA list (caller must create it).
    """
    out_path = Path(out_tsv)
    db_path_p = Path(db_path)
    db_dir = db_path_p.parent.resolve()
    db_basename = db_path_p.name
    if len(sra_accessions) == 1:
        cmd = (
            [magicblast_cmd]
            + ["-sra", sra_accessions[0]]
            + ["-db", db_basename, "-out", out_path.name, "-outfmt", "tabular", "-num_threads", str(num_threads)]
        )
        return cmd, str(db_dir), None
    batch_file = str(out_path.parent / f"_sra_batch_{out_path.stem}.txt")
    cmd = (
        [magicblast_cmd]
        + ["-sra_batch", Path(batch_file).name]
        + ["-db", db_basename, "-out", out_path.name, "-outfmt", "tabular", "-num_threads", str(num_threads)]
    )
    return cmd, str(db_dir), batch_file


def run_magic_blast(
    sra_accessions: List[str],
    db_path: str,
    out_tsv: str,
    magicblast_cmd: str = "magicblast",
    num_threads: int = 4,
    evalue: float = 1e-5,
) -> bool:
    """
    Run magic-blast with -db (BLAST database) and -sra / -sra_batch per Magic-BLAST docs.
    Uses -outfmt tabular. Runs with cwd=DB directory and -db as basename so the DB is found.
    """
    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd, cwd, batch_file = build_magicblast_command(
        sra_accessions, db_path, out_tsv, magicblast_cmd=magicblast_cmd, num_threads=num_threads
    )
    try:
        if batch_file:
            with open(batch_file, "w") as f:
                f.write("\n".join(sra_accessions))
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, cwd=cwd)
        if r.returncode != 0 and r.stderr:
            print(f"[!] magic-blast stderr: {r.stderr[:500]}")
        if r.returncode == 0 and Path(cwd).resolve() != out_path.parent.resolve():
            (Path(cwd) / out_path.name).rename(out_tsv)
        return r.returncode == 0
    finally:
        if batch_file and Path(batch_file).exists():
            try:
                Path(batch_file).unlink()
            except OSError:
                pass


def _translate_hit(sseq: str, sstart: int, send: int) -> str:
    """Translate subject nucleotide hit to protein. sstart/send 1-based; reverse if sstart > send."""
    s = sseq.replace("-", "").upper()
    if not s:
        return ""
    if sstart > send:
        s = str(Seq(s).reverse_complement())
    # Frame from start (1-based -> 0-based frame)
    frame = (abs(sstart) - 1) % 3
    s = s[frame:]
    return str(Seq(s).translate(to_stop=False))


def load_ref_seqs(fasta_path: str) -> Dict[str, str]:
    """Load reference sequences from FASTA into id -> sequence (uppercase). Matches makeblastdb IDs (first token of header)."""
    ref_seqs = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        ref_seqs[rec.id] = str(rec.seq).upper()
    return ref_seqs


def parse_magicblast_tabular(
    tab_path: str,
    ref_seqs: Dict[str, str],
) -> List[Tuple[str, str, str, str, int, int]]:
    """
    Parse Magic-BLAST -outfmt tabular output. Subject sequence is not in the file; we look up
    ref_seqs[ref_id] and extract the aligned segment. Returns list of
    (qseqid, sseqid, protein_sequence, nucleotide_sseq, sstart, send) for filtering.
    Tabular columns: 1=query id, 2=ref id, 9=sstart, 10=send (1-based), 15=ref strand.
    """
    rows = []
    with open(tab_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 15:
                continue
            qseqid, sseqid = parts[0], parts[1]
            try:
                sstart_i, send_i = int(parts[8]), int(parts[9])
            except (ValueError, IndexError):
                continue
            ref_seq = ref_seqs.get(sseqid) or ref_seqs.get(sseqid.split("|")[-1] if "|" in sseqid else "")
            if not ref_seq:
                continue
            # Extract segment in reference order (1-based inclusive)
            lo, hi = min(sstart_i, send_i), max(sstart_i, send_i)
            segment = ref_seq[lo - 1 : hi]
            protein = _translate_hit(segment, sstart_i, send_i)
            if not protein or "*" in protein:
                protein = protein.split("*")[0]
            sseq_nt = segment.replace("-", "").upper()
            if sstart_i > send_i:
                sseq_nt = str(Seq(sseq_nt).reverse_complement())
            rows.append((qseqid, sseqid, protein, sseq_nt, sstart_i, send_i))
    return rows


# CRISPR repeat: min length and min tandem count (required for crRNA binding)
CRISPR_MIN_LEN = 20
CRISPR_MIN_COUNT = 2
# Flanking window (bp) upstream and downstream of ORF to search for CRISPR array
CRISPR_FLANK_BP = 5000


def _orf_bounds_in_hit(nucleotide: str, sstart: int, send: int, protein_len: int) -> Tuple[int, int]:
    """
    Return (orf_start, orf_end) 0-based indices in the gap-stripped nucleotide string.
    Forward: orf is at frame..frame+protein_len*3; reverse: orf is at end of alignment.
    """
    n = nucleotide.replace("-", "").upper()
    frame = (abs(sstart) - 1) % 3
    orf_nt_len = protein_len * 3
    if sstart <= send:
        orf_start = frame
        orf_end = frame + orf_nt_len
    else:
        orf_start = len(n) - (frame + orf_nt_len)
        orf_end = len(n) - frame
    return max(0, orf_start), min(len(n), orf_end)


def _find_repeats_with_positions(
    nucleotide: str, min_len: int = 20, min_count: int = 2
) -> List[Tuple[str, int]]:
    """Find CRISPR-like tandem repeats; return list of (repeat_sequence, start_position)."""
    n = nucleotide.replace("-", "").upper()
    found = []
    seen_motifs = set()
    for L in range(min_len, min(55, len(n) // min_count + 1)):
        for i in range(len(n) - L * min_count + 1):
            motif = n[i : i + L]
            if motif in seen_motifs:
                continue
            count = 0
            j = i
            while j <= len(n) - L and n[j : j + L] == motif:
                count += 1
                j += L
            if count >= min_count:
                seen_motifs.add(motif)
                found.append((motif, i))
    return found


def _nearest_crispr_in_flanking(
    nucleotide: str,
    orf_start: int,
    orf_end: int,
    flank_bp: int = CRISPR_FLANK_BP,
    min_len: int = CRISPR_MIN_LEN,
    min_count: int = CRISPR_MIN_COUNT,
) -> Optional[str]:
    """
    Search ±flank_bp upstream and downstream of the ORF for CRISPR repeats.
    Return the repeat sequence nearest to the ORF (upstream or downstream), or None if none found.
    """
    n = nucleotide.replace("-", "").upper()
    up_start = max(0, orf_start - flank_bp)
    up_end = orf_start
    down_start = orf_end
    down_end = min(len(n), orf_end + flank_bp)
    upstream_region = n[up_start:up_end]
    downstream_region = n[down_start:down_end]
    best_repeat: Optional[str] = None
    best_dist: Optional[int] = None
    for region, orf_bound, dist_from_bound in [
        (upstream_region, orf_start, lambda pos: up_end - up_start - pos),
        (downstream_region, orf_end, lambda pos: pos),
    ]:
        for repeat_seq, pos in _find_repeats_with_positions(region, min_len, min_count):
            d = dist_from_bound(pos)
            if d < 0:
                continue
            if best_dist is None or d < best_dist:
                best_dist = d
                best_repeat = repeat_seq
    return best_repeat


def _extract_crispr_repeat_sequences(nucleotide: str, min_len: int = 20, min_count: int = 2) -> List[str]:
    """
    Extract CRISPR repeat sequences from nucleotide hit: tandem repeats of length
    min_len repeated at least min_count times. Returns list of unique repeat sequences
    (for metadata; needed for crRNA binding / functional protein).
    """
    n = nucleotide.replace("-", "").upper()
    found = []
    seen_motifs = set()
    for L in range(min_len, min(55, len(n) // min_count + 1)):
        for i in range(len(n) - L * min_count + 1):
            motif = n[i : i + L]
            if motif in seen_motifs:
                continue
            count = 0
            j = i
            while j <= len(n) - L and n[j : j + L] == motif:
                count += 1
                j += L
            if count >= min_count:
                seen_motifs.add(motif)
                found.append(motif)
    return found


def _has_crispr_repeat(nucleotide: str, min_len: int = 20, min_count: int = 2) -> bool:
    """True if nucleotide has a CRISPR-like tandem repeat (required for functional crRNA binding)."""
    return len(_extract_crispr_repeat_sequences(nucleotide, min_len, min_count)) > 0


def passes_cas13_filters(
    protein: str,
    min_tail: int = 15,
    require_m: bool = True,
    nucleotide_hit: Optional[str] = None,
    skip_crispr_check: bool = False,
) -> bool:
    """
    Apply filters: length 700-1400 aa, exactly 2 HEPN, N-term M, C-term tail.
    CRISPR is mandatory: either nucleotide_hit contains a repeat (in-hit), or caller
    guarantees a repeat in flanking (skip_crispr_check=True).
    """
    if not protein or len(protein) < MIN_AA or len(protein) > MAX_AA:
        return False
    motifs = list(HEPN_REGEX.finditer(protein))
    if len(motifs) != HEPN_COUNT_EXACT:
        return False
    if not passes_n_term(protein, require_m):
        return False
    if not passes_c_term(protein, min_tail):
        return False
    if not skip_crispr_check and (not nucleotide_hit or not _has_crispr_repeat(nucleotide_hit, CRISPR_MIN_LEN, CRISPR_MIN_COUNT)):
        return False
    return True


# Delimiter for multiple CRISPR repeat sequences in metadata CSV (no commas)
CRISPR_REPEAT_DELIM = "|"


def _process_one_batch(
    batch_index: int,
    batch: List[str],
    db_path: str,
    out_dir: Path,
    ref_seqs: Dict[str, str],
    min_tail: int,
    require_m: bool,
    magicblast_cmd: str,
    num_threads: int,
) -> List[Tuple[str, str, str]]:
    """
    Run magic-blast for one batch, parse, filter. Returns list of (run_id, protein, crispr_str).
    Used by parallel workers; main thread dedupes and assigns seq_id.
    """
    tsv = str(out_dir / f"_magicblast_batch_{batch_index}.tsv")
    ok = run_magic_blast(batch, db_path, tsv, magicblast_cmd=magicblast_cmd, num_threads=num_threads)
    results = []
    if not ok or not Path(tsv).exists():
        return results
    try:
        for qseqid, sseqid, protein, sseq_nt, sstart_i, send_i in parse_magicblast_tabular(tsv, ref_seqs):
            orf_start, orf_end = _orf_bounds_in_hit(sseq_nt, sstart_i, send_i, len(protein))
            nearest_crispr = _nearest_crispr_in_flanking(
                sseq_nt, orf_start, orf_end,
                flank_bp=CRISPR_FLANK_BP,
                min_len=CRISPR_MIN_LEN,
                min_count=CRISPR_MIN_COUNT,
            )
            if not nearest_crispr:
                continue
            if not passes_cas13_filters(
                protein,
                min_tail=min_tail,
                require_m=require_m,
                skip_crispr_check=True,
            ):
                continue
            run_id = sseqid.split(".")[0] if "." in sseqid else sseqid.split("_")[0]
            results.append((run_id, protein, nearest_crispr))
    finally:
        try:
            Path(tsv).unlink()
        except OSError:
            pass
    return results


def mine_sra_with_magicblast(
    sra_runs: List[str],
    reference_fasta: str,
    output_dir: str = "data/raw_sequences",
    magicblast_cmd: str = "magicblast",
    run_batch_size: int = 50,
    num_threads: int = 4,
    max_workers: int = 8,
) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str, str, str]]]:
    """
    For each batch of SRA runs, run magic-blast (optionally in parallel), parse hits,
    filter to full Cas13-like ORFs with mandatory CRISPR repeat.
    Returns (discoveries, metadata).

    max_workers: number of batches to run in parallel (default 8 for ~32 cores with num_threads=4).
    Total CPU use ~ max_workers * num_threads.
    """
    cfg = get_full_orf_config()
    min_tail = cfg["min_tail"]
    require_m = cfg["require_m"]
    ref_nt_path = build_nucleotide_reference(reference_fasta, max_aa=600)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    db_path = str(out_dir / "_cas13_ref_db")
    makeblastdb_cmd = str(Path(magicblast_cmd).parent / "makeblastdb") if os.path.dirname(magicblast_cmd) else "makeblastdb"
    build_blast_db(ref_nt_path, db_path, makeblastdb_cmd=makeblastdb_cmd)
    ref_seqs = load_ref_seqs(ref_nt_path)
    discoveries = []
    metadata = []
    seen_seqs = set()

    batch_indices = list(range(0, len(sra_runs), run_batch_size))
    if max_workers <= 1:
        for i in batch_indices:
            batch = sra_runs[i : i + run_batch_size]
            for run_id, protein, crispr_str in _process_one_batch(
                i, batch, db_path, out_dir, ref_seqs, min_tail, require_m, magicblast_cmd, num_threads
            ):
                key = protein[:200] + "|" + protein[-200:]
                if key in seen_seqs:
                    continue
                seen_seqs.add(key)
                seq_id = f"Cas13_{run_id}_{len(discoveries)}"
                discoveries.append((seq_id, protein))
                metadata.append((seq_id, run_id, crispr_str, "0"))
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    _process_one_batch,
                    i,
                    sra_runs[i : i + run_batch_size],
                    db_path,
                    out_dir,
                    ref_seqs,
                    min_tail,
                    require_m,
                    magicblast_cmd,
                    num_threads,
                ): i
                for i in batch_indices
            }
            for future in as_completed(futures):
                for run_id, protein, crispr_str in future.result():
                    key = protein[:200] + "|" + protein[-200:]
                    if key in seen_seqs:
                        continue
                    seen_seqs.add(key)
                    seq_id = f"Cas13_{run_id}_{len(discoveries)}"
                    discoveries.append((seq_id, protein))
                    metadata.append((seq_id, run_id, crispr_str, "0"))

    if ref_nt_path.startswith(tempfile.gettempdir()):
        try:
            Path(ref_nt_path).unlink(missing_ok=True)
        except OSError:
            pass
    return discoveries, metadata


def save_discoveries(
    discoveries: List[Tuple[str, str]],
    metadata: List[Tuple[str, str, str, str]],
    output_dir: str = "data/raw_sequences",
) -> Tuple[str, str]:
    """Write deep_hits_YYYYMMDD_HHMMSS.fasta and matching _metadata.csv. Returns (fasta_path, csv_path)."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    fasta_path = str(Path(output_dir) / f"deep_hits_{ts}.fasta")
    csv_path = str(Path(output_dir) / f"deep_hits_{ts}_metadata.csv")
    with open(fasta_path, "w") as f:
        for seq_id, seq in discoveries:
            f.write(f">{seq_id}\n{seq}\n")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        import csv as csv_module
        w = csv_module.writer(f)
        # repeat_domains = CRISPR domain sequences (pipe-separated if multiple) for crRNA binding
        w.writerow(["sequence_id", "sra_accession", "repeat_domains", "score"])
        for row in metadata:
            w.writerow(row)
    return fasta_path, csv_path
