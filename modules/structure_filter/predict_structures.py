"""
Run OmegaFold (or ESMFold) on a FASTA; output one PDB per sequence.
Used by the pipeline for in-silico structure filtering.

Key design: each sequence is predicted in its **own subprocess** so that GPU
memory is fully released between predictions.  This prevents the CUDA OOM
cascade where one failure fragments VRAM and causes all subsequent predictions
to also OOM.
"""
import os
import logging
import subprocess
import argparse
import tempfile
from pathlib import Path
from typing import Optional

from Bio import SeqIO


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_omegafold_venv_python(repo_path: Path) -> Optional[str]:
    """
    Find OmegaFold venv Python executable.
    Checks common venv locations relative to the OmegaFold repo.
    Returns path to Python executable if found, None otherwise.
    """
    # Common venv names to check
    venv_names = [
        "venv_omegafold1",
        "venv_omegafold",
        "omegafold_venv",
        "venv",
        "env",
    ]
    
    # Check parent directory (e.g., /workspace/venv_omegafold1 when repo is /workspace/OmegaFold)
    parent_dir = repo_path.parent
    for venv_name in venv_names:
        venv_path = parent_dir / venv_name
        if venv_path.exists() and venv_path.is_dir():
            # Check for Python executable (Linux/Mac)
            python_exe = venv_path / "bin" / "python"
            if python_exe.exists():
                return str(python_exe)
            # Check for Python executable (Windows)
            python_exe = venv_path / "Scripts" / "python.exe"
            if python_exe.exists():
                return str(python_exe)
    
    # Also check if OMEGAFOLD_VENV is explicitly set
    omegafold_venv = os.environ.get("OMEGAFOLD_VENV")
    if omegafold_venv:
        venv_path = Path(omegafold_venv).resolve()
        python_exe = venv_path / "bin" / "python"
        if python_exe.exists():
            return str(python_exe)
        python_exe = venv_path / "Scripts" / "python.exe"
        if python_exe.exists():
            return str(python_exe)
    
    return None


def _pdb_exists_for_record(record, output_dir: Path) -> bool:
    """Check whether a PDB file already exists for this sequence (resume support)."""
    # OmegaFold names output PDBs after the FASTA header (record.id)
    pdb_path = output_dir / f"{record.id}.pdb"
    return pdb_path.exists() and pdb_path.stat().st_size > 0


def _gpu_env() -> dict:
    """Return an env dict with PYTORCH_CUDA_ALLOC_CONF set to reduce fragmentation."""
    env = os.environ.copy()
    alloc_conf = env.get("PYTORCH_CUDA_ALLOC_CONF", "")
    if "expandable_segments" not in alloc_conf:
        extra = "expandable_segments:True"
        env["PYTORCH_CUDA_ALLOC_CONF"] = f"{alloc_conf},{extra}" if alloc_conf else extra
    return env


def _run_omegafold_single_fasta(
    input_fasta: str,
    output_dir: str,
    num_cycle: int,
    device: Optional[str],
    omegafold_repo: Optional[str],
    python_exe: str,
    main_py: Path,
    repo: Optional[str],
) -> None:
    """Run OmegaFold once on a single FASTA file (one subprocess)."""
    cmd = [python_exe, str(main_py), input_fasta, output_dir]
    cmd.extend(["--num_cycle", str(num_cycle)])
    if device:
        cmd.extend(["--device", device])
    subprocess.run(cmd, check=True, cwd=repo if repo else None, env=_gpu_env())


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

def run_omegafold(
    input_fasta: str,
    output_dir: str,
    num_cycle: int = 3,
    device: Optional[str] = None,
    omegafold_repo: Optional[str] = None,
    batch_size: int = 1,
    max_residues: int = 0,
) -> None:
    """Run OmegaFold on input FASTA.  Writes PDBs to output_dir.

    Each batch is run as its own subprocess so GPU VRAM is fully released
    between batches.  Default batch_size=1 (safest for VRAM).

    Args:
        max_residues: Skip sequences longer than this many amino acids
            (0 = no limit).  750 is a reasonable cap for a 24 GB GPU.
    """
    input_path = Path(input_fasta).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    repo = omegafold_repo or os.environ.get("OMEGAFOLD_REPO")
    if repo:
        repo_path = Path(repo).resolve()
        # Check if repo path is valid (not a placeholder)
        if not repo_path.exists() or "path/to" in str(repo_path).lower() or repo_path == Path("/path/to/OmegaFold"):
            raise FileNotFoundError(
                f"Invalid OmegaFold repository path: {repo}\n"
                f"Please set OMEGAFOLD_REPO environment variable to the actual OmegaFold repository path, "
                f"or use --omegafold-repo argument, or use --skip-structure to skip structure prediction."
            )
        main_py = repo_path / "main.py"
        if not main_py.exists():
            raise FileNotFoundError(
                f"OmegaFold main.py not found at: {main_py}\n"
                f"Repository path: {repo_path}\n"
                f"Please verify that OmegaFold is installed at this location."
            )
        venv_python = _find_omegafold_venv_python(repo_path)
        python_exe = venv_python if venv_python else "python"
        if venv_python:
            print(f"[*] Using OmegaFold venv Python: {venv_python}")
        else:
            print(f"[*] Using system Python (OmegaFold venv not found, using: {python_exe})")
    else:
        python_exe = "python"
        main_py = None
        repo_path = None

    records = list(SeqIO.parse(input_path, "fasta"))
    total_input = len(records)
    if total_input == 0:
        print("[!] No sequences in input FASTA.")
        return

    # --- Sort by sequence length (shortest first) so cheaper predictions run
    #     first and we discover the VRAM ceiling gradually. ---
    records.sort(key=lambda r: len(r.seq))

    # --- Skip sequences whose PDB already exists (resume support) ---
    skipped_existing = 0
    skipped_too_long = 0
    to_predict = []
    for rec in records:
        if _pdb_exists_for_record(rec, out_path):
            skipped_existing += 1
            continue
        if max_residues > 0 and len(rec.seq) > max_residues:
            skipped_too_long += 1
            logging.info(f"Skipping {rec.id}: {len(rec.seq)} residues exceeds max_residues={max_residues}")
            continue
        to_predict.append(rec)

    if skipped_existing:
        print(f"[*] Skipping {skipped_existing} sequence(s) with existing PDB files (resume).")
    if skipped_too_long:
        print(f"[*] Skipping {skipped_too_long} sequence(s) exceeding {max_residues} residues.")

    total = len(to_predict)
    if total == 0:
        print("[*] All sequences already predicted or skipped. Nothing to do.")
        return

    if batch_size < 1:
        batch_size = 1
    n_batches = (total + batch_size - 1) // batch_size
    print(f"[*] OmegaFold: {total} sequences in {n_batches} batch(es) of up to {batch_size}.")
    print(f"[*] PYTORCH_CUDA_ALLOC_CONF will include expandable_segments:True to reduce fragmentation.")

    consecutive_oom_failures = 0
    MAX_CONSECUTIVE_OOM = 5  # stop early if this many batches fail in a row

    for i in range(0, total, batch_size):
        chunk = to_predict[i : i + batch_size]
        batch_num = i // batch_size + 1
        seq_lengths = [len(r.seq) for r in chunk]
        fd, temp_path = tempfile.mkstemp(suffix=".fasta", prefix="omegafold_batch_")
        try:
            os.close(fd)
            SeqIO.write(chunk, temp_path, "fasta")
            ids_str = ", ".join(r.id for r in chunk)
            print(f"[*] Batch {batch_num}/{n_batches} ({len(chunk)} seq, "
                  f"{min(seq_lengths)}-{max(seq_lengths)} residues): {ids_str}")
            if repo and main_py is not None:
                _run_omegafold_single_fasta(
                    temp_path,
                    str(out_path),
                    num_cycle,
                    device,
                    omegafold_repo,
                    python_exe,
                    main_py,
                    repo,
                )
            else:
                subprocess.run(
                    ["omegafold", temp_path, str(out_path), "--num_cycle", str(num_cycle)]
                    + (["--device", device] if device else []),
                    check=True,
                    env=_gpu_env(),
                )
            consecutive_oom_failures = 0  # reset on success
        except subprocess.CalledProcessError as exc:
            # Check whether any PDB was actually produced (partial success)
            produced = sum(1 for r in chunk if _pdb_exists_for_record(r, out_path))
            if produced:
                print(f"[!] Batch {batch_num} partially succeeded ({produced}/{len(chunk)} PDBs produced).")
                consecutive_oom_failures = 0
            else:
                consecutive_oom_failures += 1
                print(f"[!] Batch {batch_num} failed (exit code {exc.returncode}). "
                      f"Consecutive failures: {consecutive_oom_failures}/{MAX_CONSECUTIVE_OOM}")
            if consecutive_oom_failures >= MAX_CONSECUTIVE_OOM:
                print(f"[!!] {MAX_CONSECUTIVE_OOM} consecutive batch failures -- likely persistent OOM. "
                      f"Stopping early. Remaining {total - i - batch_size} sequences skipped. "
                      f"Consider --max-residues or a larger GPU.")
                break
        finally:
            try:
                os.unlink(temp_path)
            except OSError:
                pass

    produced_total = sum(1 for r in to_predict if _pdb_exists_for_record(r, out_path))
    print(f"[+] Structures written to {out_path}  ({produced_total}/{total_input} sequences have PDBs)")


def main():
    parser = argparse.ArgumentParser(description="Predict structures with OmegaFold for pipeline filter.")
    parser.add_argument("--input", default="data/design/drift_variants.fasta", help="Input FASTA")
    parser.add_argument("--output-dir", default="data/structure_pipeline/structures/omegafold", help="Output PDB dir")
    parser.add_argument("--num-cycle", type=int, default=3)
    parser.add_argument("--device", default=None, help="cuda or cpu")
    parser.add_argument("--omegafold-repo", default=None, help="Path to OmegaFold repo (or OMEGAFOLD_REPO)")
    parser.add_argument("--batch-size", type=int, default=1,
                        help="Sequences per subprocess batch (default 1 = safest for VRAM)")
    parser.add_argument("--max-residues", type=int, default=0,
                        help="Skip sequences longer than this (0 = no limit; 750 is reasonable for 24 GB GPU)")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"[!] Input not found: {args.input}")
        return 1
    run_omegafold(
        args.input,
        args.output_dir,
        num_cycle=args.num_cycle,
        device=args.device,
        omegafold_repo=args.omegafold_repo,
        batch_size=args.batch_size,
        max_residues=args.max_residues,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
