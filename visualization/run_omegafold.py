"""
Run OmegaFold on 2-3 HEPN filtered FASTA.
Uses PyTorch (no JAX); works on RunPod and most GPU environments.

Each sequence is predicted in its **own subprocess** so GPU VRAM is fully
released between predictions.  This prevents the CUDA OOM cascade where one
failure fragments memory and causes all subsequent predictions to fail.

OmegaFold is not on PyPI and only supports Python 3.8–3.10. On Python 3.11/3.12:
  git clone https://github.com/HeliXonProtein/OmegaFold.git
  cd OmegaFold && pip install torch biopython
  python run_omegafold.py --omegafold-repo /path/to/OmegaFold
"""
import argparse
import logging
import os
import subprocess
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


def _pdb_exists(record, output_dir: Path) -> bool:
    """Check whether a PDB file already exists for this sequence (resume support)."""
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


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

def run_omegafold(
    input_fasta: str,
    output_dir: str,
    num_cycle: int = 3,
    subbatch_size: Optional[int] = None,
    device: Optional[str] = None,
    omegafold_repo: Optional[str] = None,
    max_residues: int = 0,
) -> None:
    """
    Run OmegaFold.  Outputs one PDB per sequence in output_dir.

    Each sequence is run in its own subprocess so CUDA memory is fully
    released between predictions (prevents OOM cascades).

    Args:
        max_residues: Skip sequences longer than this (0 = no limit).
            750 is a reasonable cap for a 24 GB GPU.
    """
    input_path = Path(input_fasta).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    repo = omegafold_repo or os.environ.get("OMEGAFOLD_REPO")
    if repo:
        repo_path = Path(repo).resolve()
        main_py = repo_path / "main.py"
        if not main_py.exists():
            raise FileNotFoundError(f"OmegaFold main.py not found: {main_py}")
        
        # Try to find and use OmegaFold venv Python
        venv_python = _find_omegafold_venv_python(repo_path)
        if venv_python:
            print(f"[*] Using OmegaFold venv Python: {venv_python}")
            python_exe = venv_python
        else:
            python_exe = "python"
            print(f"[*] Using system Python (OmegaFold venv not found, using: {python_exe})")
    else:
        repo_path = None
        main_py = None
        python_exe = None

    # --- Read and sort sequences (shortest first) ---
    records = list(SeqIO.parse(input_path, "fasta"))
    records.sort(key=lambda r: len(r.seq))
    total_input = len(records)

    if total_input == 0:
        print("[!] No sequences in input FASTA.")
        return

    # --- Filter: skip already-predicted and too-long sequences ---
    skipped_existing = 0
    skipped_too_long = 0
    to_predict = []
    for rec in records:
        if _pdb_exists(rec, out_path):
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

    print(f"[*] OmegaFold: predicting {total} sequences one-at-a-time (subprocess per sequence).")
    print(f"[*] PYTORCH_CUDA_ALLOC_CONF will include expandable_segments:True to reduce fragmentation.")

    succeeded = 0
    failed = 0
    consecutive_failures = 0
    MAX_CONSECUTIVE_FAILURES = 5  # stop early after this many failures in a row

    env = _gpu_env()

    for idx, rec in enumerate(to_predict, 1):
        n_residues = len(rec.seq)
        print(f"[*] ({idx}/{total}) Predicting {rec.id} ({n_residues} residues)...")

        # Write single-sequence temp FASTA
        fd, temp_path = tempfile.mkstemp(suffix=".fasta", prefix="omegafold_single_")
        try:
            os.close(fd)
            SeqIO.write([rec], temp_path, "fasta")

            # Build command
            if repo and main_py is not None:
                cmd = [python_exe, str(main_py), temp_path, str(out_path)]
            else:
                cmd = ["omegafold", temp_path, str(out_path)]

            cmd.extend(["--num_cycle", str(num_cycle)])
            if subbatch_size is not None:
                cmd.extend(["--subbatch_size", str(subbatch_size)])
            if device:
                cmd.extend(["--device", device])

            result = subprocess.run(
                cmd,
                cwd=repo if repo else None,
                env=env,
            )

            if result.returncode == 0 and _pdb_exists(rec, out_path):
                succeeded += 1
                consecutive_failures = 0
                print(f"    [OK] {rec.id}")
            else:
                failed += 1
                consecutive_failures += 1
                print(f"    [FAIL] {rec.id} (exit code {result.returncode}, "
                      f"consecutive failures: {consecutive_failures}/{MAX_CONSECUTIVE_FAILURES})")

        except Exception as exc:
            failed += 1
            consecutive_failures += 1
            print(f"    [ERROR] {rec.id}: {exc}")

        finally:
            try:
                os.unlink(temp_path)
            except OSError:
                pass

        # Early termination: sequences are sorted by length, so if several
        # consecutive predictions fail (likely OOM), longer ones will also fail.
        if consecutive_failures >= MAX_CONSECUTIVE_FAILURES:
            remaining = total - idx
            print(f"[!!] {MAX_CONSECUTIVE_FAILURES} consecutive failures -- likely persistent OOM at "
                  f"{n_residues}+ residues. Stopping early ({remaining} sequences remaining). "
                  f"Consider --max-residues {n_residues - 50} or a larger GPU.")
            break

    print(f"[DONE] OmegaFold: {succeeded} succeeded, {failed} failed, "
          f"{skipped_existing} skipped (existing), {skipped_too_long} skipped (too long).")
    print(f"[SUCCESS] OmegaFold output: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run OmegaFold on 2-3 HEPN FASTA (one subprocess per sequence for VRAM safety)."
    )
    parser.add_argument(
        "--input",
        default="data/structure_pipeline/input_2-3_hepn.fasta",
        help="Input FASTA",
    )
    parser.add_argument(
        "--output-dir",
        default="data/structure_pipeline/structures/omegafold",
        help="Output directory",
    )
    parser.add_argument(
        "--num-cycle",
        type=int,
        default=3,
        help="Number of optimization cycles (default: 3)",
    )
    parser.add_argument(
        "--subbatch-size",
        type=int,
        default=None,
        help="Subbatch size to reduce GPU memory (lower = less VRAM, slower)",
    )
    parser.add_argument(
        "--device",
        default=None,
        help="Device (cuda, cpu, or cuda:0)",
    )
    parser.add_argument(
        "--omegafold-repo",
        default=None,
        help="Path to cloned OmegaFold repo (use when omegafold pip install fails, e.g. Python 3.12). Or set OMEGAFOLD_REPO.",
    )
    parser.add_argument(
        "--max-residues",
        type=int,
        default=0,
        help="Skip sequences longer than this (0 = no limit; 750 is reasonable for 24 GB GPU)",
    )
    args = parser.parse_args()

    run_omegafold(
        args.input,
        args.output_dir,
        num_cycle=args.num_cycle,
        subbatch_size=args.subbatch_size,
        device=args.device,
        omegafold_repo=args.omegafold_repo,
        max_residues=args.max_residues,
    )


if __name__ == "__main__":
    main()
