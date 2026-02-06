"""
Run OmegaFold on 2-3 HEPN filtered FASTA.
Uses PyTorch (no JAX); works on RunPod and most GPU environments.

OmegaFold is not on PyPI and only supports Python 3.8–3.10. On Python 3.11/3.12:
  git clone https://github.com/HeliXonProtein/OmegaFold.git
  cd OmegaFold && pip install torch biopython
  python run_omegafold.py --omegafold-repo /path/to/OmegaFold
"""
import argparse
import os
import subprocess
from pathlib import Path
from typing import Optional


def run_omegafold(
    input_fasta: str,
    output_dir: str,
    num_cycle: int = 3,
    subbatch_size: Optional[int] = None,
    device: Optional[str] = None,
    omegafold_repo: Optional[str] = None,
) -> None:
    """
    Run OmegaFold. Outputs one PDB per sequence in output_dir.

    If omegafold_repo is set (or OMEGAFOLD_REPO env), runs python main.py
    from that clone instead of the 'omegafold' CLI (needed for Python 3.12).
    """
    input_path = Path(input_fasta).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    repo = omegafold_repo or os.environ.get("OMEGAFOLD_REPO")
    if repo:
        main_py = Path(repo).resolve() / "main.py"
        if not main_py.exists():
            raise FileNotFoundError(f"OmegaFold main.py not found: {main_py}")
        cmd = ["python", str(main_py), str(input_path), str(out_path)]
    else:
        cmd = ["omegafold", str(input_path), str(out_path)]

    cmd.extend(["--num_cycle", str(num_cycle)])
    if subbatch_size is not None:
        cmd.extend(["--subbatch_size", str(subbatch_size)])
    if device:
        cmd.extend(["--device", device])

    subprocess.run(cmd, check=True, cwd=repo if repo else None)
    print(f"[SUCCESS] OmegaFold output: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run OmegaFold on 2-3 HEPN FASTA."
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
    args = parser.parse_args()

    run_omegafold(
        args.input,
        args.output_dir,
        num_cycle=args.num_cycle,
        subbatch_size=args.subbatch_size,
        device=args.device,
        omegafold_repo=args.omegafold_repo,
    )


if __name__ == "__main__":
    main()
