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
    Automatically detects and uses OmegaFold venv Python if available.
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
        
        cmd = [python_exe, str(main_py), str(input_path), str(out_path)]
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
