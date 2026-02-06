"""
Run ColabFold batch on 2-3 HEPN filtered FASTA.
Requires: Docker with ColabFold image, or local colabfold_batch install.
"""
import argparse
import os
import subprocess
from pathlib import Path


def run_colabfold_batch(
    input_fasta: str,
    output_dir: str,
    use_docker: bool = True,
    num_models: int = 1,
    num_recycle: int = 3,
) -> None:
    """
    Run colabfold_batch. If use_docker, invokes:
      docker run ... ghcr.io/.../colabfold:latest colabfold_batch ...
    Otherwise assumes colabfold_batch is on PATH.
    """
    input_path = Path(input_fasta).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    if use_docker:
        # ColabFold Docker: mount input/output, run colabfold_batch
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{input_path.parent}:/input",
            "-v", f"{out_path}:/output",
            "-v", f"{out_path}:/tmp",  # temp dir for ColabFold
            "ghcr.io/sokrypton/colabfold:latest",
            "colabfold_batch",
            f"/input/{input_path.name}",
            "/output",
            "--num-models", str(num_models),
            "--num-recycle", str(num_recycle),
        ]
    else:
        cmd = [
            "colabfold_batch",
            str(input_path),
            str(out_path),
            "--num-models", str(num_models),
            "--num-recycle", str(num_recycle),
        ]

    subprocess.run(cmd, check=True)
    print(f"[SUCCESS] ColabFold output: {out_path}")


def list_pdb_paths(output_dir: str) -> dict:
    """
    Scan ColabFold output and return {seq_id: pdb_path}.
    ColabFold writes to output_dir/<job_name>/<job_name>_relaxed_rank_001*.pdb
    or <job_name>_unrelaxed_rank_001*.pdb
    """
    result = {}
    out_path = Path(output_dir)
    for sub in out_path.iterdir():
        if not sub.is_dir():
            continue
        for pdb in sub.glob("*_rank_001*.pdb"):
            # job name often matches first part of seq_id
            seq_id = sub.name.replace("_unrelaxed", "").replace("_relaxed", "").split("_rank")[0]
            result[seq_id] = str(pdb.resolve())
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run ColabFold batch on 2-3 HEPN FASTA."
    )
    parser.add_argument(
        "--input",
        default="data/structure_pipeline/input_2-3_hepn.fasta",
        help="Input FASTA",
    )
    parser.add_argument(
        "--output-dir",
        default="data/structure_pipeline/structures/colabfold",
        help="Output directory",
    )
    parser.add_argument(
        "--no-docker",
        action="store_true",
        help="Use local colabfold_batch instead of Docker",
    )
    parser.add_argument(
        "--num-models",
        type=int,
        default=1,
        help="Number of models to predict",
    )
    parser.add_argument(
        "--num-recycle",
        type=int,
        default=3,
        help="Number of recycle steps",
    )
    args = parser.parse_args()

    run_colabfold_batch(
        args.input,
        args.output_dir,
        use_docker=not args.no_docker,
        num_models=args.num_models,
        num_recycle=args.num_recycle,
    )


if __name__ == "__main__":
    main()
