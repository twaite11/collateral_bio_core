"""
Run InterProScan 5 via Docker on 2-3 HEPN filtered FASTA.
Parses TSV output into domain annotations JSON for structure dashboard.
"""
import argparse
import json
import os
import subprocess
from pathlib import Path
from typing import Optional

INTERPRO_IMAGE = "interpro/interproscan:5.76-107.0"


def run_interproscan(
    input_fasta: str,
    output_dir: str,
    temp_dir: str,
    ips_data_dir: Optional[str] = None,
    cpu: int = 4,
) -> str:
    """
    Run InterProScan via Docker. Returns path to output TSV.
    Requires: Docker, InterProScan image, and data dir (mount at /opt/interproscan/data).
    """
    input_path = Path(input_fasta).resolve()
    out_path = Path(output_dir).resolve()
    tmp_path = Path(temp_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    tmp_path.mkdir(parents=True, exist_ok=True)

    # Default data dir: interproscan-5.x/data next to script or in project
    if ips_data_dir is None:
        proj = Path(__file__).resolve().parent.parent
        ips_data_dir = proj / "data" / "interproscan_data"
    data_path = Path(ips_data_dir).resolve()

    cmd = [
        "docker", "run", "--rm",
        "-v", f"{data_path}:/opt/interproscan/data",
        "-v", f"{input_path.parent}:/input",
        "-v", f"{out_path}:/output",
        "-v", f"{tmp_path}:/temp",
        INTERPRO_IMAGE,
        "--input", f"/input/{input_path.name}",
        "--output-dir", "/output",
        "--tempdir", "/temp",
        "--formats", "tsv",
        "--cpu", str(cpu),
    ]
    subprocess.run(cmd, check=True)
    # InterProScan writes <basename>.tsv to output dir
    tsv_name = input_path.stem + ".tsv"
    return str(out_path / tsv_name)


def parse_tsv(tsv_path: str) -> dict:
    """
    Parse InterProScan TSV into {seq_id: [{domain, signature, start, end, evalue}, ...]}.
    TSV columns: 1=accession, 4=analysis, 5=signature, 6=description, 7=start, 8=stop, 9=score.
    """
    result = {}
    with open(tsv_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seq_id = parts[0]
            analysis = parts[3]
            signature = parts[4]
            description = parts[5] if len(parts) > 5 else ""
            try:
                start = int(parts[6])
                stop = int(parts[7])
            except (ValueError, IndexError):
                continue
            score = parts[8] if len(parts) > 8 else "-"
            domain = f"{signature}" + (f" ({description})" if description else "")
            if seq_id not in result:
                result[seq_id] = []
            result[seq_id].append({
                "domain": domain,
                "signature": signature,
                "start": start,
                "end": stop,
                "evalue": score,
            })
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Run InterProScan on 2-3 HEPN FASTA and emit domain JSON."
    )
    parser.add_argument(
        "--input",
        default="data/structure_pipeline/input_2-3_hepn.fasta",
        help="Input FASTA",
    )
    parser.add_argument(
        "--output-dir",
        default="data/structure_pipeline/interpro",
        help="InterProScan output directory",
    )
    parser.add_argument(
        "--tempdir",
        default="data/structure_pipeline/interpro_temp",
        help="Temporary directory for InterProScan",
    )
    parser.add_argument(
        "--ips-data",
        default=None,
        help="InterProScan data directory (required for Docker run)",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=4,
        help="CPU cores",
    )
    parser.add_argument(
        "--skip-docker",
        action="store_true",
        help="Skip Docker run; only parse existing TSV",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / (Path(args.input).stem + ".tsv")
    json_path = out_dir / "domains.json"

    if not args.skip_docker:
        if not os.path.exists(args.input):
            print(f"[!] Input not found: {args.input}")
            return
        print("[*] Running InterProScan via Docker...")
        run_interproscan(
            args.input,
            args.output_dir,
            args.tempdir,
            ips_data_dir=args.ips_data,
            cpu=args.cpu,
        )
        print(f"[SUCCESS] TSV: {tsv_path}")
    else:
        if not tsv_path.exists():
            print(f"[!] TSV not found: {tsv_path}. Run without --skip-docker first.")
            return

    domains = parse_tsv(str(tsv_path))
    with open(json_path, "w") as f:
        json.dump(domains, f, indent=2)
    print(f"[SUCCESS] Domains JSON: {json_path} ({len(domains)} sequences)")


if __name__ == "__main__":
    main()
