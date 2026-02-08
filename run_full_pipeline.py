#!/usr/bin/env python3
"""
Full 4-step pipeline: Embed -> Mutate for drift -> Structure filter -> Identity filter.
Optionally run matchmaker on final FASTA.

Usage:
  python run_full_pipeline.py --input data/raw_sequences/deep_hits_latest.fasta
  python run_full_pipeline.py --input data/mined_sequences/family_grouped_20260205.fasta --skip-structure
  python run_full_pipeline.py --input data/design/drift_variants.fasta --skip-design --run-matchmaker

Expects:
  - data/references/known_cas13.fasta (reference Cas13 sequences for identity/drift)
  - OMEGAFOLD_REPO or --omegafold-repo for structure step (or --skip-structure)
  - data/structure_pipeline/references/ with 5W1H.pdb, 6DTD.pdb, 6IV9.pdb (or run visualization/run_tmscore.py once to download)
"""
import os
import sys
import argparse
import glob
from pathlib import Path

# Project root
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def _latest_fasta(pattern: str) -> str:
    """Return path to latest FASTA matching pattern."""
    files = sorted(glob.glob(str(ROOT / pattern)), key=os.path.getmtime, reverse=True)
    return files[0] if files else ""


def main():
    parser = argparse.ArgumentParser(description="Run full pipeline: embed -> mutate -> structure -> identity [-> matchmaker]")
    parser.add_argument("--input", default="", help="Input FASTA. If empty, use latest data/raw_sequences/deep_hits_*.fasta")
    parser.add_argument("--skip-design", action="store_true", help="Skip embed + mutate; use input as pre-designed pool")
    parser.add_argument("--skip-structure", action="store_true", help="Skip OmegaFold + bi-lobed filter (e.g. no GPU)")
    parser.add_argument("--skip-identity", action="store_true", help="Skip identity filter")
    parser.add_argument("--run-matchmaker", action="store_true", help="Run matchmaker on final FASTA")
    parser.add_argument("--omegafold-repo", default=os.environ.get("OMEGAFOLD_REPO"), help="Path to OmegaFold repo")
    parser.add_argument("--max-identity", type=float, default=0.85, help="Drift goal: max identity to known Cas13")
    parser.add_argument("--tm-threshold", type=float, default=0.4, help="Min TM-score for bi-lobed pass")
    parser.add_argument("--use-trans-cleavage-prompt", action="store_true", help="Use Gemini to suggest mutations that may increase trans-cleavage (maintain stability vs RfxCas13d/PspCas13a)")
    args = parser.parse_args()

    input_fasta = args.input
    if not input_fasta:
        input_fasta = _latest_fasta("data/raw_sequences/deep_hits_*.fasta")
    if not input_fasta:
        input_fasta = _latest_fasta("data/mined_sequences/deep_hits_*.fasta")
    if not input_fasta or not Path(input_fasta).exists():
        print("[!] No input FASTA. Provide --input or place deep_hits_*.fasta in data/raw_sequences/ or data/mined_sequences/")
        return 1

    print(f"[*] Input FASTA: {input_fasta}")

    current_fasta = input_fasta

    # Step 2a: Embed pool (optional; for interpolation or logging)
    if not args.skip_design:
        from modules.design.embed_pool import embed_fasta
        embed_dir = str(ROOT / "data" / "design" / "embeddings")
        embed_fasta(current_fasta, embed_dir)
        # Step 2b: Mutate for drift
        from modules.design.mutate_for_drift import main as mutate_main
        mutate_argv = [
            "mutate_for_drift.py",
            "--input", current_fasta,
            "--output", str(ROOT / "data" / "design" / "drift_variants.fasta"),
            "--max-identity", str(args.max_identity),
        ]
        if os.environ.get("ESM_REFERENCE_FASTA"):
            mutate_argv.extend(["--stability-refs", os.environ.get("ESM_REFERENCE_FASTA")])
        if getattr(args, "use_trans_cleavage_prompt", False) or os.environ.get("USE_TRANS_CLEAVAGE_PROMPT", "").lower() in ("1", "true", "yes"):
            mutate_argv.append("--use-trans-cleavage-prompt")
        sys.argv = mutate_argv
        mutate_main()
        current_fasta = str(ROOT / "data" / "design" / "drift_variants.fasta")
        if not Path(current_fasta).exists():
            print("[!] mutate_for_drift produced no output.")
            return 1

    # Step 3: Structure filter
    if not args.skip_structure:
        if not args.omegafold_repo and not os.environ.get("OMEGAFOLD_REPO"):
            print("[!] OMEGAFOLD_REPO not set. Use --skip-structure or set OMEGAFOLD_REPO.")
            return 1
        from modules.structure_filter.run_filter import run_full_filter
        n = run_full_filter(
            current_fasta,
            str(ROOT / "data" / "structure_pipeline" / "passed_structures.fasta"),
            str(ROOT / "data" / "structure_pipeline" / "failed_structures.log"),
            str(ROOT / "data" / "structure_pipeline" / "structures" / "omegafold"),
            str(ROOT / "data" / "structure_pipeline" / "references"),
            tm_threshold=args.tm_threshold,
            omegafold_repo=args.omegafold_repo,
        )
        current_fasta = str(ROOT / "data" / "structure_pipeline" / "passed_structures.fasta")
        if n == 0 and Path(current_fasta).exists() and os.path.getsize(current_fasta) == 0:
            print("[!] No sequences passed structure filter.")
            return 1
    else:
        # If skipping structure, copy current FASTA to passed_structures for identity step
        pass_fasta = str(ROOT / "data" / "structure_pipeline" / "passed_structures.fasta")
        import shutil
        Path(pass_fasta).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(current_fasta, pass_fasta)
        current_fasta = pass_fasta

    # Step 4: Identity filter
    if not args.skip_identity:
        from modules.identity_filter import run_identity_filter
        out_passed = str(ROOT / "data" / "identity_filtered" / "passed.fasta")
        out_meta = str(ROOT / "data" / "identity_filtered" / "identity_metadata.csv")
        n =         run_identity_filter(
            current_fasta,
            str(ROOT / "data" / "references" / "known_cas13.fasta"),
            out_passed,
            out_meta,
            max_identity=args.max_identity,
        )
        current_fasta = out_passed
        if n == 0:
            print("[!] No sequences passed identity filter.")
            return 1

    # Optional: Matchmaker
    if args.run_matchmaker:
        from modules.matchmaker import MasterMatchmaker
        target_file = str(ROOT / "data" / "targets" / "high_specificity_targets.csv")
        if not Path(target_file).exists():
            target_file = str(ROOT / "data" / "targets" / "known_fusions.csv")
        mm = MasterMatchmaker(enzyme_fasta=current_fasta, target_file=target_file)
        mm.run_matching()
        mm.save_leads("lead_candidates.csv")
        print("[+] Matchmaker wrote lead_candidates.csv")

    print(f"[SUCCESS] Pipeline complete. Final FASTA: {current_fasta}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
