#!/usr/bin/env python3
"""
Inspect ARCHS4 sample metadata (source_name_ch1) to verify normal vs cancer
classification. Run from project root. Requires human_matrix.h5 in data/expression_data/.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modules.targeting.archs4_loader import _classify_tissue
import h5py


def main() -> None:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(root)
    h5_path = os.getenv("ARCHS4_H5_PATH", "data/expression_data/human_matrix.h5")
    if not os.path.isabs(h5_path):
        h5_path = os.path.join(root, h5_path)
    if not os.path.isfile(h5_path):
        print(f"[!] Not found: {h5_path}")
        print("    Set ARCHS4_H5_PATH or place human_matrix.h5 in data/expression_data/")
        return
    print(f"[*] Reading {h5_path} ...")
    with h5py.File(h5_path, "r") as f:
        if "meta/samples/source_name_ch1" not in f:
            print("[!] meta/samples/source_name_ch1 not found.")
            return
        sample_dset = f["meta"]["samples"]["source_name_ch1"]
        tissues = [t.decode("utf-8") for t in sample_dset[:]]
    uniq: dict[str, int] = {}
    for t in tissues:
        u = (t or "").strip()
        uniq[u] = uniq.get(u, 0) + 1
    normal, cancer, unknown = [], [], []
    for name, count in sorted(uniq.items(), key=lambda x: -x[1]):
        if not name:
            continue
        ctx = _classify_tissue(name)
        if ctx == "normal":
            normal.append((name, count))
        elif ctx == "cancer":
            cancer.append((name, count))
        else:
            unknown.append((name, count))
    print("\n--- NORMAL (sample source_name_ch1) ---")
    for name, count in normal[:40]:
        print(f"  {count:6d}  {name[:80]}")
    print(f"  ... {len(normal)} unique normal labels total")
    print("\n--- CANCER ---")
    for name, count in cancer[:40]:
        print(f"  {count:6d}  {name[:80]}")
    print(f"  ... {len(cancer)} unique cancer labels total")
    print("\n--- UNKNOWN ---")
    for name, count in unknown[:40]:
        print(f"  {count:6d}  {name[:80]}")
    print(f"  ... {len(unknown)} unique unknown labels total")
    print("\n[*] Use these to tune NORMAL_KEYWORDS / CANCER_KEYWORDS in archs4_loader.py if needed.")


if __name__ == "__main__":
    main()
