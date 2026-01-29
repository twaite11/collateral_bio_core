from __future__ import annotations

from typing import Any, List, Optional, Tuple

import h5py
import numpy as np
import pandas as pd

# Keywords to classify ARCHS4 samples as normal vs cancer (source_name_ch1)
# Broadened to match common GEO/ARCHS4 annotations
NORMAL_KEYWORDS = (
    "normal", "healthy", "non-tumor", "non-tumour", "adjacent normal", "nat ",
    "nat,", "control", "wild type", "non-malignant", "benign ", "non-cancer",
    "non-diseased", "healthy tissue", "uninvolved", "matched normal",
)
CANCER_KEYWORDS = (
    "tumor", "tumour", "cancer", "carcinoma", "adenocarcinoma", "squamous",
    "melanoma", "lymphoma", "leukemia", "sarcoma", "glioblastoma", "gbm",
    "metastasis", "metastatic", "malignant", "neoplasm", "tcga ",
)


def _classify_tissue(tissue: str) -> str:
    """Return 'normal', 'cancer', or 'unknown'."""
    t = (tissue or "").lower()
    is_normal = any(k in t for k in NORMAL_KEYWORDS)
    is_cancer = any(k in t for k in CANCER_KEYWORDS)
    if is_cancer and not is_normal:
        return "cancer"
    if is_normal and not is_cancer:
        return "normal"
    if is_normal and is_cancer:
        return "normal"  # e.g. "normal adjacent to tumor" – treat as normal for safety
    return "unknown"


def _tissue_matches_organ(tissue: str, organ_keywords: List[str]) -> bool:
    """True if tissue string matches any organ keyword."""
    t = (tissue or "").lower()
    return any(kw.lower() in t for kw in organ_keywords)


class ARCHS4Loader:
    def __init__(self, h5_path="data/expression_data/human_matrix.h5"):
        self.h5_path = h5_path

    def get_gene_expression(self, gene_symbol):
        """
        Extracts expression data for a single gene, auto-detecting matrix orientation.
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                # 1. Locate Gene Symbols
                # Try 'meta/genes/gene_symbol' (v2) first, then 'meta/genes/genes' (v1)
                if 'meta/genes/gene_symbol' in f:
                    gene_dset = f['meta']['genes']['gene_symbol']
                else:
                    gene_dset = f['meta']['genes']['genes']
                
                # Decode all genes to strings for searching
                all_genes = [g.decode('utf-8') for g in gene_dset[:]]
                
                if gene_symbol not in all_genes:
                    print(f"[!] Gene {gene_symbol} not found in database.")
                    return None
                
                # Get the index of the gene
                gene_idx = all_genes.index(gene_symbol)
                
                # 2. Get Sample Metadata (Tissues)
                # Try 'meta/samples/source_name_ch1'
                if 'meta/samples/source_name_ch1' in f:
                    sample_dset = f['meta']['samples']['source_name_ch1']
                else:
                    print("[!] Could not find sample source metadata.")
                    return None
                    
                tissues = [t.decode('utf-8') for t in sample_dset[:]]
                
                # 3. Smart Extraction (Check Shape)
                expression_dset = f['data']['expression']
                shape = expression_dset.shape
                
                # Case A: Matrix is (Genes, Samples) -> Standard for Bio
                if shape[0] == len(all_genes):
                    # print(f"[*] Detected (Genes, Samples) orientation: {shape}")
                    expression_values = expression_dset[gene_idx, :]
                    
                # Case B: Matrix is (Samples, Genes) -> Compressed format
                elif shape[1] == len(all_genes):
                    # print(f"[*] Detected (Samples, Genes) orientation: {shape}")
                    expression_values = expression_dset[:, gene_idx]
                    
                else:
                    print(f"[!] Shape mismatch! Matrix: {shape}, Genes: {len(all_genes)}, Samples: {len(tissues)}")
                    return None

                # 4. Return Data
                return pd.DataFrame({
                    'Tissue': tissues,
                    'Expression': expression_values
                })
                
        except FileNotFoundError:
            print(f"[!] Database not found at {self.h5_path}.")
            return None
        except Exception as e:
            print(f"[!] Error reading H5 file: {e}")
            return None

    def get_gene_expression_normal_vs_cancer(self, gene_symbol: str, organ_keywords: Optional[List[str]] = None):
        """
        Returns expression split by context: {'normal': df, 'cancer': df}.
        If organ_keywords is given, restrict to samples whose Tissue matches any keyword (organ-specific).
        """
        df = self.get_gene_expression(gene_symbol)
        if df is None or df.empty:
            return {"normal": pd.DataFrame(), "cancer": pd.DataFrame()}
        df = df.copy()
        if organ_keywords:
            mask = df["Tissue"].astype(str).apply(lambda t: _tissue_matches_organ(t, organ_keywords))
            df = df.loc[mask].reset_index(drop=True)
            if df.empty:
                return {"normal": pd.DataFrame(), "cancer": pd.DataFrame()}
        df["Context"] = df["Tissue"].astype(str).map(_classify_tissue)
        out = {
            "normal": df.loc[df["Context"] == "normal", ["Tissue", "Expression"]].reset_index(drop=True),
            "cancer": df.loc[df["Context"] == "cancer", ["Tissue", "Expression"]].reset_index(drop=True),
        }
        return out

    def fusion_absent_in_normal_present_in_cancer(
        self,
        fusion_name: str,
        normal_max_tpm: float = 1.0,
        cancer_min_tpm: float = 1.0,
        fusion_cancers: Optional[List[str]] = None,
        enrichment_factor: float = 2.0,
        use_organ_specific: bool = True,
    ) -> Tuple[bool, str, List[str]]:
        """
        Check parent genes: (1) low/absent in normal, (2) present in cancer.
        When fusion_cancers (TCGA codes) are provided and use_organ_specific:
        - Restrict to organ(s) for those cancers; require enrichment (cancer > enrichment_factor * normal)
          instead of strict "absent in normal", so cancer-specific signal can pass.
        Otherwise: global normal/cancer, strict normal_max threshold.
        Returns (ok, reason, details).
        """
        try:
            from modules.targeting.fusion_metadata import organs_for_cancers
        except Exception:
            def organs_for_cancers(_: List[str]) -> List[str]:
                return []

        parts = [p.strip() for p in fusion_name.replace("--", "-").split("-") if len(p.strip()) > 2]
        if not parts:
            return False, "No valid genes", []
        details: List[str] = []
        organ_kw: List[str] = []
        if use_organ_specific and fusion_cancers:
            organ_kw = organs_for_cancers(fusion_cancers)
        if organ_kw:
            details.append(f"Organ context: {organ_kw} (cancers: {fusion_cancers})")

        for gene in parts:
            ctx = self.get_gene_expression_normal_vs_cancer(gene, organ_keywords=organ_kw if organ_kw else None)
            n_df, c_df = ctx["normal"], ctx["cancer"]

            n_mean = float(n_df["Expression"].mean(skipna=True)) if not n_df.empty else 0.0
            n_max = float(n_df["Expression"].max(skipna=True)) if not n_df.empty else 0.0
            c_mean = float(c_df["Expression"].mean(skipna=True)) if not c_df.empty else 0.0
            if np.isnan(n_mean):
                n_mean = 0.0
            if np.isnan(n_max):
                n_max = 0.0
            if np.isnan(c_mean):
                c_mean = 0.0

            n_samples = len(n_df)
            c_samples = len(c_df)

            if organ_kw and n_samples > 0 and c_samples > 0:
                # Organ-specific: require enrichment (cancer > factor * normal)
                enriched = c_mean >= enrichment_factor * n_mean
                above_min = c_mean >= cancer_min_tpm
                low_normal = n_max <= max(normal_max_tpm, 50.0)
                line = (
                    f"{gene}: organ n={n_samples} c={c_samples} "
                    f"normal_mean={n_mean:.1f} normal_max={n_max:.1f} cancer_mean={c_mean:.1f} "
                    f"enriched={enriched} above_min={above_min} low_normal={low_normal}"
                )
                details.append(line)
                if not above_min:
                    return False, f"{gene} cancer mean {c_mean:.1f} < {cancer_min_tpm} (organ-specific)", details
                if not enriched and n_mean > 0:
                    return False, f"{gene} not enriched in cancer vs normal (organ) {c_mean:.1f} vs {n_mean:.1f}", details
                continue

            if not organ_kw or (n_samples == 0 and c_samples == 0):
                # Global or no organ data: strict absent-in-normal + present-in-cancer
                if n_df.empty:
                    details.append(f"{gene}: no normal-tissue data (global)")
                    return False, f"{gene}: no normal-tissue data (cannot verify)", details
                n_ok = n_max <= normal_max_tpm
                c_ok = c_df.empty or c_mean >= cancer_min_tpm
                details.append(
                    f"{gene}: normal_max={n_max:.1f} (ok={n_ok}), cancer_mean={c_mean:.1f} (ok={c_ok})"
                )
                if not n_ok:
                    return False, f"{gene} expressed in normal (max {n_max:.1f} > {normal_max_tpm})", details
                if not c_df.empty and not c_ok:
                    return False, f"{gene} not expressed in cancer (mean {c_mean:.1f} < {cancer_min_tpm})", details
                continue

            if n_samples == 0 and c_samples > 0:
                # Organ-specific, no normal in organ: require cancer above min only
                details.append(f"{gene}: no normal in organ, cancer_mean={c_mean:.1f}")
                if c_mean < cancer_min_tpm:
                    return False, f"{gene} cancer mean {c_mean:.1f} < {cancer_min_tpm} (organ, no normal)", details
                continue

            if n_samples > 0 and c_samples == 0:
                details.append(f"{gene}: organ normal n={n_samples} but no cancer samples")
                return False, f"{gene}: no cancer samples in organ", details

        return True, "Pass (organ-specific enriched or global absent-in-normal)", details