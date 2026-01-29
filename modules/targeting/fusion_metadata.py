"""
Load fusion → associated cancers from KB_and_Pub and novel_Recur CSVs.
Used for context-specific ARCHS4 filtering (organ-specific normal/cancer).
"""
from __future__ import annotations

import os

import pandas as pd

# TCGA cancer code -> organ/tissue keywords for ARCHS4 source_name_ch1 matching
TCGA_TO_ORGAN: dict[str, list[str]] = {
    "ACC": ["adrenal"],
    "BLCA": ["bladder", "urothelial"],
    "BRCA": ["breast"],
    "CESC": ["cervix", "cervical"],
    "CHOL": ["bile duct", "cholangiocarcinoma", "liver"],
    "COAD": ["colon", "colorectal"],
    "DLBC": ["lymph", "lymphoma", "diffuse large b-cell"],
    "ESCA": ["esophagus", "esophageal"],
    "GBM": ["brain", "glioblastoma", "gbm", "cerebral"],
    "HNSC": ["head", "neck", "oral", "pharynx", "larynx"],
    "KICH": ["kidney"],
    "KIRC": ["kidney", "renal"],
    "KIRP": ["kidney", "renal"],
    "LAML": ["blood", "bone marrow", "leukemia", "acute myeloid"],
    "LGG": ["brain", "glioma", "astrocytoma", "oligodendroglioma"],
    "LIHC": ["liver", "hepatocellular"],
    "LUAD": ["lung", "adenocarcinoma"],
    "LUSC": ["lung", "squamous"],
    "MESO": ["mesothelioma", "pleura"],
    "OV": ["ovary", "ovarian"],
    "PAAD": ["pancreas", "pancreatic"],
    "PCPG": ["adrenal", "pheochromocytoma"],
    "PRAD": ["prostate"],
    "READ": ["rectum", "rectal", "colorectal"],
    "SARC": ["soft tissue", "sarcoma", "bone"],
    "SKCM": ["skin", "melanoma"],
    "STAD": ["stomach", "gastric"],
    "TGCT": ["testis", "testicular", "germ cell"],
    "THCA": ["thyroid"],
    "THYM": ["thymus"],
    "UCEC": ["uterus", "endometrium", "uterine"],
    "UCS": ["uterus", "uterine"],
    "UVM": ["uveal", "eye", "melanoma"],
}


def parse_associated_disease(disease_str: str) -> list[str]:
    """Parse 'PRAD, LGG' or 'Unknown' into list of TCGA codes."""
    if not disease_str or (isinstance(disease_str, str) and disease_str.strip().lower() == "unknown"):
        return []
    s = str(disease_str).strip()
    return [c.strip() for c in s.replace(";", ",").split(",") if c.strip()]


def load_fusion_to_cancers(
    kb_matrix_path: str = "data/KB_and_Pub_Recur_per_cancer.csv",
    novel_matrix_path: str | None = "data/novel_Recur_per_cancer.csv",
    data_dir: str | None = None,
) -> dict[str, list[str]]:
    """
    Build fusion → [TCGA codes] from matrix CSVs (rows=cancer, cols=fusion).
    Returns dict[str, list[str]].
    """
    out: dict[str, list[str]] = {}
    root = data_dir or os.getcwd()

    def add_from_path(path: str) -> None:
        p = os.path.join(root, path) if not os.path.isabs(path) else path
        if not os.path.isfile(p):
            return
        try:
            df = pd.read_csv(p)
            if df.empty or "Cancer" not in df.columns:
                return
            df = df.set_index("Cancer")
            known_tcga = set(TCGA_TO_ORGAN.keys())
            for fusion in df.columns:
                if not isinstance(fusion, str) or "Unnamed" in fusion:
                    continue
                try:
                    vals = pd.to_numeric(df[fusion], errors="coerce").fillna(0)
                except Exception:
                    continue
                cancers = [c for c in vals[vals > 0].index.astype(str).tolist() if c.strip() and c.upper() in known_tcga]
                if not cancers:
                    continue
                existing = out.get(fusion, [])
                for c in cancers:
                    if c not in existing:
                        existing.append(c)
                out[fusion] = existing
        except Exception as e:
            print(f"[!] fusion_metadata: error reading {p}: {e}")

    add_from_path(kb_matrix_path)
    if novel_matrix_path:
        add_from_path(novel_matrix_path)
    return out


def organs_for_cancers(tcga_codes: list[str]) -> list[str]:
    """Return unique organ keywords for given TCGA codes."""
    seen: set[str] = set()
    for c in tcga_codes:
        for kw in TCGA_TO_ORGAN.get(c.upper(), []):
            seen.add(kw)
    return list(seen)
