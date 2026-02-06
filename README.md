# Senary Bio

> **Computational Discovery Platform for Cas13d Therapeutics**

---

<p align="center">
  <marquee behavior="scroll" direction="left" scrollamount="4" style="background: linear-gradient(90deg, #0d1117 0%, #161b22 50%, #0d1117 100%); padding: 10px 0; border-radius: 4px; font-family: monospace; font-size: 14px;">
    <strong>Python 3.8+</strong> · <strong>BioPython</strong> · <strong>Pandas</strong> · <strong>NumPy</strong> · <strong>h5py</strong> · <strong>ARCHS4</strong> · <strong>NCBI Entrez</strong> · <strong>ESM-2</strong> · <strong>Transformers</strong> · <strong>PyTorch</strong> · <strong>Google Gemini</strong> · <strong>JupyterLab</strong> · <strong>SQLite</strong> · <strong>Ollama/Llama</strong> · <strong>ChimerDB</strong> · <strong>TCGA</strong>
  </marquee>
</p>

<div align="center">

**Tech Stack Bar:** Python · BioPython · Pandas · NumPy · h5py · ARCHS4 · NCBI Entrez · ESM-2 · Transformers · PyTorch · Google Gemini · JupyterLab · SQLite · Ollama/Llama · ChimerDB · TCGA

</div>

---

## 🎯 Mission

**Senary Bio** is pioneering the development of precision cancer therapeutics through novel Cas13d enzyme discovery. Our mission is to identify and validate high-collateral Cas13d variants capable of targeting cancer-specific fusion RNAs, creating a "suicide switch" mechanism that induces apoptosis specifically in tumor cells while preserving healthy tissue.

We combine computational biology, machine learning, and high-throughput screening to accelerate the discovery of next-generation RNA-guided therapeutics.

---

## 🗺️ Full Filtering Map

See **[FILTERING_MAP.md](FILTERING_MAP.md)** for a complete map from NCBI scraping → enzyme filters → target filters → matchmaker → expert agent → final Type VI Cas13d candidates.

---

## 🏗️ System Architecture & Logic Flow

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         COLLATERAL BIO PLATFORM                                   │
└─────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────┐   ┌─────────────────────────────────────────────┐
│         THE ARMORY (Enzyme Mining)           │   │        THE VAULT (Target Discovery)         │
├─────────────────────────────────────────────┤   ├─────────────────────────────────────────────┤
│  ncbi_miner.py      → NCBI Protein search    │   │  fusion_metadata.py  → Fusion → TCGA cancers│
│  sra_scout.py       → WGS + BioProject       │   │  specificity_filter  → Tissue specificity   │
│  autonomous_prospector → AI-driven deep mine │   │  mutation_loader     → VCF mutation mining  │
│  deep_miner_utils   → ESM-2 + CRISPR detect  │   │  archs4_loader       → Expression + safety  │
│  hepn_filter        → HEPN motif validation  │   │  run_targeting       → ARCHS4 query test    │
│  debug_sra          → NCBI connectivity test │   │                                             │
└──────────────────────┬──────────────────────┘   └──────────────────────┬──────────────────────┘
                       │                                                  │
                       ▼                                                  ▼
              ┌────────────────┐                                 ┌────────────────┐
              │ Cas13d Variants│                                 │ Fusion Targets │
              │ (FASTA / deep) │                                 │ (Specificity   │
              │                │                                 │  filtered)     │
              └────────┬───────┘                                 └────────┬───────┘
                       │                                                  │
                       └────────────────────┬─────────────────────────────┘
                                            │
                                            ▼
                               ┌────────────────────────┐
                               │   THE MATCHMAKER       │
                               │   (modules/matchmaker) │
                               │                        │
                               │  • Enzyme × Fusion     │
                               │  • PFS rule (cut sites)│
                               │  • Disease mapping     │
                               └────────────┬───────────┘
                                            │
                                            ▼
                               ┌────────────────────────┐
                               │  lead_candidates.csv   │
                               └────────────┬───────────┘
                                            │
                                            ▼
                               ┌────────────────────────┐
                               │   EXPERT AGENT         │
                               │   (analysis/expert_    │
                               │    agent.py)           │
                               │                        │
                               │  • ARCHS4 safety       │
                               │  • Organ-specific      │
                               │  • Gemini AI verdict   │
                               └────────────┬───────────┘
                                            │
                                            ▼
                               ┌────────────────────────┐
                               │ lead_candidates_       │
                               │ filtered.csv           │
                               │ dashboard.html         │
                               └────────────────────────┘
```

---

## 📋 Detailed Logic & Module Flow

### 1. Enzyme Mining (The Armory)

| Module | Purpose | Logic |
|--------|---------|-------|
| **ncbi_miner** | Annotated Cas13d from NCBI Protein | `Entrez.esearch(db="protein")` → fetch FASTA → save `search_YYYYMMDD.fasta` |
| **sra_scout** | Unannotated metagenomes (WGS) | Normalizes query, tries `wgs[Prop]` → fallback broader search → BioProject elink; 6-frame translate, HEPN `R.{4,6}H` + topology (100–600 aa spacing); saves `undiscovered_cas13d_*.fasta` |
| **autonomous_prospector** | AI-driven continuous mining | LLM formulates env query → SRAScout.search_wgs → semantic filter (LLM picks top datasets) → DeepEngine (ESM-2) + NeighborhoodWatch (CRISPR) → deep_mine ORFs 800–1100 aa; SQLite `visited_ids` to avoid re-processing; saves `deep_hits_*.fasta` |
| **deep_miner_utils** | Deep learning engine | **DeepEngine**: ESM-2 35M, cosine similarity vs Cas13d reference; **NeighborhoodWatch**: CRISPR array detection (24/28/32 bp chunks, 2–3 repeats) |
| **hepn_filter** | HEPN motif validation | Scans FASTA for ≥2 `R.{4}H` motifs → retains valid enzymes |
| **debug_sra** | Connectivity check | Tests NCBI fetch with known ID (E. coli) to verify network + translation |

### 2. Target Discovery (The Vault)

| Module | Purpose | Logic |
|--------|---------|-------|
| **fusion_metadata** | Fusion → cancers mapping | Loads `KB_and_Pub_Recur_per_cancer.csv` + novel matrix; builds `fusion → [TCGA]`; `TCGA_TO_ORGAN` maps cancer codes to ARCHS4 keywords |
| **specificity_filter** | High-specificity targets | Loads disease matrix (rows=cancer, cols=fusion); keeps fusions in ≤`max_tissue_types` (default 3); outputs `high_specificity_targets.csv` |
| **mutation_loader** | VCF mutation mining | Parses VCF for gene-specific mutations (e.g. KRAS G12C) for validation |
| **archs4_loader** | Expression & safety | HDF5 human_matrix; `get_gene_expression`, `get_gene_expression_normal_vs_cancer`, `fusion_absent_in_normal_present_in_cancer`; organ-specific mode uses enrichment factor |

### 3. Matchmaker

- Loads enzymes (FASTA or mock) and targets (`high_specificity_targets.csv` or `known_fusions.csv`)
- Disease map from `KB_and_Pub_Recur_per_cancer.csv` or `disease_matrix_*.csv` if no `Primary_Disease`
- Screens enzyme × target; PFS rule (no G at 3′); outputs `lead_candidates.csv`

### 4. Expert Agent

- Loads `lead_candidates.csv`, filters by `Associated_Disease`
- Groups by (Target_Fusion, Associated_Disease) to minimize API calls
- ARCHS4: organ-specific enrichment or global absent-in-normal
- Gemini AI verdict (GO / NO-GO / HOLD), screening strategy
- Outputs `lead_candidates_filtered.csv`

---

## 🛠️ Technology Stack

### Core

- **Python 3.8+** – Primary language
- **BioPython** – NCBI Entrez, SeqIO, FASTA
- **Pandas** – Fusion targets, expression, matrix ops
- **NumPy** – Numerical ops
- **h5py** – ARCHS4 HDF5
- **JupyterLab** – Exploration
- **openpyxl** – ChimerDB Excel
- **python-dotenv** – `.env` config

### Deep Learning (Autonomous Prospector)

- **PyTorch** – ESM-2
- **Transformers** – `facebook/esm2_t12_35M_UR50D`
- **requests** – Ollama/local LLM API

### AI & Data

- **Google Gemini** – Expert agent
- **Ollama / Llama** – Local LLM for prospector (optional)
- **SQLite** – Prospector history & visited IDs

### Data Sources

- **NCBI** – Protein, Nucleotide, BioProject
- **ARCHS4** – Human expression
- **ChimerDB** – Fusion RNAs
- **TCGA** – Cancer codes

---

## 🚀 Quick Start

### Prerequisites

```bash
python -m venv venv
# Windows: venv\Scripts\activate
# macOS/Linux: source venv/bin/activate
pip install -r requirements.txt
```

For **Autonomous Prospector** (optional):

```bash
pip install torch transformers requests
```

### Data Setup

| File | Description |
|------|-------------|
| `known_fusions.csv` | Validation targets |
| `novel_fusions.csv` | Discovery targets |
| `disease_matrix_known.csv` / `KB_and_Pub_Recur_per_cancer.csv` | Fusion × cancer matrix |
| `data/expression_data/human_matrix.h5` | ARCHS4 (download from [ARCHS4](https://maayanlab.cloud/archs4/)) |

Regenerate CSVs from Excel:

```bash
python utils/split_excel.py
```

---

## 📂 Workflow Commands

### Step 1: Mine Enzymes

**NCBI Protein (annotated):**
```bash
python -c "from modules.mining.ncbi_miner import EnzymeMiner; EnzymeMiner().search_and_fetch('Cas13d')"
```

**SRA Scout (WGS metagenomes):**
```bash
python -c "
from modules.mining.sra_scout import SRAScout
scout = SRAScout()
ids = scout.search_wgs('hydrothermal vent metagenome', max_records=50)
candidates = scout.fetch_and_mine(ids)
scout.save_discoveries(candidates)
"
```

**Autonomous Prospector (AI + ESM-2):**
```bash
# Requires: torch, transformers, LLM_LOCAL_URL (Ollama) or LLM_PROVIDER
python modules/mining/autonomous_prospector.py
```

**Debug NCBI:**
```bash
python modules/mining/debug_sra.py
```

### Step 2: Family Grouping (mined sequences)

Group mined sequences from `data/mined_sequences/` into families by homology (ESM-2) and HEPN count. Naming: `SN01_001`, `SN01_002`, etc. Output: `data/fam_fasta.fasta` (used by matchmaker).

```bash
# Requires: torch, transformers, scipy
python modules/mining/family_grouper.py
```

Options: `--input-dir`, `--output-dir`, `--threshold 0.7`, `--prefix SN`, `--glob "*.fasta"`

Env: `FAMILY_DEVICE=cuda` (or `cpu`) to force device; `EMBED_BATCH_SIZE=50` for memory tuning; AMD GPU requires PyTorch ROCm build (Linux).

### Step 3: Specificity Filter (optional)

```bash
python modules/targeting/specificity_filter.py
# Uses disease_matrix_novel.csv → data/high_specificity_targets.csv
```

### Step 4: Matchmaker

```bash
python modules/matchmaker.py
# Uses high_specificity_targets.csv if present, else known_fusions.csv
```

### Step 5: Expert Agent

```bash
# .env: GEMINI_API_KEY
python modules/analysis/expert_agent.py
# → lead_candidates_filtered.csv
```

### Step 6: ARCHS4 Query Test

```bash
python run_targeting.py
```

### Structure Pipeline (OmegaFold)

For 2-3 HEPN filtered enzymes: structure prediction, TM-score homology vs Cas13a/b, and 3D dashboard.

```bash
# 1. Filter to 2-3 HEPN sequences
python visualization/filter_23_hepn.py

# 2. Run OmegaFold (PyTorch; works on RunPod, local GPU)
python visualization/run_omegafold.py

# 3. Compute TM-score vs Cas13 references
python visualization/run_tmscore.py

# 4. Generate dashboard
python visualization/structure_dashboard.py

# 5. Serve and view
python -m http.server 8000
# Open: http://localhost:8000/visualization/structure_dashboard.html
```

Requires: `tmtools` (or US-align), `biopython`. OmegaFold is not on PyPI and only supports Python 3.8–3.10. On Python 3.11/3.12 (e.g. RunPod), clone and run from source:

```bash
git clone https://github.com/HeliXonProtein/OmegaFold.git
cd OmegaFold && pip install torch biopython
```

Then run the structure pipeline with:

```bash
python visualization/run_omegafold.py --omegafold-repo /path/to/OmegaFold
```

(Or set `OMEGAFOLD_REPO=/path/to/OmegaFold`.) OmegaFold outputs to `data/structure_pipeline/structures/omegafold/`.

---

## 📂 Project Structure

```
collateral_bio_core/
├── README.md
├── PIPELINE.md
├── requirements.txt
├── run_targeting.py           # ARCHS4 loader test
├── main.py
├── data/
│   ├── raw_sequences/         # Mined FASTA
│   ├── mined_sequences/       # Novel hits (deep_hits_*.fasta) → family_grouped_*.fasta
│   ├── expression_data/       # human_matrix.h5
│   ├── high_specificity_targets.csv
│   ├── known_fusions.csv, novel_fusions.csv
│   ├── disease_matrix_known.csv, disease_matrix_novel.csv
│   └── KB_and_Pub_Recur_per_cancer.csv
├── modules/
│   ├── mining/
│   │   ├── ncbi_miner.py
│   │   ├── sra_scout.py
│   │   ├── autonomous_prospector.py
│   │   ├── family_grouper.py   # Groups mined seqs by homology (SN01_001, etc.)
│   │   ├── deep_miner_utils.py
│   │   ├── hepn_filter.py
│   │   └── debug_sra.py
│   ├── targeting/
│   │   ├── archs4_loader.py
│   │   ├── fusion_metadata.py
│   │   ├── specificity_filter.py
│   │   └── mutation_loader.py
│   ├── discovery/
│   │   └── fusion_caller.py
│   ├── analysis/
│   │   └── expert_agent.py
│   └── matchmaker.py
├── utils/
│   ├── split_excel.py
│   ├── logger.py
│   └── inspect_archs4_metadata.py
├── prompts/
│   └── expert_persona.txt
├── lead_candidates.csv
└── lead_candidates_filtered.csv
```

---

## ⚙️ Environment Variables

| Variable | Purpose |
|----------|---------|
| `GEMINI_API_KEY` | Expert agent AI |
| `TARGET_FUSIONS_CSV` | `novel_fusions.csv` for novel run |
| `NORMAL_MAX_TPM`, `CANCER_MIN_TPM` | ARCHS4 filter thresholds |
| `ENRICHMENT_FACTOR` | Organ-specific enrichment (default 2.0) |
| `USE_ORGAN_SPECIFIC` | 1 = organ-specific, 0 = global |
| `LLM_PROVIDER`, `LLM_LOCAL_URL`, `LLM_MODEL` | Prospector LLM (e.g. Ollama) |
| `DEEP_MINE_MAX`, `ESM_THRESHOLD`, `REQUIRE_CRISPR` | Prospector tuning |

---

## ⚠️ Troubleshooting

- **"Enzyme file not found"** – Matchmaker falls back to mock enzymes.
- **"ARCHS4 file not found"** – Download `human_matrix.h5` into `data/expression_data/`.
- **"Column not found"** – Matchmaker handles `fusionsss` typo; verify CSV headers.
- **Prospector import error** – Install `torch`, `transformers`, `requests`.

---

## 📄 License

Proprietary – Collateral Bio © 2026

---

*Built with precision for precision medicine.*
