# Collateral Bio

> **Computational Discovery Platform for Cas13d Therapeutics**

---

## 🎯 Mission

**Collateral Bio** is pioneering the development of precision cancer therapeutics through novel Cas13d enzyme discovery. Our mission is to identify and validate high-collateral Cas13d variants capable of targeting cancer-specific fusion RNAs, creating a "suicide switch" mechanism that induces apoptosis specifically in tumor cells while preserving healthy tissue.

We combine computational biology, machine learning, and high-throughput screening to accelerate the discovery of next-generation RNA-guided therapeutics.

---

## 🏗️ System Architecture

The Collateral Bio platform operates as an integrated discovery pipeline:

```
┌─────────────────────────────────────────────────────────────┐
│                    COLLATERAL BIO PLATFORM                  │
└─────────────────────────────────────────────────────────────┘

┌──────────────────────┐         ┌──────────────────────┐
│   THE ARMORY         │         │   THE VAULT          │
│   (Enzyme Mining)    │         │   (Target Discovery) │
│                      │         │                      │
│  • NCBI Database     │         │  • ChimerDB          │
│  • Metagenomes       │         │  • TCGA Data         │
│  • Public Sequences  │         │  • Fusion RNAs       │
└──────────┬───────────┘         └──────────┬───────────┘
           │                                 │
           v                                 v
    ┌──────────────┐                 ┌──────────────┐
    │ Cas13d       │                 │ Fusion RNA   │
    │ Variants     │                 │ Targets      │
    └──────┬───────┘                 └──────┬───────┘
           │                                 │
           └────────────┬────────────────────┘
                        │
                        v
            ┌───────────────────────┐
            │   THE MATCHMAKER      │
            │   (Virtual Wet Lab)   │
            │                       │
            │  • Druggability       │
            │  • Safety Profiling   │
            │  • Market Analysis    │
            └───────────┬───────────┘
                        │
                        v
            ┌───────────────────────┐
            │  lead_candidates.csv  │
            │  (Series Seed Assets) │
            └───────────────────────┘
```

---

## 🛠️ Technology Stack

### Core Technologies
- **Python 3.8+** - Primary development language
- **BioPython** - NCBI Entrez API integration, sequence parsing, FASTA handling
- **Pandas** - Data manipulation and analysis for fusion targets and expression data
- **NumPy** - Numerical computations and matrix operations
- **h5py** - ARCHS4 expression database access (HDF5 format)
- **JupyterLab** - Interactive data exploration and visualization
- **openpyxl** - Excel file processing for ChimerDB data extraction

### Data Sources
- **NCBI Protein Database** - Enzyme sequence mining
- **ChimerDB** - Fusion RNA target database
- **ARCHS4** - Human tissue expression atlas (safety profiling)
- **TCGA** - Cancer genomics data integration

### Infrastructure
- **Git** - Version control and IP tracking
- **Virtual Environments** - Dependency isolation
- **HDF5** - Efficient storage of large expression matrices

---

## 🚀 Quick Start

### Prerequisites

Ensure you have **Python 3.8+** installed on your system.

```bash
# Create virtual environment
python -m venv venv

# Activate virtual environment
# On macOS/Linux:
source venv/bin/activate

# On Windows:
venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Data Setup

The platform requires specific data files in the `data/` directory:

| File Name | Description | Source |
|-----------|-------------|--------|
| `known_fusions.csv` | Validation set of proven targets (e.g., TMPRSS2-ERG) | Extracted from ChimerDB (HRS_Recurrent_known) |
| `novel_fusions.csv` | Discovery set of novel targets (e.g., SPSB4-ACPL2) | Extracted from ChimerDB (HRS_Recurrent_novel) |
| `disease_matrix_known.csv` | Disease mapping: fusions to cancer types | Extracted from ChimerDB (KB_and_Pub_Recur) |
| `human_matrix.h5` | Safety map: expression across 50k+ tissues | Download from [ARCHS4](https://maayanlab.cloud/archs4/) |

> **Note:** If starting fresh with the raw Excel file (`Recurrent_table.xlsx`), run:
> ```bash
> python utils/split_excel.py
> ```
> This will auto-generate the required CSV files.

---

## 📋 Workflow

**Full command walkthrough** (matrix + fusion CSVs): see **[PIPELINE.md](PIPELINE.md)** for exact commands to run the full pipeline step-by-step or via `python run_pipeline.py`.

---

### Step 1: Mine the Warheads (Enzymes)

Discover novel Cas13d variants from public sequence databases:

```bash
python -c "from modules.mining.ncbi_miner import EnzymeMiner; EnzymeMiner().search_and_fetch('Cas13d')"
```

**Output:** Timestamped FASTA files saved to `data/raw_sequences/`

### Step 2: The Matchmaker (Simulation)

The core discovery engine screens enzymes against fusion targets, evaluating:
- **Patient Count** - Market size validation
- **Disease Link** - Cancer type association
- **Druggability** - PFS rule validation (cut site analysis)

**Known fusions (default):** uses `data/known_fusions.csv` and `data/KB_and_Pub_Recur_per_cancer.csv`.
```bash
python modules/matchmaker.py
```

**Novel fusions:** set `TARGET_FUSIONS_CSV=novel_fusions.csv`, then run as above.

**Output:** `lead_candidates.csv` containing ranked therapeutic enzyme-target pairs

### Step 3: AI CSO Agent & Cancer-Only Filter

The expert agent **filters** lead candidates so the **final output includes only fusions that are absent in normal healthy tissues and present in cancer samples**. It then runs safety checks and AI commentary.

**Requirements:** `data/expression_data/human_matrix.h5` (ARCHS4) and `GEMINI_API_KEY` in `.env`.

```bash
python modules/analysis/expert_agent.py
```

**Output:**
- `dashboard.html` – Ranked candidates (cancer-only, absent in normal)
- `lead_candidates_filtered.csv` – Filtered leads in CSV form

**Env options:**  
- `NORMAL_MAX_TPM` (default `1.0`), `CANCER_MIN_TPM` (default `1.0`): exclusion thresholds.  
- `MAX_CANDIDATES` (default: unset = **all**): max transcript rows to search; increase or leave unset to avoid missing hits.  
- `BATCH_SIZE` (default `200`): rows per batch when scanning.  

Fusions are **excluded** if either parent gene is expressed above `NORMAL_MAX_TPM` in normal tissues or below `CANCER_MIN_TPM` in cancer (when cancer data exists). The agent processes candidates in batches; if none in a batch pass the no-normal filter, it **keeps searching** the next batch until it finds hits or exhausts the file.

### Step 4: Safety Audit (optional)

Validate a specific gene with the safety profiler:

```bash
# Edit the script to audit your specific gene
python modules/targeting/safety_profiler.py
```

---

## 📂 Project Structure

```
collateral-bio/
├── .gitignore                   # Ignore large data files (crucial for git)
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── main.py                      # Central command CLI
├── data/
│   ├── raw_sequences/           # Mined FASTA files (IP Origin)
│   ├── processed_enzymes/       # Filtered Cas13d candidates
│   └── expression_data/         # ARCHS4 human_matrix.h5 (Excluded from Git)
├── modules/
│   ├── __init__.py
│   ├── mining/                  # ENZYME DISCOVERY
│   │   ├── __init__.py
│   │   ├── ncbi_miner.py        # The Clean Room scraper
│   │   └── hepn_filter.py       # The "Scissors" detector
│   ├── targeting/               # TARGET DISCOVERY
│   │   ├── __init__.py
│   │   ├── archs4_loader.py     # The Binary Expression Finder
│   │   ├── chimerdb_loader.py   # Fusion RNA loader
│   │   └── ccle_filter.py       # Cancer cell line filter
│   ├── discovery/               # DISCOVERY MODULES
│   │   └── fusion_caller.py     # Fusion detection
│   └── matchmaker.py            # CORE: Enzyme-Target Matching Engine
├── utils/
│   ├── __init__.py
│   ├── split_excel.py           # ChimerDB data processor
│   └── logger.py                # Logging utilities
├── notebooks/                   # For visual exploration
│   └── 01_candidate_viz.ipynb
├── lead_candidates.csv          # Matchmaker output (pre-filter)
└── lead_candidates_filtered.csv # Final output: Cancer-only, absent-in-normal leads
```

---

## ⚠️ Troubleshooting

### Common Issues

**"Enzyme file not found"**
- The system defaults to "Mock Enzyme" simulation if the miner hasn't been run yet
- This is acceptable for testing pipeline logic
- Run Step 1 to generate real enzyme data

**"Column not found"**
- ChimerDB CSVs may have typos (e.g., `fusionsss`)
- The matchmaker includes built-in typo handling
- Verify CSV headers if data loading fails

**"ARCHS4 file not found"**
- Download `human_matrix.h5` from [ARCHS4](https://maayanlab.cloud/archs4/)
- Place in `data/expression_data/` directory
- File is large (~10GB) and excluded from git

---

## 📄 License

Proprietary - Collateral Bio © 2026

---

## 🤝 Contact

For inquiries about Collateral Bio's Cas13d discovery platform, please contact the development team.

---

*Built with precision for precision medicine.*
