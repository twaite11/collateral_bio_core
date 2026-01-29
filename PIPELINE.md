# Collateral Bio – Full Pipeline Walkthrough

This guide walks you through the **exact commands** to run the full pipeline using your **matrix file** (`KB_and_Pub_Recur_per_cancer.csv`) and **fusion CSVs** (`known_fusions.csv`, `novel_fusions.csv`).

**Quick start:** If those files are in `data/`, ARCHS4 is in `data/expression_data/human_matrix.h5`, and `.env` has `GEMINI_API_KEY`, run:

```bash
python run_pipeline.py
```

Otherwise, follow the steps below.

---

## 1. Prerequisites

- **Python 3.8+**
- **Data files** in `data/`:
  - **Matrix:** `data/KB_and_Pub_Recur_per_cancer.csv` (or `disease_matrix_known.csv`) – fusions × cancer types
  - **Fusion CSVs:** `data/known_fusions.csv` and/or `data/novel_fusions.csv` – fusion lists for targeting
  - **ARCHS4 (Expert Agent):** `data/expression_data/human_matrix.h5` – normal vs cancer expression
- **`.env`** in project root with `GEMINI_API_KEY=...` (for Expert Agent AI)

---

## 2. One-time setup

From the project root (`collateral_bio_core/`):

```bash
# Create and activate venv (Windows)
python -m venv venv
venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

---

## 3. (Optional) Refresh fusion CSVs from Excel

If you use `Recurrent_table.xlsx` and want to regenerate `known_fusions.csv` / `novel_fusions.csv`:

```bash
python utils/split_excel.py
```

This reads `data/Recurrent_table.xlsx` and writes the fusion CSVs and matrix CSVs (e.g. `KB_and_Pub_Recur_per_cancer.csv`) into `data/`.

---

## 4. Place ARCHS4 matrix (Expert Agent)

Download **human_matrix.h5** from [ARCHS4](https://maayanlab.cloud/archs4/) and put it here:

```
data/expression_data/human_matrix.h5
```

Create the folder if needed:

```bash
mkdir data\expression_data
# then copy human_matrix.h5 into it
```

---

## 5. Run the full pipeline

**Option A – Single command (all steps):**

```bash
python run_pipeline.py
```

This runs, in order:

1. **Mine** – NCBI Cas13d search → `data/raw_sequences/search_YYYYMMDD_HHMMSS.fasta`
2. **Matchmaker** – enzymes × fusions (using your matrix + fusion CSVs) → `lead_candidates.csv`
3. **Expert Agent** – cancer-only filter, safety, AI → `dashboard.html`, `lead_candidates_filtered.csv`

---

**Option B – Step-by-step commands:**

**Step 1 – Mine Cas13d enzymes (NCBI):**

```bash
python -c "from modules.mining.ncbi_miner import EnzymeMiner; m = EnzymeMiner('data/raw_sequences'); m.search_and_fetch(query='Cas13d NOT synthetic construct', max_results=50)"
```

**Step 2 – Matchmaker (enzymes vs fusions):**

Uses latest `data/raw_sequences/search_*.fasta`, `data/known_fusions.csv`, and `data/KB_and_Pub_Recur_per_cancer.csv`.

```bash
python modules/matchmaker.py
```

**Step 3 – Expert Agent (filter + safety + AI):**

Reads `lead_candidates.csv`, uses ARCHS4 for normal/cancer filtering, then AI commentary.

```bash
python modules/analysis/expert_agent.py
```

---

## 6. Choosing fusion CSV (known vs novel)

- **Known fusions (default):** `data/known_fusions.csv`
- **Novel fusions:** set env before running matchmaker / pipeline:

```bash
set TARGET_FUSIONS_CSV=novel_fusions.csv
python modules/matchmaker.py
```

Or for the full pipeline:

```bash
set TARGET_FUSIONS_CSV=novel_fusions.csv
python run_pipeline.py
```

On macOS/Linux use `export TARGET_FUSIONS_CSV=novel_fusions.csv` instead of `set`.

---

## 7. Outputs

| Step         | Output                        | Description                                        |
|-------------|-------------------------------|----------------------------------------------------|
| Mine        | `data/raw_sequences/search_*.fasta` | Cas13d sequences from NCBI                         |
| Matchmaker  | `lead_candidates.csv`         | Enzyme–fusion pairs (ranked)                       |
| Expert Agent| `lead_candidates_filtered.csv`| Cancer-only, absent-in-normal filtered leads       |
| Expert Agent| `dashboard.html`              | Dashboard of filtered candidates (open in browser) |

---

## 8. Inspect ARCHS4 normal/cancer classification

To verify that ARCHS4 samples are correctly classified as normal vs cancer (and tune keywords if needed):

```bash
python utils/inspect_archs4_metadata.py
```

This prints unique `source_name_ch1` labels and their classification (normal / cancer / unknown).

## 9. Env options for Expert Agent filter

- **`ENRICHMENT_FACTOR`** (default `2.0`): organ-specific mode requires `cancer_mean ≥ enrichment_factor × normal_mean`.
- **`USE_ORGAN_SPECIFIC`** (default `1`): use organ-specific normal/cancer and enrichment when fusions have associated cancers (from KB_and_Pub / novel_Recur). Set to `0` to use global normal/cancer only.

## 10. Troubleshooting

- **"No module named 'Bio'"** → `pip install biopython` (or `pip install -r requirements.txt`).
- **"ARCHS4 required"** → Add `data/expression_data/human_matrix.h5` (see step 4).
- **"Disease matrix file not found"** → Ensure `data/KB_and_Pub_Recur_per_cancer.csv` exists (or `disease_matrix_known.csv`).
- **"Target file not found"** → Ensure `data/known_fusions.csv` (or `novel_fusions.csv`) exists; run `utils/split_excel.py` if needed.
- **No fusions pass filter** → Run `utils/inspect_archs4_metadata.py` to check normal/cancer labels. Ensure fusions have associated cancers (matchmaker uses `KB_and_Pub_Recur_per_cancer`). Try lowering `ENRICHMENT_FACTOR` or raising `CANCER_MIN_TPM`.
