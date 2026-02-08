# Running the Full Pipeline on a VPS

This guide describes how to run the full 4-step discovery pipeline (mine → embed & mutate → structure filter → identity filter → matchmaker) on a VPS (e.g. RunPod, Lambda, or a bare Ubuntu server).

---

## 1. VPS Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| **CPU** | 4 cores | 8+ cores |
| **RAM** | 16 GB | 32 GB |
| **GPU** | Optional for mining/design; **required for OmegaFold** | 1× NVIDIA (e.g. A10, RTX 3090) |
| **Disk** | 50 GB | 100+ GB (for structures + embeddings) |
| **OS** | Ubuntu 20.04 / 22.04 | Ubuntu 22.04 |

- **Mining + design (ESM-2):** Runs on CPU or GPU; GPU speeds up embedding and mutation scoring.
- **Structure filter (OmegaFold):** Needs a GPU; otherwise use `--skip-structure` and run OmegaFold elsewhere.

---

## 2. One-Time Setup on the VPS

### 2.1 System dependencies

```bash
sudo apt-get update
sudo apt-get install -y python3.10 python3.10-venv git
# For US-align (TM-score) if you run structure filter locally
# sudo apt-get install -y unzip wget
# Then install US-align from https://zhanggroup.org/US-align/
```

### 2.2 Clone repo and venv

```bash
cd /opt  # or your preferred path
git clone https://github.com/your-org/collateral_bio_core.git
cd collateral_bio_core
python3.10 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 2.3 OmegaFold (only if running structure step on this VPS)

OmegaFold supports Python 3.8–3.10. If the project uses Python 3.11+, use a separate venv for OmegaFold:

```bash
# Option A: Separate OmegaFold venv (Python 3.10)
python3.10 -m venv venv_omegafold
source venv_omegafold/bin/activate
pip install torch biopython
git clone https://github.com/HeliXonProtein/OmegaFold.git /opt/OmegaFold
# Download weights per OmegaFold README
export OMEGAFOLD_REPO=/opt/OmegaFold
deactivate
```

Then in the main pipeline, set `OMEGAFOLD_REPO=/opt/OmegaFold` and call OmegaFold via subprocess (the pipeline script uses `python main.py` from that repo).

### 2.4 Reference data

- **Known Cas13 (for identity/drift):** Create `data/references/known_cas13.fasta` with sequences of LwaCas13a, RfxCas13d, etc. (one FASTA per reference). If missing, identity filter treats all as 0% identity.
- **Reference PDBs (for structure filter):** Download once:
  ```bash
  python visualization/run_tmscore.py --skip-download  # still creates dir
  # Then download 5W1H.pdb and 6DTD.pdb into data/structure_pipeline/references/
  # Or run run_tmscore without --skip-download to fetch from RCSB.
  ```
- **Targets for matchmaker:** Ensure `data/targets/high_specificity_targets.csv` or `data/targets/known_fusions.csv` and `data/matrices/` disease matrix CSVs exist (see README).

### 2.5 Environment file

```bash
cp config/pipeline.env.example .env
# Edit .env: set ESM_SIMILARITY_CEILING=0.82, OMEGAFOLD_REPO, FAMILY_DEVICE=cuda, etc.
# Full enzyme + CRISPR: REQUIRE_START_M=1, MIN_CTERM_TAIL=15, REQUIRE_FULL_STRUCTURE=1, MIN_REPEAT_COUNT=1
```

### 2.6 Full enzyme and CRISPR repeats (optional)

To keep only **full-enzyme** hits (N-term M, C-term tail after HEPN, not truncated at contig boundary) and **full locus** (CRISPR array + at least one repeat domain saved for synthesis):

- `REQUIRE_START_M=1` – ORF must start with Met (default on).
- `MIN_CTERM_TAIL=15` – Minimum residues after last HEPN motif (default 15).
- `CONTIG_BOUNDARY_MARGIN=30` – Reject ORFs that end within 30 nt of contig end (default 30).
- `REQUIRE_FULL_STRUCTURE=1` – When set with `REQUIRE_CRISPR=1`, also require `len(repeat_domains) >= MIN_REPEAT_COUNT` so every hit has CRISPR repeats saved.
- `MIN_REPEAT_COUNT=1` – Minimum repeat sequences to keep a hit (default 1).

### 2.7 Mining for closest similarity to RfxCas13d / PspCas13a

Set **`ESM_REFERENCE_FASTA=data/references/mining_refs.fasta`** (or path to a FASTA with **RfxCas13d** and **PspCas13a**). The miner scores each ORF against both refs and keeps the **closest** match (max ESM similarity). Leave **`ESM_SIMILARITY_CEILING`** unset so high-similarity hits are not filtered out. Replace placeholder sequences in `mining_refs.fasta` with full sequences from UniProt (e.g. RfxCas13d from PDB 6IV9) for best results.

### 2.8 Design: trans-cleavage mutation suggestions

In **mutate_for_drift** (or **run_full_pipeline**), use **`--use-trans-cleavage-prompt`** (or **`USE_TRANS_CLEAVAGE_PROMPT=1`**) and set **`GEMINI_API_KEY`**. The pipeline will ask Gemini for mutations that maintain structural stability to RfxCas13d or PspCas13a but might increase trans-cleavage activity; suggested mutants are validated with ESM-2 vs the same refs.

---

## 3. Pipeline Run Modes

### Mode A: Mining only (continuous)

Run the autonomous prospector to fill `data/raw_sequences/deep_hits_*.fasta` and metadata. No GPU required for mining; GPU speeds up ESM-2 scoring.

```bash
source venv/bin/activate
export ESM_SIMILARITY_CEILING=0.82   # diversity mode: favor distant homologs
python modules/mining/autonomous_prospector.py
# Runs until stopped. Check data/raw_sequences/ and data/prospector.db.
```

### Mode B: Full pipeline (post-mine) on latest FASTA

After you have at least one `deep_hits_*.fasta` (or any pool FASTA), run embed → mutate → structure → identity, then optionally matchmaker.

**With GPU (OmegaFold on this machine):**

```bash
source venv/bin/activate
python run_full_pipeline.py \
  --input data/raw_sequences/deep_hits_20260204_010346.fasta \
  --max-identity 0.85 \
  --tm-threshold 0.4 \
  --run-matchmaker
```

**Without GPU (skip structure filter):**

```bash
python run_full_pipeline.py \
  --input data/raw_sequences/deep_hits_20260204_010346.fasta \
  --skip-structure \
  --max-identity 0.85 \
  --run-matchmaker
```

**Using latest deep_hits automatically:**

```bash
python run_full_pipeline.py --run-matchmaker
# Picks latest data/raw_sequences/deep_hits_*.fasta or data/mined_sequences/deep_hits_*.fasta
```

### Mode C: Design + identity only (no structure)

If you already have a FASTA (e.g. merged mined sequences) and want to run design and identity filter only:

```bash
python run_full_pipeline.py \
  --input data/mined_sequences/family_grouped_20260205.fasta \
  --skip-structure \
  --run-matchmaker
```

### Mode D: Identity filter only (pre-filtered FASTA)

If you already ran structure filter elsewhere and have `passed_structures.fasta`:

```bash
python -c "
from modules.identity_filter import run_identity_filter
run_identity_filter(
  'data/structure_pipeline/passed_structures.fasta',
  'data/references/known_cas13.fasta',
  'data/identity_filtered/passed.fasta',
  'data/identity_filtered/identity_metadata.csv',
  max_identity=0.85
)
"
```

---

## 4. Recommended VPS Workflow

1. **Start mining in background (tmux/screen):**
   ```bash
   tmux new -s mining
   source venv/bin/activate
   export ESM_SIMILARITY_CEILING=0.82
   python modules/mining/autonomous_prospector.py
   # Detach: Ctrl+B, D
   ```

2. **Periodically run the full pipeline on latest hits:**
   ```bash
   python run_full_pipeline.py --run-matchmaker
   ```
   Or run nightly via cron:
   ```cron
   0 2 * * * cd /opt/collateral_bio_core && source venv/bin/activate && python run_full_pipeline.py --run-matchmaker >> logs/pipeline.log 2>&1
   ```

3. **If OmegaFold is on another machine:** Run mining + design + identity on the VPS with `--skip-structure`. Copy `data/design/drift_variants.fasta` to the GPU box, run OmegaFold and TM-score there, copy back `passed_structures.fasta`, then run identity filter locally.

---

## 5. Key Paths and Outputs

| Step | Input | Output |
|------|--------|--------|
| Mining | NCBI WGS (via autonomous_prospector) | `data/raw_sequences/deep_hits_*.fasta`, `*_metadata.csv` |
| Embed | Pool FASTA | `data/design/embeddings/embeddings.npy`, `sequence_ids.txt` |
| Mutate | Pool FASTA | `data/design/drift_variants.fasta` |
| Structure | `drift_variants.fasta` | `data/structure_pipeline/structures/omegafold/*.pdb`, `passed_structures.fasta` |
| Identity | `passed_structures.fasta` | `data/identity_filtered/passed.fasta`, `identity_metadata.csv` |
| Matchmaker | `passed.fasta` + targets | `lead_candidates.csv` |

---

## 6. Troubleshooting

- **Out of memory (ESM-2):** Set `EMBED_BATCH_SIZE=25` or `FAMILY_DEVICE=cpu`.
- **OmegaFold not found:** Set `OMEGAFOLD_REPO` to the clone path or use `--skip-structure`.
- **No reference PDBs:** Run `visualization/run_tmscore.py` once (without `--skip-download`) to fetch 5W1H and 6DTD.
- **Identity filter keeps all:** Add sequences to `data/references/known_cas13.fasta`; otherwise max identity is 0 and all pass.
- **Mining returns no hits:** Lower `ESM_THRESHOLD` (e.g. 0.65) or set `REQUIRE_CRISPR=0`; check NCBI/network access.
