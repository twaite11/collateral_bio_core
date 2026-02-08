# Collateral Bio: Full Filtering Map

**From NCBI Scraping → Final Type VI Cas13d Enzyme Candidates**

Wireframe map with filtering actions and criteria at each step.

---

## End-to-End Flow Overview

```
┌─────────────────────────────────────────────────────────────────────────────────────────────────┐
│  NCBI (Raw Data)  →  Enzyme Filters  →  Target Filters  →  Matchmaker  →  Expert Agent  →  FINAL │
└─────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Phase 1: NCBI Scraping & Raw Collection

```
                                    ┌─────────────────────────────────────────┐
                                    │           NCBI DATABASES                 │
                                    └───────────────────┬─────────────────────┘
                                                        │
            ┌───────────────────────────────────────────┼───────────────────────────────────────────┐
            │                                           │                                           │
            ▼                                           ▼                                           ▼
┌───────────────────────────────┐     ┌───────────────────────────────┐     ┌───────────────────────────────┐
│       ncbi_miner              │     │       sra_scout               │     │   autonomous_prospector       │
├───────────────────────────────┤     ├───────────────────────────────┤     ├───────────────────────────────┤
│ FILTER: Search term match     │     │ FILTER: Environment keywords  │     │ FILTER: LLM picks strategy;   │
│ "Cas13d" in NCBI Protein DB   │     │ e.g. "hydrothermal vent       │     │ semantic filter picks top N   │
│ NCBI returns only annotated   │     │ metagenome". Try wgs[Prop]    │     │ datasets most likely to have   │
│ Cas13d proteins.              │     │ first; fallback to broader    │     │ uncultured CRISPR / viral      │
│                               │     │ search or BioProject elink.   │     │ defense / diverse microbes.    │
└───────────────┬───────────────┘     └───────────────┬───────────────┘     └───────────────┬───────────────┘
                │                                     │                                     │
                ▼                                     ▼                                     ▼
┌───────────────────────────────┐     ┌───────────────────────────────┐     ┌───────────────────────────────┐
│ search_YYYYMMDD.fasta         │     │ undiscovered_cas13d_*.fasta   │     │ deep_hits_*.fasta             │
│ (annotated proteins)          │     │ (unannotated WGS hits)        │     │ (AI-validated deep hits)      │
└───────────────┬───────────────┘     └───────────────┬───────────────┘     └───────────────┬───────────────┘
                │                                     │                                     │
                └─────────────────────────────────────┼─────────────────────────────────────┘
                                                      │
                                                      ▼
                                    ┌─────────────────────────────────────────┐
                                    │   data/raw_sequences/                    │
                                    │   Raw FASTA files (input to Phase 2)    │
                                    └───────────────────┬─────────────────────┘
                                                        │
                                                        ▼
```

---

## Phase 2: Enzyme-Side Filtering (Type VI Cas13d Criteria)

```
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                           ENZYME FILTER PIPELINE                                                          │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘

   RAW SEQUENCES (FASTA)
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 1: SIZE (600–1400 aa)                                                                             │
│ Keep: Proteins between 600–1400 amino acids.                                                             │
│ Throw out: Too short (fragments) or too long (wrong enzyme family).                                      │
│ Why: Cas13d family diversity; broader range captures novel variants.                                     │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 2: HEPN MOTIF (R...H pattern)                                                                     │
│ Keep: 2–3 HEPN motifs (Arginine–4to6 aa–Histidine). Motifs must be 100–1200 aa apart.                    │
│ Throw out: Proteins missing the double-domain structure or with 4+ HEPN domains.                         │
│ Why: HEPN domains are the RNA-cutting "scissors"; Cas13d needs two or three, properly spaced.            │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ├───────────────────────────────────────────────────────────────────────┐
            │  [Prospector path only – sra_scout/ncbi_miner skip to Filter 5]       │
            ▼                                                                       │
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 3: CRISPR ARRAY (Prospector only)                                                                 │
│ Keep: Contigs with 2–3+ repeats (24–32 bp chunks).                                                        │
│ Throw out: Proteins without CRISPR context.                                                               │
│ Why: Cas13d usually sits next to CRISPR arrays in nature.                                                 │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 4: ESM-2 SIMILARITY (Prospector only)                                                             │
│ Keep: Cosine similarity to known Cas13d ref > 0.75.                                                      │
│ Throw out: Sequences that don't "look like" Cas13d.                                                       │
│ Why: Deep learning detects subtle sequence patterns.                                                      │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 5: Structure filter (run_full_pipeline)                                                          │
│ Keep: 2–3 HEPN motifs (R.{4,6}H) in sequence; TM-score vs Cas13 refs above threshold.                  │
│ Throw out: Wrong HEPN count or low structural similarity.                                                 │
│ Why: Bi-lobed/HEPN check in structure_filter (bi_lobed_hepn_check.py).                                    │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTERED CAS13D ENZYMES → data/raw_sequences/*.fasta                                                    │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            │  [Optional: Family Grouping]
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FAMILY GROUPING (mined_sequences)                                                                       │
│ Partition by HEPN count → cluster by ESM-2 homology → name SN01_001, SN01_002, ...                      │
│ Output: family_grouped_*.fasta, family_manifest_*.csv                                                   │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Phase 3: Target-Side Filtering

```
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                           TARGET FILTER PIPELINE                                                          │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘

   ChimerDB / Excel:
   data/targets/known_fusions.csv, novel_fusions.csv
   data/matrices/disease_matrix_known.csv, disease_matrix_novel.csv
   KB_and_Pub_Recur_per_cancer.csv
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 1: SPECIFICITY                                                                                    │
│ Keep: Fusions present in ≤3 cancer types (configurable).                                                │
│ Throw out: Fusions in many cancer types (promiscuous / artifact).                                       │
│ Why: Tissue-specific fusions are real drivers; pan-cancer = noise.                                      │
│ Output: data/targets/high_specificity_targets.csv                                                                     │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTER 2: FUSION METADATA (Mapping)                                                                      │
│ Action: Build fusion → [TCGA codes] and TCGA → organ keywords (e.g. PRAD→prostate, BRCA→breast).        │
│ Why: Needed for organ-specific ARCHS4 safety checks later.                                              │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FILTERED FUSION TARGETS (cancer-specific, tissue-limited)                                               │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Phase 4: Matchmaker (Enzyme × Target Screening)

```
   FILTERED ENZYMES                         FILTERED TARGETS
   (data/raw_sequences/*.fasta)             (high_specificity_targets.csv or known_fusions.csv)
            │                                               │
            └───────────────────┬───────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ MATCHMAKER FILTER 1: VALID CUT SITES (PFS RULE)                                                         │
│ Keep: Pairs where fusion junction has valid cut sites (no G at position 23 of 22-base spacer).          │
│ Throw out: No valid cut sites.                                                                          │
│ Why: Cas13d won't cut well with G at 3′; more valid sites = more druggable.                             │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ MATCHMAKER FILTER 2: DISEASE ASSOCIATION                                                                │
│ Keep: Pairs where fusion has known Associated_Disease (from disease map).                               │
│ Throw out: Unknown or missing disease.                                                                  │
│ Why: We need to know which cancer we're targeting.                                                      │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ lead_candidates.csv (pre-safety, pre-AI)                                                                │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Phase 5: Expert Agent (Safety & AI Filter)

```
   lead_candidates.csv
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ EXPERT FILTER 1: ASSOCIATED DISEASE CHECK                                                               │
│ Keep: Candidates with valid Associated_Disease (cancer type).                                           │
│ Throw out: Blanks, "Unknown", "nan".                                                                    │
│ Why: Only pursue targets we can tie to a cancer.                                                        │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ EXPERT FILTER 2: ARCHS4 EXPRESSION (SAFETY)                                                             │
│ For each fusion's parent genes:                                                                         │
│ Keep: Low in normal tissue (max ≤ threshold), present in cancer (mean ≥ min).                           │
│       Organ-specific: cancer expression ≥ 2× normal (enrichment).                                       │
│ Throw out: High expression in healthy organs (risky) or absent in cancer.                               │
│ Why: Target should be "off" in normal, "on" in tumor.                                                   │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ EXPERT FILTER 3: GEMINI AI VERDICT                                                                      │
│ Action: AI returns GO / NO-GO / HOLD for each (fusion, disease) pair.                                   │
│ Considers: Target biology, patient count, cut sites, safety risks.                                      │
│ Why: First-pass prioritization before human review.                                                     │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│ FINAL OUTPUT                                                                                            │
│ lead_candidates_filtered.csv  │  dashboard.html                                                         │
│ Type VI Cas13d enzyme–fusion pairs (cancer-specific, safety-validated)                                  │
└─────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Single-Page Wireframe Summary

```
╔══════════════════════════════════════════════════════════════════════════════════════════════════════════╗
║  COLLATERAL BIO: NCBI → TYPE VI CAS13D CANDIDATES                                                        ║
╚══════════════════════════════════════════════════════════════════════════════════════════════════════════╝

  ┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐
  │    NCBI      │     │   ENZYME     │     │   TARGET     │     │  MATCHMAKER  │     │   EXPERT     │
  │   SCRAPE     │────▶│   FILTERS    │────▶│   FILTERS    │────▶│   FILTERS    │────▶│   FILTERS    │
  └──────────────┘     └──────────────┘     └──────────────┘     └──────────────┘     └──────────────┘
       │                      │                     │                     │                     │
       ▼                      ▼                     ▼                     ▼                     ▼
  Search term          Size 600–1400         Specificity            PFS cut sites        ARCHS4 safety
  match "Cas13d"       HEPN 2–3 motifs       ≤3 cancer types        Disease known       Absent normal
  Env keywords         CRISPR array          Fusion→organ           Rank by score       Present cancer
  BioProject elink     ESM-2 sim >0.75       mapping                                     Gemini GO/HOLD

  ═══════════════════════════════════════════════════════════════════════════════════════════════════════
  FINAL: lead_candidates_filtered.csv  │  Type VI Cas13d enzyme–fusion pairs (cancer-specific)
  ═══════════════════════════════════════════════════════════════════════════════════════════════════════
```

---

## Filter Summary Table (Quick Reference)

| Stage | Filter | Keep | Throw out |
|-------|--------|------|-----------|
| NCBI | Search | Matches query / env keywords | Rest of DB |
| Enzyme | Size | 600–1400 aa | Too short/long |
| Enzyme | HEPN | 2–3 motifs, 100–1200 aa apart | Missing or 4+ domains |
| Enzyme | full_orf_checks / structure filter | ≥15 aa after last HEPN, 2–3 HEPN | C-terminal fragments, wrong HEPN count |
| Enzyme | CRISPR | Contigs with repeats | No CRISPR context |
| Enzyme | ESM-2 | Similarity > 0.75 | Doesn't look like Cas13d |
| Target | Specificity | ≤3 cancer types | Promiscuous fusions |
| Matchmaker | PFS | Valid cut sites | No cut sites |
| Matchmaker | Disease | Known cancer | Unknown disease |
| Expert | ARCHS4 | Low normal, high cancer | Unsafe profile |
| Expert | Gemini | GO / HOLD | NO-GO |
