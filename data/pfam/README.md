# Pfam HEPN HMM (PF05168)

This directory holds the Pfam HEPN domain HMM used by the Cas13/Type VI post-filter.

## Setup

Run once to fetch the HMM:

```bash
python utils/fetch_pfam_hepn.py
```

This downloads Pfam-A.hmm.gz, extracts PF05168, and saves `PF05168.hmm`.

## Manual Setup (if fetch script fails)

If you have HMMER installed:

```bash
# Download Pfam (requires wget/curl)
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz

# Extract PF05168
hmmfetch Pfam-A.hmm PF05168 > data/pfam/PF05168.hmm
```
