# Pfam HEPN (PF05168)

This directory can hold the Pfam HEPN domain HMM (PF05168) for optional external use. The main pipeline uses sequence-based HEPN motif checks (2–3× R.{4,6}H) in the structure filter (`bi_lobed_hepn_check.py`), not the Pfam HMM.

## Manual setup (optional)

If you need PF05168.hmm for other tools:

```bash
# With HMMER installed:
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmfetch Pfam-A.hmm PF05168 > data/pfam/PF05168.hmm
```
