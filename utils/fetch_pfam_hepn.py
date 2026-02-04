"""
Fetch Pfam HEPN (PF05168) HMM for Cas13/Type VI filtering.

Downloads Pfam-A.hmm.gz, extracts PF05168, saves to data/pfam/PF05168.hmm.
Requires network access. Uses built-in parsing (no HMMER hmmfetch needed).
"""
import gzip
import os
import sys
import urllib.request

PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "pfam")
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "PF05168.hmm")
TARGET = "PF05168"


def extract_model_from_stream(stream, target_acc: str) -> str | None:
    """Extract a single HMM model from Pfam HMM stream. Returns model text or None."""
    buffer = []
    for line in stream:
        line_str = line.decode() if isinstance(line, bytes) else line
        buffer.append(line_str)
        if line_str.strip() == "//":
            text = "".join(buffer)
            if f"ACC   {target_acc}" in text or f"ACC  {target_acc}" in text:
                return text
            buffer = []
    return None


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    if os.path.exists(OUTPUT_PATH):
        print(f"[*] {OUTPUT_PATH} already exists. Delete to re-fetch.")
        return 0

    print(f"[*] Downloading Pfam-A.hmm.gz from {PFAM_URL}...")
    print("    (This may take a few minutes; file is ~100MB)")
    try:
        req = urllib.request.Request(PFAM_URL, headers={"User-Agent": "collateral_bio/1.0"})
        with urllib.request.urlopen(req, timeout=300) as resp:
            with gzip.GzipFile(fileobj=resp) as gz:
                model = extract_model_from_stream(gz, TARGET)
    except Exception as e:
        print(f"[!] Fetch failed: {e}")
        return 1

    if not model:
        print(f"[!] Model {TARGET} not found in Pfam-A.hmm.gz")
        return 1

    with open(OUTPUT_PATH, "w") as f:
        f.write(model)
    print(f"[SUCCESS] Saved {OUTPUT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
