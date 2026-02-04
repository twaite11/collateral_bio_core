"""Quick search test - verify >= 5 datasets for processing (no fetch/mining)."""
from modules.mining.sra_scout import SRAScout

scout = SRAScout()
queries = [
    "hydrothermal vent metagenome",
    "hypersaline metagenome",
    "rumen metagenome",
]

for q in queries:
    ids = scout.search_wgs(q, max_records=20)
    ok = "PASS" if len(ids) >= 5 else "FAIL"
    print(f'{ok}: "{q}" -> {len(ids)} IDs')
    if ids:
        print(f"    Sample: {ids[:3]}")

# Final assertion
ids = scout.search_wgs("hydrothermal vent metagenome", max_records=20)
assert len(ids) >= 5, f"Expected >= 5 IDs, got {len(ids)}"
print()
print(f"[OK] Search test passed: {len(ids)} datasets ready for processing on VPS")
