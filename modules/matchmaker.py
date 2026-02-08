import pandas as pd
from Bio import SeqIO
import random
import os
import sys

class MasterMatchmaker:
    def __init__(self, enzyme_fasta, target_file, disease_matrix_file=None):
        self.enzyme_fasta = enzyme_fasta
        self.target_file = target_file
        self.disease_matrix_file = disease_matrix_file
        self.results = []
        self.disease_map = {}

    def _load_enzymes(self):
        enzymes = []
        if not self.enzyme_fasta or not os.path.exists(self.enzyme_fasta):
            print(f"[!] Warning: Enzyme file '{self.enzyme_fasta}' not found.")
            print("    -> Using 'Mock Enzymes' for simulation.")
            return [
                {'id': 'Cas13d_Yellowstone_V1', 'pfs_rule': 'NOT_G'},
                {'id': 'Cas13d_DeepSea_V4', 'pfs_rule': 'NOT_G'}
            ]
            
        try:
            for record in SeqIO.parse(self.enzyme_fasta, "fasta"):
                enzymes.append({'id': record.id, 'pfs_rule': 'NOT_G'})
        except Exception as e:
            print(f"[!] Error parsing enzyme FASTA: {e}")
            return []

        if not enzymes:
            return [{'id': 'Cas13d_Mock_V1', 'pfs_rule': 'NOT_G'}]
        return enzymes

    def _load_disease_map(self):
        """Loads Matrix CSV to map Fusions -> Diseases (Only needed for raw files)."""
        if not self.disease_matrix_file or not os.path.exists(self.disease_matrix_file):
            return

        print(f"[*] Loading Disease Map from {self.disease_matrix_file}...")
        try:
            df = pd.read_csv(self.disease_matrix_file)
            # Set index to Cancer Type (first col)
            df.set_index(df.columns[0], inplace=True)
            
            for fusion in df.columns:
                associated_cancers = df.index[df[fusion] > 0].tolist()
                if associated_cancers:
                    self.disease_map[fusion] = ", ".join(associated_cancers)
        except Exception as e:
            print(f"[!] Error loading disease map: {e}")

    def _load_targets(self):
        print(f"[*] Loading Targets from {self.target_file}...")
        
        if not os.path.exists(self.target_file):
            print(f"[!] Critical: Target file '{self.target_file}' not found.")
            return pd.DataFrame()
        
        try:
            df = pd.read_csv(self.target_file)
            
            # 1. Handle "fusionsss" typo from ChimerDB raw files
            if 'fusionsss' in df.columns:
                df.rename(columns={'fusionsss': 'Fusion_Name'}, inplace=True)
                
            # 2. Normalize "Total_Patients" vs "count"
            if 'Total_Patients' in df.columns:
                df.rename(columns={'Total_Patients': 'count'}, inplace=True)

            print(f"   [+] Loaded {len(df)} targets.")
            return df.drop_duplicates(subset=['Fusion_Name'])
            
        except Exception as e:
            print(f"[!] Error reading file: {e}")
            return pd.DataFrame()

    def _get_target_sequence(self, fusion_name):
        # SIMULATION: 60nt Junction
        bases = ['A', 'C', 'T', 'G']
        return "".join(random.choice(bases) for _ in range(60))

    def run_matching(self):
        enzymes = self._load_enzymes()
        targets_df = self._load_targets()
        
        # Only load map if we need it
        if 'Primary_Disease' not in targets_df.columns:
            self._load_disease_map()

        if targets_df.empty or not enzymes:
            print("[!] Aborting: Missing Data.")
            return

        print(f"[*] Screening {len(enzymes)} Enzymes against top 50 Targets...")
        subset = targets_df.head(50)

        for index, row in subset.iterrows():
            target_name = row['Fusion_Name']
            patient_count = row.get('count', 0)
            
            # SMART DISEASE LOOKUP: 
            # 1. Check if Specificity Filter already found the disease
            # 2. Else, check the loaded map
            # 3. Else, Unknown
            if 'Primary_Disease' in row:
                disease = row['Primary_Disease']
            else:
                disease = self.disease_map.get(target_name, "Unknown")
            
            rna_seq = self._get_target_sequence(target_name)
            
            for enzyme in enzymes:
                score, cut_sites = self._calculate_cut_efficiency(enzyme, rna_seq)
                
                if score > 0:
                    self.results.append({
                        'Target_Fusion': target_name,
                        'Associated_Disease': disease,
                        'Patient_Count': patient_count,
                        'Enzyme_Variant': enzyme['id'],
                        'Valid_Cut_Sites': score,
                        'Junction_Sequence_Sim': rna_seq[:15] + "..."
                    })

    def _calculate_cut_efficiency(self, enzyme, rna_seq):
        valid_spacers = []
        spacer_len = 22
        for i in range(len(rna_seq) - spacer_len):
            spacer = rna_seq[i : i+spacer_len]
            # PFS Rule: No G at 3' end
            if rna_seq[i+spacer_len] != 'G':
                valid_spacers.append(spacer)
        return len(valid_spacers), valid_spacers

    def save_leads(self, filename="lead_candidates.csv"):
        if not self.results:
            print("[!] No matches found.")
            return
        
        df = pd.DataFrame(self.results)
        df = df.sort_values(by=['Associated_Disease', 'Valid_Cut_Sites'], ascending=[True, False])
        df.to_csv(filename, index=False)
        print(f"\n[SUCCESS] Generated '{filename}' with {len(df)} candidates.")
        print(df[['Target_Fusion', 'Associated_Disease', 'Enzyme_Variant', 'Valid_Cut_Sites']].head())

def _find_latest_enzyme_file(directory="data/raw_sequences"):
    """Find the most recent enzyme FASTA file"""
    if not os.path.exists(directory):
        return None
    
    fasta_files = [f for f in os.listdir(directory) if f.startswith("search_") and f.endswith(".fasta")]
    if not fasta_files:
        return None
    
    # Sort by filename (timestamp-based) to get latest
    fasta_files.sort(reverse=True)
    return os.path.join(directory, fasta_files[0])


def _find_enzyme_file():
    """Prefer fam_fasta (family-grouped), then search_*.fasta, else None."""
    # 1. Family-grouped enzymes (from family_grouper)
    fam_path = "data/mined_sequences/fam_fasta.fasta"
    if os.path.exists(fam_path):
        return fam_path
    # 2. Raw mined sequences
    return _find_latest_enzyme_file()


if __name__ == "__main__":
    # --- CONFIGURATION ---
    # Prefer fam_fasta.fasta (family-grouped), else latest search_*.fasta
    ENZYME_FILE = _find_enzyme_file()
    if not ENZYME_FILE:
        print("[!] Warning: No fam_fasta.fasta (data/mined_sequences/) or search_*.fasta found. Using mock enzymes.")
        ENZYME_FILE = None
    else:
        print(f"[*] Using enzyme file: {ENZYME_FILE}")
    
    # CASE 1: Run on Specificity Filtered Data (Recommended)
    # Check if filtered targets exist, otherwise use raw fusion files
    if os.path.exists("data/targets/high_specificity_targets.csv"):
        TARGET_FILE = "data/targets/high_specificity_targets.csv"
        DISEASE_FILE = None  # Not needed, specificity filter provides disease
        print("[*] Using specificity-filtered targets")
    else:
        # CASE 2: Run on Raw Validation Data
        TARGET_FILE = "data/targets/known_fusions.csv"
        DISEASE_FILE = "data/matrices/KB_and_Pub_Recur_per_cancer.csv"
        if not os.path.exists(DISEASE_FILE):
            DISEASE_FILE = "data/matrices/disease_matrix_known.csv"
        print("[*] Using raw fusion files")
    
    matcher = MasterMatchmaker(ENZYME_FILE, TARGET_FILE, DISEASE_FILE)
    matcher.run_matching()
    matcher.save_leads()