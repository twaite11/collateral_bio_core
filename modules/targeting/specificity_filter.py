import pandas as pd
import os

class SpecificityFilter:
    def __init__(self, matrix_file, output_file):
        self.matrix_file = matrix_file
        self.output_file = output_file

    def run_filter(self, max_tissue_types=3):
        print(f"[*] Loading Disease Matrix: {self.matrix_file}...")
        
        if not os.path.exists(self.matrix_file):
            print(f"[!] Error: File {self.matrix_file} not found.")
            return

        # Load Matrix (Rows = Cancers, Cols = Fusions)
        try:
            df = pd.read_csv(self.matrix_file)
            
            # Set 'Cancer' or 'Unnamed: 0' as index so we only calculate on numeric data
            if 'Cancer' in df.columns:
                df.set_index('Cancer', inplace=True)
            elif 'Unnamed: 0' in df.columns:
                df.set_index('Unnamed: 0', inplace=True)
            
            print(f"   [+] Analyzing specificity across {len(df)} cancer types...")
            
            clean_targets = []
            
            # Iterate through each Fusion (Column)
            for fusion in df.columns:
                # Get all cancer types where this fusion count > 0
                active_cancers = df[df[fusion] > 0]
                tissue_count = len(active_cancers)
                total_patient_count = active_cancers[fusion].sum()
                
                # --- THE FILTER LOGIC ---
                # 1. Must be present in at least 1 cancer
                # 2. Must NOT be present in more than 'max_tissue_types' (Artifact/Promiscuity Check)
                
                if 0 < tissue_count <= max_tissue_types:
                    primary_cancer = active_cancers[fusion].idxmax()
                    
                    clean_targets.append({
                        'Fusion_Name': fusion,
                        'Primary_Disease': primary_cancer,
                        'Tissue_Types_Count': tissue_count,
                        'Total_Patients': total_patient_count,
                        'Specificity_Score': 1.0 / tissue_count # 1.0 = Perfect Specificity
                    })
            
            # Save results
            result_df = pd.DataFrame(clean_targets)
            
            # Sort: Most Specific -> Most Abundant
            result_df = result_df.sort_values(
                by=['Specificity_Score', 'Total_Patients'], 
                ascending=[False, False]
            )
            
            result_df.to_csv(self.output_file, index=False)
            print(f"\n[SUCCESS] Filtered down to {len(result_df)} highly specific targets.")
            print(f"   [+] Saved to: {self.output_file}")
            print("--- TOP SPECIFIC TARGETS ---")
            print(result_df.head())

        except Exception as e:
            print(f"[!] Error processing matrix: {e}")

if __name__ == "__main__":
    # CONFIGURATION
    # Input: The Matrix File (Cancer x Fusion)
    # Ensure you have renamed your uploaded file to this:
    INPUT_MATRIX = "data/matrices/disease_matrix_novel.csv" 
    
    # Output: The Clean List for the Matchmaker
    OUTPUT_LIST = "data/targets/high_specificity_targets.csv"
    
    filter_tool = SpecificityFilter(INPUT_MATRIX, OUTPUT_LIST)
    filter_tool.run_filter(max_tissue_types=3)