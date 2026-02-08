import pandas as pd
import os

class ExcelSplitter:
    def __init__(self, excel_path, output_dir="data", targets_dir=None, matrices_dir=None):
        self.excel_path = excel_path
        self.output_dir = output_dir
        self.targets_dir = targets_dir or output_dir
        self.matrices_dir = matrices_dir or output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.targets_dir, exist_ok=True)
        os.makedirs(self.matrices_dir, exist_ok=True)

    def run_split(self):
        print(f"[*] Opening Excel file: {self.excel_path}...")
        
        try:
            xls = pd.ExcelFile(self.excel_path)
            sheet_names = xls.sheet_names
            print(f"[*] Found sheets: {sheet_names}")
            
            rename_map = {
                'HRS_Recurrent_known': 'known_fusions.csv',
                'HRS_Recurrent_novel': 'novel_fusions.csv',
            }

            for sheet in sheet_names:
                print(f"   -> Processing sheet: '{sheet}'...")
                df = pd.read_excel(xls, sheet_name=sheet)
                
                if sheet in rename_map:
                    filename = rename_map[sheet]
                    out_dir = self.targets_dir
                    print(f"      [!] Renaming to critical file: {filename}")
                else:
                    safe_name = sheet.replace(" ", "_").replace("&", "and")
                    filename = f"{safe_name}.csv"
                    out_dir = self.matrices_dir
                
                output_path = os.path.join(out_dir, filename)
                df.to_csv(output_path, index=False)
                print(f"      [+] Saved to {output_path}")

            print("\n[SUCCESS] Split complete. Your data folder is ready.")

        except FileNotFoundError:
            print(f"[!] Error: Could not find file '{self.excel_path}'")
            print("    Make sure the Excel file is in data/source/ or provide the full path.")
        except Exception as e:
            print(f"[!] Error processing file: {e}")

if __name__ == "__main__":
    EXCEL_FILE = "data/source/Recurrent_table.xlsx"
    splitter = ExcelSplitter(
        EXCEL_FILE,
        output_dir="data",
        targets_dir="data/targets",
        matrices_dir="data/matrices",
    )
    splitter.run_split()