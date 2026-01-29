import os
from Bio import Entrez, SeqIO
from datetime import datetime

Entrez.email = "founder@collateralbio.com"

class EnzymeMiner:
    def __init__(self, output_dir="data/raw_sequences"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
    
    def search_and_fetch(self, query="Cas13d", max_results=50):
        print(f"[*] Searching NCBI for: '{query}'...")
        try:
            handle = Entrez.esearch(db="protein", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            if not id_list:
                print("[!] No candidates found.")
                return []

            print(f"[*] Fetching {len(id_list)} sequences...")
            handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
            sequences = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"{self.output_dir}/search_{timestamp}.fasta"
            SeqIO.write(sequences, filename, "fasta")
            
            print(f"[SUCCESS] Saved {len(sequences)} candidates to {filename}")
            return filename
        except Exception as e:
            print(f"[!] Mining Error: {e}")
            return None

if __name__ == "__main__":
    EnzymeMiner().search_and_fetch("Cas13d NOT synthetic")