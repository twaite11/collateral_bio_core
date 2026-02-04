from Bio import Entrez, SeqIO
import sys
import socket

# Configure Email (REQUIRED by NCBI, but not verified)
# Use a generic placeholder if you don't want to use your own
Entrez.email = "tyco711@gmail.com" 

def check_internet():
    try:
        # Check connection to NCBI
        socket.create_connection(("www.ncbi.nlm.nih.gov", 443))
        return True
    except OSError:
        return False

def debug_fetch():
    print("[*] TEST 1: Checking Network Connection...")
    if not check_internet():
        print("[!] TEST FAILED: No internet connection to NCBI.")
        return
    print("[+] Internet Connection OK.")

    # Use a SPECIFIC, KNOWN Accession ID to test fetching directly
    # This is E. coli K-12 MG1655, a model organism genome
    known_id = "U00096.3" 
    
    print(f"[*] TEST 2: Fetching known record {known_id} directly...")
    
    try:
        # Fetch directly, skipping search
        handle = Entrez.efetch(db="nucleotide", id=known_id, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        
        if not records:
            print("[!] TEST FAILED: Downloaded file was empty.")
            return
            
        first_record = records[0]
        dna_sample = str(first_record.seq)[:50]
        print(f"[+] TEST SUCCESS: Downloaded {len(records)} contig(s).")
        print(f"    -> Sample DNA: {dna_sample}...")
        
        # Verify Translation (checking standard genetic code)
        protein_sample = str(first_record.seq.translate(to_stop=False))[:50]
        print(f"    -> Sample Protein (Frame 1): {protein_sample}...")
        
        print("\n[VERDICT] The pipeline plumbing is WORKING.")
        
    except Exception as e:
        print(f"[!] TEST FAILED: Error during fetch/parse: {e}")

if __name__ == "__main__":
    debug_fetch()