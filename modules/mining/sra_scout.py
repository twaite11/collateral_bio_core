import os
import re
from Bio import Entrez, Seq, SeqIO
from datetime import datetime

# Configure your email for NCBI
Entrez.email = "founder@senarybio.com"

class SRAScout:
    def __init__(self, output_dir="data/raw_sequences"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        
        # --- CAS13d PHYSICS ENGINE (UPDATED) ---
        # 1. HEPN Motif: Arginine, 4-6 spacers, Histidine
        # Relaxed from R.{4}H to R.{4,6}H to catch novel diversity
        self.hepn_regex = re.compile(r"R.{4,6}H")
        
        # 2. Size Constraints (Cas13d family: 600-1400 aa to capture diversity)
        self.min_size = 600
        self.max_size = 1400

    def _normalize_query(self, query):
        """Strip NCBI filter suffixes from query so caller controls filters."""
        if not query or not isinstance(query, str):
            return query or ""
        q = query.strip()
        for suffix in [" wgs[Prop]", " wgs[prop]", "[Prop]", "[prop]", "[Title]", "[title]", "[All Fields]"]:
            if q.lower().endswith(suffix.lower()):
                q = q[:-len(suffix)].strip()
        return q

    def _search_nucleotide(self, term, max_records):
        """Internal: run esearch on nucleotide with given term."""
        try:
            handle = Entrez.esearch(db="nucleotide", term=term, retmax=max_records)
            record = Entrez.read(handle)
            handle.close()
            return record.get("IdList", [])
        except Exception as e:
            print(f"[!] Nucleotide search error: {e}")
            return []

    def search_bioproject(self, environment_query, max_projects=10):
        """
        Search BioProject for metagenome studies (extremophile environments).
        Returns list of BioProject IDs.
        """
        keywords = self._normalize_query(environment_query)
        if not keywords:
            keywords = "hot spring metagenome"
        term = f'{keywords} AND "scope metagenome"[Properties]'
        print(f"[*] Scouting BioProject for: '{term}'...")
        try:
            handle = Entrez.esearch(db="bioproject", term=term, retmax=max_projects)
            record = Entrez.read(handle)
            handle.close()
            return record.get("IdList", [])
        except Exception as e:
            print(f"[!] BioProject search error: {e}")
            return []

    def bioproject_to_nucleotide_ids(self, bioproject_ids, max_per_project=20):
        """
        Use elink to get nucleotide UIDs linked to BioProject IDs.
        Returns combined list of nucleotide IDs for fetch_and_mine / deep_mine.
        """
        all_nuc_ids = []
        for bp_id in bioproject_ids:
            try:
                handle = Entrez.elink(dbfrom="bioproject", db="nucleotide", id=bp_id)
                record = Entrez.read(handle)
                handle.close()
                for link_set in record:
                    if "LinkSetDb" in link_set:
                        for link_db in link_set["LinkSetDb"]:
                            ids = [str(x) for x in link_db.get("Id", [])[:max_per_project]]
                            all_nuc_ids.extend(ids)
            except Exception as e:
                print(f"[!] elink BioProject {bp_id} error: {e}")
        seen = set()
        unique = []
        for uid in all_nuc_ids:
            if uid not in seen:
                seen.add(uid)
                unique.append(uid)
        return unique

    def search_wgs(self, environment_query="hot spring metagenome", max_records=20):
        """
        Searches NCBI Nucleotide for unannotated metagenomic contigs.
        Normalizes query (strips wgs[Prop] etc.), tries WGS first, then broader search if no hits.
        """
        keywords = self._normalize_query(environment_query)
        if not keywords:
            keywords = "hot spring metagenome"

        # Try 1: environment terms + wgs[Prop] (only once)
        full_query = f'{keywords} AND wgs[Prop]'
        print(f"[*] Scouting Nucleotide for: '{full_query}'...")
        id_list = self._search_nucleotide(full_query, max_records)

        # Try 2: if no hits, retry without wgs[Prop] for broader matching
        if len(id_list) < 5:
            fallback_query = f'{keywords}'
            print(f"[*] Few hits with wgs[Prop]. Retrying broader: '{fallback_query}'...")
            id_list = self._search_nucleotide(fallback_query, max_records)

        # Try 3: fall back to BioProject -> elink -> nucleotide
        if len(id_list) < 5:
            bp_ids = self.search_bioproject(environment_query, max_projects=10)
            if bp_ids:
                id_list = self.bioproject_to_nucleotide_ids(bp_ids, max_per_project=max_records // max(1, len(bp_ids)))
                print(f"   [+] BioProject fallback: {len(id_list)} nucleotide IDs from {len(bp_ids)} projects.")

        print(f"   [+] Found {len(id_list)} metagenomic contigs/scaffolds.")
        return id_list

    def fetch_and_mine(self, id_list):
        """
        Downloads DNA, translates to Protein (6 frames), and hunts for Cas13d.
        """
        if not id_list:
            return

        print(f"[*] Mining {len(id_list)} contigs for hidden Cas13d signals...")
        
        # Fetch DNA sequences
        # using 'fasta' rettype
        handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
        
        candidates = []
        
        for record in SeqIO.parse(handle, "fasta"):
            dna_seq = record.seq
            
            # --- 6-FRAME TRANSLATION ---
            # DNA can be read in 3 frames forward, 3 frames backward.
            # We check all of them.
            
            frames = []
            # Forward strand frames
            frames.append(dna_seq.translate(to_stop=False))
            frames.append(dna_seq[1:].translate(to_stop=False))
            frames.append(dna_seq[2:].translate(to_stop=False))
            
            # Reverse complement frames
            rc_seq = dna_seq.reverse_complement()
            frames.append(rc_seq.translate(to_stop=False))
            frames.append(rc_seq[1:].translate(to_stop=False))
            frames.append(rc_seq[2:].translate(to_stop=False))
            
            for i, protein_seq in enumerate(frames):
                # Only look at proteins in the Cas13d size range
                # We split by stop codons (*) to find Open Reading Frames (ORFs)
                orfs = str(protein_seq).split("*")
                
                for orf in orfs:
                    if self.min_size <= len(orf) <= self.max_size:
                        if self._is_cas13d_candidate(orf):
                            # FOUND ONE!
                            cand_id = f"NewCas13d_{record.id}_Frame{i}"
                            candidates.append((cand_id, orf))
                            print(f"   [!!!] DISCOVERY: Potential Cas13d found in {record.id} (Frame {i})")

        handle.close()
        return candidates

    def _is_cas13d_candidate(self, protein_seq):
        """
        Apply the HEPN logic.
        Cas13d needs 2-3 HEPN domains (RxxxxH) separated by space.
        """
        matches = [m.start() for m in self.hepn_regex.finditer(protein_seq)]
        
        if len(matches) < 2 or len(matches) > 3:
            return False
        
        # Check distance between HEPN domains (Cas13d topology)
        # Separated by 100-1200 aa to capture diverse Type VI variants
        valid_topology = False
        for i in range(len(matches)):
            for j in range(i+1, len(matches)):
                distance = matches[j] - matches[i]
                if 100 < distance < 1200:
                    valid_topology = True
                    break
        
        return valid_topology

    def save_discoveries(self, candidates):
        if not candidates:
            print("[*] No hidden enzymes found in this batch.")
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{self.output_dir}/undiscovered_cas13d_{timestamp}.fasta"
        
        with open(filename, "w") as f:
            for cand_id, seq in candidates:
                f.write(f">{cand_id}\n{seq}\n")
        
        print(f"\n[SUCCESS] Saved {len(candidates)} NEW enzymes to {filename}")
        print("These sequences are likely unannotated. You own this IP.")

if __name__ == "__main__":
    scout = SRAScout()
    
    # Example search: "Hot Spring" metagenomes often hold Cas13d ancestors
    # We look for 'WGS' entries which are raw assemblies
    ids = scout.search_wgs("hydrothermal vent metagenome", max_records=50)
    
    new_enzymes = scout.fetch_and_mine(ids)
    scout.save_discoveries(new_enzymes)