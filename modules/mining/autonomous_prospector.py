import time
import json
import os
import re
import sqlite3
import requests
import logging
import random
from datetime import datetime
import sys
from Bio import SeqIO, Entrez

# Ensure imports work for VPS environment
sys.path.append(os.getcwd())
try:
    from modules.mining.sra_scout import SRAScout
    from modules.mining.deep_miner_utils import DeepEngine, NeighborhoodWatch
except ImportError:
    print("[!] Critical: Modules missing. Ensure 'sra_scout.py' and 'deep_miner_utils.py' exist.")
    sys.exit(1)

# --- CONFIGURATION ---
DB_PATH = "data/prospector.db"
LOG_FILE = "prospector.log"
Entrez.email = "founder@senarybio.com" 

# Setup Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_FILE), logging.StreamHandler()]
)

class LLMClient:
    """Abstracts the AI backend. Optimized for Local Llama-3-Bio."""
    def __init__(self):
        self.provider = os.getenv("LLM_PROVIDER", "local").lower() 
        self.local_url = os.getenv("LLM_LOCAL_URL", "http://localhost:11434/api/chat")
        self.model_name = os.getenv("LLM_MODEL", "llama3-bio") 
        logging.info(f"Initialized LLM Backend: {self.provider} | Model: {self.model_name}")

    def generate(self, system_prompt, user_prompt):
        try:
            payload = {
                "model": self.model_name,
                "messages": [
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                "temperature": 0.8, # Slightly higher temp for more creative search terms
                "stream": False
            }
            response = requests.post(self.local_url, json=payload, timeout=120)
            
            if response.status_code == 200:
                return response.json()['message']['content'].strip()
            else:
                logging.error(f"LLM Error {response.status_code}: {response.text}")
                return ""
        except Exception as e:
            logging.error(f"Local LLM Connection Failed: {e}")
            return ""

class AutonomousProspector:
    def __init__(self):
        self.llm = LLMClient()
        self.scout = SRAScout()
        
        logging.info("Initializing Deep Learning Engines (ESM-2 & CRISPR-Context)...")
        self.deep_engine = DeepEngine()       
        self.context_engine = NeighborhoodWatch() 
        
        self._setup_db()
        self.failure_streak = 0

    def _setup_db(self):
        os.makedirs("data", exist_ok=True)
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute('''CREATE TABLE IF NOT EXISTS history 
                     (id INTEGER PRIMARY KEY, query TEXT, hits INTEGER, strategy TEXT, timestamp DATETIME)''')
        c.execute('''CREATE TABLE IF NOT EXISTS visited_ids 
                     (nucleotide_id TEXT PRIMARY KEY, timestamp DATETIME)''')
        conn.commit()
        conn.close()

    def get_visited_ids(self):
        """Return set of nucleotide IDs already mined to avoid re-processing."""
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("SELECT nucleotide_id FROM visited_ids")
        ids = {row[0] for row in c.fetchall()}
        conn.close()
        return ids

    def mark_ids_visited(self, id_list):
        """Persist nucleotide IDs as visited after processing."""
        if not id_list:
            return
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        now = datetime.now()
        c.executemany(
            "INSERT OR IGNORE INTO visited_ids (nucleotide_id, timestamp) VALUES (?, ?)",
            [(str(uid), now) for uid in id_list]
        )
        conn.commit()
        conn.close()

    def log_history(self, query, hits, strategy):
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("INSERT INTO history (query, hits, strategy, timestamp) VALUES (?, ?, ?, ?)",
                  (query, hits, strategy, datetime.now()))
        conn.commit()
        conn.close()

    def get_recent_history(self, limit=5):
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("SELECT query, hits, strategy FROM history ORDER BY id DESC LIMIT ?", (limit,))
        data = c.fetchall()
        conn.close()
        return "\n".join([f"- Strategy: {row[2]} | Query: '{row[0]}' | Hits: {row[1]}" for row in data])

    def fetch_metadata(self, id_list):
        if not id_list: return {}
        try:
            handle = Entrez.esummary(db="nucleotide", id=",".join(str(x) for x in id_list))
            records = Entrez.read(handle)
            handle.close()
            return {str(r['Id']): r.get('Title', r.get('Description', 'Unknown')) for r in records}
        except Exception as e:
            logging.error(f"Metadata fetch failed: {e}")
            return {}

    def semantic_filter(self, metadata_dict):
        if not metadata_dict: return []

        deep_mine_max = int(os.getenv("DEEP_MINE_MAX", "15"))
        meta_keys = list(metadata_dict.keys())
        descriptions = "\n".join([f"ID {k}: {v}" for k, v in list(metadata_dict.items())[:20]])

        system = "You are an expert computational biologist specializing in metagenomics."
        prompt = f"""
        We are mining for **Type VI CRISPR-Cas13d** systems.
        Review these dataset descriptions. Select the Top {min(deep_mine_max, 20)} that are most likely to contain **uncultured CRISPR loci**, **viral defense islands**, or **diverse microbial communities**.

        Look for high-value terms: 'microbial mat', 'biofilm', 'hot spring', 'hypersaline', 'uncultured', 'phage defense', 'rumen', 'permafrost', 'deep sea'.
        Avoid: 'human', 'clinical', 'mitochondrion', 'chloroplast'.

        DATASETS:
        {descriptions}

        OUTPUT FORMAT:
        Return ONLY a JSON list of IDs. Do not write anything else.
        Example: ["12345", "67890"]
        """

        response = self.llm.generate(system, prompt)
        selected_ids = []

        try:
            if response and '[' in response and ']' in response:
                clean_json = response[response.find('['):response.rfind(']')+1]
                parsed = json.loads(clean_json)
                selected_ids = [str(i) for i in parsed if str(i) in metadata_dict]
        except (json.JSONDecodeError, TypeError):
            pass

        if not selected_ids and response:
            raw_ids = re.findall(r'\d+', response)
            selected_ids = [uid for uid in raw_ids if uid in metadata_dict]
            if selected_ids:
                logging.info("Semantic filter: extracted IDs via regex fallback.")

        if len(selected_ids) < 3:
            fallback_count = min(deep_mine_max, len(meta_keys))
            selected_ids = meta_keys[:fallback_count]
            logging.warning(f"Semantic filter returned {len(selected_ids) if selected_ids else 0} valid IDs. Using first {fallback_count} from metadata.")

        return selected_ids[:deep_mine_max]

    def formulate_strategy(self):
        """The 'Ralph' Loop: Adapts strategy based on success/failure"""
        history = self.get_recent_history()
        
        system = "You are the Autonomous Director of Discovery. Plan the next NCBI mining operation."
        prompt = f"""
        GOAL: Discover novel Cas13d enzymes in NCBI Whole Genome Shotgun (WGS) data.
        
        STATUS:
        - Failure Streak: {self.failure_streak}
        - Recent History:
        {history}
        
        TASK:
        Generate the next specific NCBI Query.
        1. If Failure Streak > 2: PIVOT to a new environment type.
        2. If succeeding: REFINE the current query.
        3. DIVERSITY IS KEY. Do not stick to just hot springs.
        
        ENVIRONMENT IDEAS:
        - Hypersaline (Salt lakes, brines)
        - Digestive (Termite gut, Rumen)
        - Marine (Deep sea vents, Sediment)
        - Soil (Permafrost, Desert varnish)
        - Symbionts (Sponge symbiont, Coral microbiome)
        
        OUTPUT FORMAT:
        JSON with keys: "strategy_name", "query".
        Return ONLY environmental keywords for "query" - no filters like wgs[Prop].
        Example: {{ "strategy_name": "Pivot to Hypersaline", "query": "hypersaline microbial mat metagenome" }}
        """
        
        try:
            response = self.llm.generate(system, prompt)
            clean_json = response[response.find('{'):response.rfind('}')+1]
            data = json.loads(clean_json)
            return data.get('query'), data.get('strategy_name')
        except:
            # BROADENED FALLBACK LIST (keywords only, no filters)
            fallbacks = [
                ("hydrothermal vent metagenome", "Fallback: Thermal"),
                ("hypersaline metagenome", "Fallback: Saline"),
                ("permafrost metagenome", "Fallback: Cryo"),
                ("rumen metagenome", "Fallback: Gut"),
                ("deep sea sediment metagenome", "Fallback: Marine")
            ]
            query, strategy = random.choice(fallbacks)
            return query, strategy

    def deep_mine(self, id_list):
        require_crispr = os.getenv("REQUIRE_CRISPR", "1").lower() in ("1", "true", "yes")
        esm_threshold = float(os.getenv("ESM_THRESHOLD", "0.75"))
        hits = []
        for ncbi_id in id_list:
            logging.info(f"   -> Scanning ID {ncbi_id}...")
            contigs_scanned = 0
            contigs_with_crispr = 0
            orfs_scored = 0
            best_score = 0.0
            try:
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
                records = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                for record in records:
                    dna_seq = record.seq
                    contigs_scanned += 1
                    has_crispr = self.context_engine.has_crispr_array(dna_seq)
                    if has_crispr:
                        contigs_with_crispr += 1
                    if require_crispr and not has_crispr:
                        continue
                    if has_crispr:
                        logging.info(f"      [!] CRISPR Array detected in {ncbi_id}. Analyzing proteins...")

                    frames = [dna_seq[i:].translate(to_stop=False) for i in range(3)]
                    frames += [dna_seq.reverse_complement()[i:].translate(to_stop=False) for i in range(3)]
                    for frame in frames:
                        orfs = str(frame).split("*")
                        for orf in orfs:
                            if 800 < len(orf) < 1100:
                                score = self.deep_engine.score_candidate(orf)
                                orfs_scored += 1
                                if score > best_score:
                                    best_score = score
                                if score > esm_threshold:
                                    logging.info(f"      [***] DISCOVERY: Novel Enzyme (Score: {score:.4f})")
                                    hits.append((f"{ncbi_id}_ORF", orf, score))

                if contigs_scanned > 0:
                    logging.info(f"      ID {ncbi_id}: {contigs_scanned} contigs, {contigs_with_crispr} with CRISPR, {orfs_scored} ORFs scored, best={best_score:.3f}")
                if orfs_scored > 0 and best_score < esm_threshold:
                    logging.info(f"      (Best score {best_score:.3f} below threshold {esm_threshold})")
            except Exception as e:
                logging.error(f"Error mining {ncbi_id}: {e}")

        return hits

    def run_loop(self):
        require_crispr = os.getenv("REQUIRE_CRISPR", "1").lower() in ("1", "true", "yes")
        esm_threshold = float(os.getenv("ESM_THRESHOLD", "0.75"))
        deep_mine_max = int(os.getenv("DEEP_MINE_MAX", "15"))
        logging.info(f"--- Senary Bio: Deep Prospector (Model: {self.llm.model_name}) ---")
        logging.info(f"Config: DEEP_MINE_MAX={deep_mine_max}, ESM_THRESHOLD={esm_threshold}, REQUIRE_CRISPR={require_crispr}")

        while True:
            query, strategy = self.formulate_strategy()
            logging.info(f"PLAN: {strategy} | Query: {query}")
            
            raw_ids = self.scout.search_wgs(query, max_records=50)

            # Exclude already-visited IDs to avoid endless re-processing
            visited = self.get_visited_ids()
            ids = [uid for uid in raw_ids if uid not in visited]
            if visited:
                logging.info(f"Skipping {len(visited)} previously visited IDs. {len(ids)} new candidates.")

            if not ids:
                logging.warning("All IDs already visited or none found. Pivoting next iteration.")
                self.failure_streak += 1
                time.sleep(5)
                continue

            logging.info(f"Analyzing {len(ids)} candidates with {self.llm.model_name}...")
            metadata = self.fetch_metadata(ids)
            best_ids = self.semantic_filter(metadata)
            
            logging.info(f"FILTER: Selected {len(best_ids)} high-probability targets.")

            if not best_ids:
                logging.warning("No targets to mine (semantic filter returned empty). Skipping iteration.")

            hits = 0
            if best_ids:
                new_discoveries = self.deep_mine(best_ids)
                if new_discoveries:
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    with open(f"data/raw_sequences/deep_hits_{timestamp}.fasta", "w") as f:
                        for name, seq, score in new_discoveries:
                            f.write(f">{name}_Score_{score:.3f}\n{seq}\n")
                    
                    logging.info(f"SUCCESS: Saved {len(new_discoveries)} Deep Learning Hits.")
                    hits = len(new_discoveries)
                    self.failure_streak = 0
                else:
                    self.failure_streak += 1
            
            self.log_history(query, hits, strategy)

            # Mark all IDs from this search as visited to avoid re-processing
            self.mark_ids_visited(raw_ids)

            logging.info("Sleeping for 30s...")
            time.sleep(30)

if __name__ == "__main__":
    prospector = AutonomousProspector()
    prospector.run_loop()