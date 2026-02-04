import pandas as pd
import json
import os
import time
import google.generativeai as genai
import sys
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Optional Import for Safety Checks
sys.path.append(os.getcwd())
try:
    from modules.targeting.archs4_loader import ARCHS4Loader
except ImportError:
    ARCHS4Loader = None

class ExpertAgent:
    def __init__(self, candidates_file, prompt_file):
        self.candidates_file = candidates_file
        self.prompt_file = prompt_file
        self.safety_loader = None
        self.model = None
        
        # Load API key from .env file
        api_key = os.getenv("GEMINI_API_KEY")
        if not api_key:
            print("Warning: GEMINI_API_KEY not found in environment variables. AI analysis will be disabled.")
        else:
            try:
                genai.configure(api_key=api_key)
                # Use gemini-flash-latest for free tier (faster and free)
                # This is the latest flash model available in the free tier
                self.model = genai.GenerativeModel('gemini-flash-latest')
                print("[OK] Gemini API configured successfully")
            except Exception as e:
                print(f"Error configuring Gemini API: {e}")
                self.model = None
        
        if ARCHS4Loader:
            h5_path = "data/expression_data/human_matrix.h5"
            if os.path.exists(h5_path):
                self.safety_loader = ARCHS4Loader(h5_path)

    def _has_associated_disease(self, disease):
        """Return True if candidate has a valid associated disease (cancer type)."""
        if not disease or not isinstance(disease, str):
            return False
        s = str(disease).strip().lower()
        return s and s != "unknown" and s != "nan"

    def load_candidates(self, max_candidates=None):
        """Load candidates from CSV, keeping only fusion genes with associated diseases.
        Fusion transcripts are cancer-specific; we focus on disease-linked targets."""
        if not os.path.exists(self.candidates_file):
            return []
        df = pd.read_csv(self.candidates_file)
        mask = df['Associated_Disease'].apply(self._has_associated_disease)
        df = df[mask]
        if max_candidates is not None:
            df = df.head(max_candidates)
        return df.to_dict('records')

    def check_safety(self, target_name, associated_disease=None):
        """Safety check. For fusion genes with associated diseases: fusion transcript is
        cancer-specific (absent in healthy tissue), so we do not penalize parent gene
        expression in healthy organs."""
        if self._has_associated_disease(associated_disease):
            return 100, ["Fusion cancer-specific (absent in healthy tissue)"]

        if not self.safety_loader:
            return 50, ["Safety Check Skipped"]

        parts = target_name.replace("--", "-").split("-")
        genes = [p for p in parts if len(p) > 2]
        risks = []

        for gene in genes:
            df = self.safety_loader.get_gene_expression(gene)
            if df is not None:
                summary = df.groupby('Tissue')['Expression'].mean()
                for organ in ['brain', 'heart', 'liver', 'lung', 'kidney']:
                    matches = summary[summary.index.str.contains(organ, case=False, na=False)]
                    if not matches.empty and matches.mean() > 100:
                        risks.append(f"{gene} in {organ}")

        score = max(0, 100 - len(risks) * 25)
        return (score, risks) if risks else (100, ["Clean Profile"])

    def load_expert_prompt(self):
        """Load the expert persona prompt from file"""
        if not os.path.exists(self.prompt_file):
            print(f"Warning: Prompt file not found: {self.prompt_file}")
            return ""
        try:
            with open(self.prompt_file, 'r', encoding='utf-8') as f:
                prompt_content = f.read()
                if prompt_content:
                    print(f"[OK] Loaded expert persona from {self.prompt_file}")
                return prompt_content
        except Exception as e:
            print(f"Warning: Could not load prompt file: {e}")
            return ""
    
    def analyze_ai(self, candidate, risks):
        """Analyze a candidate using Gemini AI"""
        if not self.model:
            return {
                "Verdict": "N/A",
                "Rationale": "AI model not configured",
                "Strategy": "N/A"
            }
        
        # Load expert persona prompt
        expert_prompt = self.load_expert_prompt()
        
        # Build the analysis prompt
        # Include Valid Cut Sites if available
        valid_cut_sites = candidate.get('Valid_Cut_Sites', candidate.get('Cut_Sites', 'N/A'))
        
        prompt = f"""{expert_prompt}

TARGET DATA:
- Target Fusion: {candidate.get('Target_Fusion', 'Unknown')}
- Associated Disease: {candidate.get('Associated_Disease', 'Unknown')}
- Patient Count: {candidate.get('Patient_Count', 0)}
- Valid Cas13d Cut Sites (Score): {valid_cut_sites}
- Safety Risks: {', '.join(risks) if risks else 'None identified'}

Please analyze this target and return your evaluation as JSON with the following structure:
{{
    "Verdict": "GO" | "NO-GO" | "HOLD",
    "Scientific_Rationale": "Brief explanation",
    "Screening_Strategy": "Recommended approach",
    "Risk_Factors": "Key concerns or advantages"
}}
"""
        
        try:
            # Use modern Gemini API with structured output
            generation_config = {
                "response_mime_type": "application/json",
                "temperature": 0.7,
            }
            
            response = self.model.generate_content(
                prompt,
                generation_config=generation_config
            )
            
            # Parse JSON response
            if response.text:
                return json.loads(response.text)
            else:
                print(f"Warning: Empty response from Gemini API")
                return {}
                
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON response: {e}")
            print(f"Response text: {response.text if 'response' in locals() else 'N/A'}")
            return {}
        except Exception as e:
            error_msg = str(e)
            if "429" in error_msg or "quota" in error_msg.lower() or "rate limit" in error_msg.lower():
                print(f"  Rate limit hit. Waiting 60s then retrying once...")
                time.sleep(60)
                try:
                    response = self.model.generate_content(prompt, generation_config={
                        "response_mime_type": "application/json", "temperature": 0.7,
                    })
                    if response.text:
                        return json.loads(response.text)
                except Exception:
                    print(f"  Retry failed. Skipping AI for this candidate.")
                return {}
            else:
                print(f"Error calling Gemini API: {e}")
            return {}

    def run(self, max_candidates=None):
        """Run analysis on candidates. Deduplicates by fusion+disease for API calls:
        AI verdict is identical per target, so we call once per unique pair and map to all enzyme variants."""
        print("--- Collateral Bio: AI Analysis ---")
        all_candidates = self.load_candidates(max_candidates=max_candidates)

        if not all_candidates:
            print("[!] No candidates found. Ensure lead_candidates.csv has fusion genes with Associated_Disease.")
            return

        print(f"[*] Focus: fusion genes with associated diseases (cancer-specific, absent in healthy tissue)")

        # Group by unique (Target_Fusion, Associated_Disease) - AI analysis is identical per target
        groups = {}
        for c in all_candidates:
            key = (c.get('Target_Fusion', 'Unknown'), c.get('Associated_Disease', 'Unknown'))
            groups.setdefault(key, []).append(c)

        delay_sec = float(os.getenv("API_DELAY_SECONDS", "15"))
        print(f"[*] {len(all_candidates)} candidates, {len(groups)} unique fusion-disease pairs (1 API call each)")
        print(f"[*] Delay between API calls: {delay_sec}s to avoid rate limits")
        print()

        ai_cache = {}  # (fusion, disease) -> ai_res
        results = []

        for i, ((target, disease), members) in enumerate(sorted(groups.items(), key=lambda x: (x[0][1], x[0][0])), 1):
            cand = members[0]  # representative for API call
            print(f"[{i}/{len(groups)}] Analyzing {target} ({disease})... [covers {len(members)} enzyme variants]")

            safety_score, risks = self.check_safety(target, associated_disease=disease)
            ai_res = self.analyze_ai(cand, risks)
            ai_cache[(target, disease)] = (ai_res, safety_score, risks)

            for c in members:
                ar, ss, rk = ai_cache[(target, disease)]
                results.append({
                    'Target_Fusion': target,
                    'Associated_Disease': c.get('Associated_Disease', 'Unknown'),
                    'Patient_Count': c.get('Patient_Count', 0),
                    'Enzyme_Variant': c.get('Enzyme_Variant', ''),
                    'Valid_Cut_Sites': c.get('Valid_Cut_Sites', ''),
                    'Safety_Risks': '; '.join(rk) if rk else 'Clean',
                    'AI_Verdict': ar.get('Verdict', 'N/A'),
                    'AI_Rationale': ar.get('Scientific_Rationale', ar.get('Rationale', 'Analyzing...')),
                    'Screening_Strategy': ar.get('Screening_Strategy', 'N/A'),
                    'Risk_Factors': ar.get('Risk_Factors', 'N/A'),
                    'Safety_Score': ss
                })

            if i < len(groups):
                time.sleep(delay_sec)

        # Filtered cancer list: unique fusion + disease pairs
        fusion_disease = {(r['Target_Fusion'], r['Associated_Disease']) for r in results}
        cancer_list = sorted(fusion_disease, key=lambda x: (x[1], x[0]))

        print("\n--- Filtered Cancer List (Fusion + Disease) ---")
        for fusion, disease in cancer_list:
            print(f"  {fusion}  ->  {disease}")
        print(f"  ({len(cancer_list)} unique fusion-disease pairs)")

        # Write lead_candidates_filtered.csv
        out_path = "lead_candidates_filtered.csv"
        out_df = pd.DataFrame(results)
        out_df.to_csv(out_path, index=False)
        print(f"\n[SUCCESS] Analysis Complete. Wrote {out_path} ({len(results)} candidates).")
        return results

if __name__ == "__main__":
    agent = ExpertAgent("lead_candidates.csv", "prompts/expert_persona.txt")
    # max_candidates=None = all; set e.g. max_candidates=10 to test on a subset
    agent.run(max_candidates=None)