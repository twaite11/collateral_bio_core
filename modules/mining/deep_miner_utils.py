import torch
import re
from transformers import AutoTokenizer, AutoModel
from Bio.Seq import Seq

class DeepEngine:
    def __init__(self, model_name="facebook/esm2_t12_35M_UR50D"):
        """
        Initializes the ESM-2 Model. 
        
        UPDATED FOR RUNPOD (GPU):
        - Upgraded model to 35M parameters (better accuracy).
        - Automatically moves computation to CUDA (GPU).
        """
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        print(f"[*] Loading Deep Learning Model: {model_name} on {self.device.upper()}...")
        
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModel.from_pretrained(model_name).to(self.device)
        self.model.eval() 
        
        # --- THE GOLDEN REFERENCE (Cas13d) ---
        self.ref_seq_segment = "RHYLDEIIEQISEFSKRVILADANLDKVLSAYNKHRDKPIREQAENIIHLFTLTNLGAPAAFKYFDTTIDRKRYTSTKEVLDATLIHQSITGLYETRIDLSQLGGD"
        
        # Calculate reference vector once on GPU
        self.ref_vector = self._get_embedding(self.ref_seq_segment)
        print(f"[+] Reference Vector Calculated on {self.device}.")

    def _get_embedding(self, sequence):
        """Converts Amino Acid sequence into Vector on GPU."""
        # Move inputs to the same device as model
        inputs = self.tokenizer(sequence, return_tensors="pt").to(self.device)
        
        with torch.no_grad():
            outputs = self.model(**inputs)
            
        # Mean pooling
        return outputs.last_hidden_state.mean(dim=1)

    def score_candidate(self, candidate_seq):
        """Returns Similarity Score (0.0 to 1.0)."""
        if len(candidate_seq) < 300: return 0.0
        process_seq = candidate_seq[:1000] # Truncate for speed
        
        try:
            cand_vector = self._get_embedding(process_seq)
            cosine_sim = torch.nn.functional.cosine_similarity(self.ref_vector, cand_vector)
            return cosine_sim.item()
        except Exception:
            return 0.0

class NeighborhoodWatch:
    """Finds CRISPR Arrays. Uses multiple chunk sizes (24-32bp) and accepts 2+ repeats."""

    def has_crispr_array(self, dna_sequence):
        seq_str = str(dna_sequence)
        length = len(seq_str)
        if length < 400:
            return False

        min_repeats = 3
        if length < 1500:
            min_repeats = 2

        for chunk_size in (24, 28, 32):
            seen = {}
            for i in range(0, length - chunk_size):
                chunk = seq_str[i : i + chunk_size]
                if chunk in seen:
                    seen[chunk] += 1
                    if seen[chunk] >= min_repeats:
                        return True
                else:
                    seen[chunk] = 1
        return False