"""
Propose point mutations, score with ESM-2 (stability vs RfxCas13d/PspCas13a), and keep variants
that stay below a max identity vs reference (drift goal). Optional: ask Gemini for mutations
that maintain structural stability but might increase trans-cleavage activity.
"""
import os
import sys
import re
import json
import time
import argparse
import random
import heapq
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

sys.path.insert(0, os.getcwd())
from modules.mining.deep_miner_utils import DeepEngine

# Standard amino acids
AA = list("ACDEFGHIKLMNPQRSTVWY")


def pairwise_identity(seq1: str, seq2: str) -> float:
    """Global alignment identity (fraction of aligned positions that match)."""
    # #region agent log
    log_path = os.path.join(os.getcwd(), ".cursor", "debug.log")
    try:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(json.dumps({"id": "log_pairwise_entry", "timestamp": time.time() * 1000, "location": "mutate_for_drift.py:24", "message": "pairwise_identity entry", "data": {"seq1_len": len(seq1), "seq2_len": len(seq2), "rough_identity": sum(1 for a, b in zip(seq1[:100], seq2[:100]) if a == b) / max(100, 1) if seq1 and seq2 else 0.0}, "runId": "debug", "hypothesisId": "A,B"}) + "\n")
    except: pass
    # #endregion
    aligner = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"))
    # Limit to 1 alignment to prevent overflow when counting optimal alignments
    aligner.max_alignments = 1
    # #region agent log
    try:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(json.dumps({"id": "log_before_align", "timestamp": time.time() * 1000, "location": "mutate_for_drift.py:33", "message": "before aligner.align with max_alignments=1", "data": {}, "runId": "debug", "hypothesisId": "C"}) + "\n")
    except: pass
    # #endregion
    try:
        alns = aligner.align(seq1, seq2)
        # #region agent log
        try:
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(json.dumps({"id": "log_after_align", "timestamp": time.time() * 1000, "location": "mutate_for_drift.py:40", "message": "after aligner.align", "data": {"alns_type": str(type(alns))}, "runId": "debug", "hypothesisId": "C"}) + "\n")
        except: pass
        # #endregion
        # Use next(iter()) to safely get first alignment without triggering len() calculation
        a = next(iter(alns), None)
        if a is None:
            return 0.0
    except OverflowError as e:
        # #region agent log
        try:
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(json.dumps({"id": "log_overflow_caught", "timestamp": time.time() * 1000, "location": "mutate_for_drift.py:47", "message": "OverflowError caught, using fallback", "data": {"error": str(e)}, "runId": "debug", "hypothesisId": "C"}) + "\n")
        except: pass
        # #endregion
        # Fallback: use local alignment which is faster and less likely to overflow
        aligner_fallback = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"))
        aligner_fallback.max_alignments = 1
        alns_fallback = aligner_fallback.align(seq1, seq2)
        a = next(iter(alns_fallback), None)
        if a is None:
            return 0.0
    matches = sum(1 for i, j in zip(a[0], a[1]) if i == j and i != "-")
    length = max(len(seq1), len(seq2))
    result = matches / length if length else 0.0
    # #region agent log
    try:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(json.dumps({"id": "log_pairwise_exit", "timestamp": time.time() * 1000, "location": "mutate_for_drift.py:58", "message": "pairwise_identity exit", "data": {"result": result, "matches": matches, "length": length}, "runId": "debug", "hypothesisId": "A,B"}) + "\n")
    except: pass
    # #endregion
    return result


def max_identity_to_refs(seq: str, ref_seqs: list) -> float:
    return max(pairwise_identity(seq, r) for r in ref_seqs) if ref_seqs else 0.0


def propose_mutations(sequence: str, num_mutants: int = 5, num_sites: int = 3):
    """Propose num_mutants variants, each with up to num_sites random AA substitutions."""
    seq = list(sequence)
    L = len(seq)
    mutants = []
    seen = set()
    for _ in range(num_mutants * 4):  # allow retries
        if len(mutants) >= num_mutants:
            break
        positions = random.sample(range(L), min(num_sites, L)) if L else []
        mut = seq[:]
        desc = []
        for pos in positions:
            wild = mut[pos]
            choices = [a for a in AA if a != wild]
            if not choices:
                continue
            mut[pos] = random.choice(choices)
            desc.append(f"{pos}{wild}{mut[pos]}")
        mut_str = "".join(mut)
        key = (tuple(positions), mut_str)
        if key not in seen:
            seen.add(key)
            mutants.append((mut_str, "_".join(desc)))
    return mutants


def suggest_trans_cleavage_mutations(seq: str, ref_name: str, prompt_path: str, max_len_preview: int = 400):
    """
    Call Gemini with trans-cleavage prompt; return list of (mut_seq, desc) from suggested mutations.
    Returns [] if no API key or parse failure.
    """
    api_key = os.getenv("GEMINI_API_KEY")
    if not api_key:
        return []
    try:
        import google.generativeai as genai
        genai.configure(api_key=api_key)
        model = genai.GenerativeModel("gemini-flash-latest")
    except Exception:
        return []
    path = Path(prompt_path)
    if not path.exists():
        return []
    prompt_text = path.read_text(encoding="utf-8")
    prompt_text = prompt_text.replace("{ref_name}", ref_name).replace("{seq_preview}", seq[:max_len_preview])
    try:
        response = model.generate_content(prompt_text)
        text = (response.text or "").strip()
        # Extract JSON array
        start = text.find("[")
        end = text.rfind("]") + 1
        if start < 0 or end <= start:
            return []
        arr = json.loads(text[start:end])
        seq_list = list(seq)
        mut_seqs = []
        for item in arr:
            pos_1 = item.get("position_1based")
            wild = item.get("wild", "")
            mutant = item.get("mutant", "")
            if pos_1 is None or not wild or not mutant or len(wild) != 1 or len(mutant) != 1:
                continue
            pos_0 = int(pos_1) - 1
            if pos_0 < 0 or pos_0 >= len(seq_list) or seq_list[pos_0] != wild:
                continue
            mut = seq_list[:]
            mut[pos_0] = mutant
            mut_str = "".join(mut)
            desc = f"{pos_0}{wild}{mutant}"
            mut_seqs.append((mut_str, desc))
        return mut_seqs
    except (json.JSONDecodeError, ValueError, KeyError, AttributeError) as e:
        return []


def main():
    parser = argparse.ArgumentParser(description="Mutate for drift: ESM-2 stability vs RfxCas13d/PspCas13a + identity < threshold. Optional: Gemini suggests trans-cleavage mutations.")
    parser.add_argument("--input", default="data/raw_sequences/deep_hits_latest.fasta", help="Input FASTA")
    parser.add_argument("--references", default="data/references/known_cas13.fasta",
                        help="FASTA of known Cas13 for identity/drift check")
    parser.add_argument("--stability-refs", default=None,
                        help="FASTA for ESM stability (RfxCas13d, PspCas13a). Default: ESM_REFERENCE_FASTA or data/references/mining_refs.fasta")
    parser.add_argument("--output", default="data/design/drift_variants.fasta", help="Output FASTA")
    parser.add_argument("--max-identity", type=float, default=0.85, help="Max identity to any reference (drift goal)")
    parser.add_argument("--min-esm-score", type=float, default=0.5, help="Min ESM-2 similarity to stability ref (RfxCas13d/PspCas13a)")
    parser.add_argument("--num-mutants-per-seq", type=int, default=5, help="Max mutants to try per input sequence")
    parser.add_argument("--num-mutations", type=int, default=3, help="Number of positions to mutate per variant (random)")
    parser.add_argument("--use-trans-cleavage-prompt", action="store_true",
                        help="Ask Gemini for mutations that maintain stability but might increase trans-cleavage; validate with ESM")
    parser.add_argument("--trans-cleavage-prompt", default="prompts/trans_cleavage_mutations.txt", help="Path to trans-cleavage prompt")
    parser.add_argument("--max-seqs-per-enzyme", type=int, default=None,
                        help="Limit to top N sequences per enzyme (by ESM score) to reduce memory usage. Default: no limit")
    args = parser.parse_args()

    input_path = Path(args.input)
    ref_path = Path(args.references)
    stability_fasta = args.stability_refs or os.getenv("ESM_REFERENCE_FASTA", "").strip() or "data/references/mining_refs.fasta"
    stability_path = Path(stability_fasta)

    if not input_path.exists():
        print(f"[!] Input not found: {input_path}")
        return 1

    ref_seqs = []
    if ref_path.exists():
        ref_seqs = [str(r.seq) for r in SeqIO.parse(ref_path, "fasta")]
    else:
        print(f"[!] No reference FASTA at {ref_path}; identity check will be 0.")

    engine = DeepEngine(reference_fasta=stability_path if stability_path.exists() else None)
    records = list(SeqIO.parse(input_path, "fasta"))
    
    # Cache for enzyme assignments (populated during batch filtering)
    enzyme_cache = {}
    
    # Filter to top N sequences per enzyme if requested
    if args.max_seqs_per_enzyme and args.max_seqs_per_enzyme > 0:
        print(f"[*] Scoring {len(records)} sequences to select top {args.max_seqs_per_enzyme} per enzyme...")
        
        # Pre-filter by length
        valid_records = [(rec, str(rec.seq)) for rec in records if len(str(rec.seq)) >= 300]
        
        if not valid_records:
            print("[!] No sequences >= 300 AA. Exiting.")
            return 1
        
        # Use batch scoring for efficiency
        sequences = [seq for _, seq in valid_records]
        batch_size = int(os.getenv("EMBED_BATCH_SIZE", "50"))
        
        try:
            import torch
            embeddings, valid_indices = engine.get_embeddings_batch(sequences, batch_size=batch_size)
            
            if embeddings is None or len(valid_indices) == 0:
                print("[!] Batch scoring failed, falling back to individual scoring")
                raise ValueError("Batch scoring failed")
            
            # Score all sequences against all reference enzymes in batch
            scores_by_enzyme = defaultdict(list)  # enzyme -> [(score, rec_idx), ...]
            
            # Normalize embeddings for cosine similarity
            embeddings_norm = embeddings / embeddings.norm(dim=1, keepdim=True).clamp(min=1e-9)
            
            for ref_idx, (ref_name, ref_vec) in enumerate(zip(engine.ref_names, engine.ref_vectors)):
                # Normalize reference vector
                # ref_vec is shape (1, embed_dim) from _get_embedding, squeeze to (embed_dim,)
                ref_vec_squeezed = ref_vec.squeeze(0) if ref_vec.dim() > 1 else ref_vec
                ref_vec_norm = ref_vec_squeezed / ref_vec_squeezed.norm().clamp(min=1e-9)
                # Expand to (1, embed_dim) for broadcasting with (batch_size, embed_dim)
                ref_vec_norm = ref_vec_norm.unsqueeze(0)
                
                # Compute cosine similarity for all sequences at once
                # cosine_similarity with dim=1: compares along embedding dimension
                # Returns shape (batch_size,) when comparing (1, embed_dim) with (batch_size, embed_dim)
                similarities = torch.nn.functional.cosine_similarity(
                    ref_vec_norm, embeddings_norm, dim=1
                )
                
                # Convert to numpy and ensure it's 1D
                similarities_np = similarities.cpu().detach().numpy()
                if similarities_np.ndim > 1:
                    similarities_np = similarities_np.flatten()
                
                # Store scores for this enzyme
                for emb_idx, sim_score in enumerate(similarities_np):
                    rec_idx = valid_indices[emb_idx]
                    scores_by_enzyme[ref_name].append((float(sim_score), rec_idx))
            
            # For each sequence, find the best enzyme match
            best_by_sequence = {}  # rec_idx -> (best_score, best_enzyme)
            for enzyme, scored_list in scores_by_enzyme.items():
                for score, rec_idx in scored_list:
                    if rec_idx not in best_by_sequence or score > best_by_sequence[rec_idx][0]:
                        best_by_sequence[rec_idx] = (score, enzyme)
            
            # Group by best enzyme and keep top N per enzyme using heaps
            enzyme_groups = defaultdict(list)  # enzyme -> [(score, rec_idx), ...]
            for rec_idx, (score, enzyme) in best_by_sequence.items():
                enzyme_groups[enzyme].append((score, rec_idx))
            
            # Keep top N per enzyme
            filtered_indices = set()
            for enzyme, scored_list in enzyme_groups.items():
                if len(scored_list) > args.max_seqs_per_enzyme:
                    # Use heap to efficiently get top N
                    top_n = heapq.nlargest(args.max_seqs_per_enzyme, scored_list, key=lambda x: x[0])
                else:
                    top_n = scored_list
                
                for score, rec_idx in top_n:
                    filtered_indices.add(rec_idx)
                
                top_scores = [f"{s:.3f}" for s, _ in top_n[:5]]
                print(f"[+] Selected {len(top_n)}/{len(scored_list)} sequences for {enzyme} (top scores: {top_scores})")
            
            # Reconstruct filtered records and cache enzyme assignments
            filtered_records_list = []
            for rec_idx in sorted(filtered_indices):
                rec, seq = valid_records[rec_idx]
                filtered_records_list.append(rec)
                # Cache the enzyme assignment from batch processing
                if rec_idx in best_by_sequence:
                    _, enzyme = best_by_sequence[rec_idx]
                    enzyme_cache[seq] = enzyme
            records = filtered_records_list
            print(f"[+] Filtered to {len(records)} total sequences across {len(enzyme_groups)} enzymes")
            
        except (ValueError, ImportError, AttributeError) as e:
            # Fallback to individual scoring if batch fails
            print(f"[*] Batch scoring not available ({e}), using individual scoring...")
            scored_by_enzyme = {}
            for rec, seq in valid_records:
                try:
                    score, enzyme = engine.score_candidate_with_ref(seq)
                    if enzyme not in scored_by_enzyme:
                        scored_by_enzyme[enzyme] = []
                    scored_by_enzyme[enzyme].append((score, rec))
                    # Cache enzyme assignment
                    enzyme_cache[seq] = enzyme
                except Exception as ex:
                    print(f"[!] Failed to score {rec.id}: {ex}")
                    continue
            
            # Sort by score (descending) and take top N per enzyme
            filtered_records = []
            for enzyme, scored_list in scored_by_enzyme.items():
                scored_list.sort(key=lambda x: x[0], reverse=True)
                top_n = scored_list[:args.max_seqs_per_enzyme]
                filtered_records.extend([rec for _, rec in top_n])
                top_scores = [f"{s:.3f}" for s, _ in top_n[:5]]
                print(f"[+] Selected {len(top_n)}/{len(scored_list)} sequences for {enzyme} (top scores: {top_scores})")
            
            records = filtered_records
            print(f"[+] Filtered to {len(records)} total sequences across {len(scored_by_enzyme)} enzymes")
    
    kept = []
    for rec in records:
        seq = str(rec.seq)
        if len(seq) < 300:
            continue
        base_identity = max_identity_to_refs(seq, ref_seqs)
        
        # Use cached enzyme from batch processing if available, otherwise score
        if seq in enzyme_cache:
            closest_ref = enzyme_cache[seq]
        else:
            _, closest_ref = engine.score_candidate_with_ref(seq)

        mutants = []
        if args.use_trans_cleavage_prompt and os.getenv("GEMINI_API_KEY"):
            suggested = suggest_trans_cleavage_mutations(seq, closest_ref, args.trans_cleavage_prompt)
            mutants.extend(suggested)
        if len(mutants) < args.num_mutants_per_seq:
            mutants.extend(propose_mutations(seq, num_mutants=args.num_mutants_per_seq - len(mutants), num_sites=args.num_mutations))

        added_any = False
        for mut_seq, desc in mutants[: args.num_mutants_per_seq]:
            score = engine.score_candidate(mut_seq)
            ident = max_identity_to_refs(mut_seq, ref_seqs)
            if score >= args.min_esm_score and ident < args.max_identity:
                mut_id = f"{rec.id}_mut_{desc}_id{ident:.2f}"
                kept.append((mut_id, mut_seq, score, ident))
                added_any = True
        if not added_any and base_identity < args.max_identity:
            kept.append((rec.id, seq, engine.score_candidate(seq), base_identity))

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        for mid, mseq, score, ident in kept:
            f.write(f">{mid}_esm{score:.3f}\n{mseq}\n")
    print(f"[+] Wrote {len(kept)} drift variants to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
