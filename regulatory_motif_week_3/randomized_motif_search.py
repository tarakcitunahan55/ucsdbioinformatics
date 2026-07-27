from collections import Counter
import random

def read_file(file):
    with open(file, "r") as f:
        args = f.readline().rstrip().split()
    k = int(args[0])
    t = int(args[1])
    dna_list = args[2:]
    return k, t, dna_list

def build_profile_laplace(motifs, k):
    """Build profile matrix with Laplace pseudocounts to avoid zero probabilities"""
    t = len(motifs)
    profile = {"A": [], "C": [], "G": [], "T": []}
    for column in zip(*motifs):
        counts = Counter(column)
        for nuc in "ACGT":
            # add 1 to count (pseudocount) and 4 to denominator (one for each nucleotide)
            profile[nuc].append((counts.get(nuc, 0) + 1) / (t + 4))
    return profile

def profile_most_probable(sequence, k, profile):
    """Find the most probable kmer in sequence given profile matrix"""
    best_prob = -1
    best_kmer = sequence[:k]
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        prob = 1
        for pos, nuc in enumerate(kmer):
            prob *= profile[nuc][pos]
        if prob > best_prob:
            best_prob = prob
            best_kmer = kmer
    return best_kmer

def score(motifs):
    """Score = total mismatches from consensus at each position"""
    t = len(motifs)
    total = 0
    for column in zip(*motifs):
        counts = Counter(column)
        total += t - max(counts.values())
    return total

def randomized_motif_search(dna_list, k, t):
    """Single run of randomized motif search"""
    # randomly select one k-mer from each sequence
    motifs = [seq[random.randint(0, len(seq) - k):random.randint(0, len(seq) - k) + k]
              for seq in dna_list]
    
    # fix: ensure random selection is valid
    motifs = [seq[i:i+k] for seq in dna_list 
              for i in [random.randint(0, len(seq) - k)]]
    
    best_motifs = motifs[:]  # copy current motifs as best
    
    while True:
        profile = build_profile_laplace(motifs, k)  # build profile from current motifs
        motifs = [profile_most_probable(seq, k, profile) for seq in dna_list]  # update motifs
        
        if score(motifs) < score(best_motifs):
            best_motifs = motifs[:]  # update best if improved
        else:
            return best_motifs  # return if no improvement

def run_randomized_motif_search(dna_list, k, t, runs=1000):
    """Find longer motifs (important in real data) than greedy search but slower
    Run many (1000) times and return best result
    Although the motifs returned by RandomizedMotifSearch are slightly less conserved than the motifs returned by MedianString,
    RandomizedMotifSearch has the advantage of being able to find longer motifs (since MedianString becomes too slow for longer motifs)."""
    
    best_motifs = None
    best_score = float('inf')
    
    for _ in range(runs):
        motifs = randomized_motif_search(dna_list, k, t)
        current_score = score(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs
    
    return best_motifs

k, t, dna_list = read_file("regulatory_motif_week_3/randomized_motif_search.txt")
result = run_randomized_motif_search(dna_list, k, t, runs=1000)
print(" ".join(result))