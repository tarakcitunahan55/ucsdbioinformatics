from collections import Counter
import random

def read_file(file):
    with open(file, "r") as f:
        args = f.readline().rstrip().split()
    k = int(args[0])
    t = int(args[1])
    N = int(args[2]) #number of iterations per run
    dna_list = args[3:]
    return k, t, N, dna_list

def build_profile_laplace(motifs, k):
    """Build profile matrix with Laplace pseudocounts to avoid zero probabilities"""
    t = len(motifs)
    profile = {"A": [], "C": [], "G": [], "T": []}
    for column in zip(*motifs):
        counts = Counter(column)
        for nuc in "ACGT":
            profile[nuc].append((counts.get(nuc, 0) + 1) / (t + 4))
    return profile

def profile_randomly_generated_kmer(sequence, k, profile):
    """Select a k-mer randomly based on profile probabilities — not the most probable"""
    probabilities = []
    kmers = []
    
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        prob = 1
        for pos, nuc in enumerate(kmer):
            prob *= profile[nuc][pos]
        kmers.append(kmer)
        probabilities.append(prob)
    
    # normalize probabilities so they sum to 1
    total = sum(probabilities)
    probabilities = [p / total for p in probabilities]
    
    # randomly select a kmer weighted by its probability, high probability kmers chosen more often but not always
    return random.choices(kmers, weights=probabilities, k=1)[0]

def score(motifs):
    """Score = total mismatches from consensus at each position"""
    t = len(motifs)
    total = 0
    for column in zip(*motifs):
        counts = Counter(column)
        total += t - max(counts.values())
    return total

def gibbs_sampler(dna_list, k, t, N):
    """
    Random one motif update per iteration → occasional bad moves are intentional → 
            stopping at plateau would miss future improvements after bad moves
            so must run fixed N to allow recovery from bad random choices"""
    
    # randomly select one k-mer from each sequence
    motifs = [seq[random.randint(0, len(seq) - k):random.randint(0, len(seq) - k) + k]
              for seq in dna_list]
    motifs = [seq[i:i+k] for seq in dna_list
              for i in [random.randint(0, len(seq) - k)]]
    
    best_motifs = motifs[:]
    
    for j in range(N): #always runs exactly N times, no early stopping regardless of whether it's improving or not
        i = random.randint(0, t - 1)  # randomly select one sequence index
        
        # build profile from all motifs EXCEPT the i-th one
        motifs_except_i = motifs[:i] + motifs[i+1:]
        profile = build_profile_laplace(motifs_except_i, k)
        
        # replace i-th motif with a profile-randomly generated kmer
        motifs[i] = profile_randomly_generated_kmer(dna_list[i], k, profile)
        
        # update best if improved
        if score(motifs) < score(best_motifs):
            best_motifs = motifs[:]
    
    return best_motifs

def run_gibbs_sampler(dna_list, k, t, N, runs=20):
    """Run Gibbs sampler 20 times and return best result

     Profile-random selection occasionally picks suboptimal but allows escaping local optima to find global optima
     One motif update per iteration → preserves good motifs found so far -> better than randomized search
"""
    best_motifs = None
    best_score = float('inf')
    
    for _ in range(runs):
        motifs = gibbs_sampler(dna_list, k, t, N)
        current_score = score(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs
    
    return best_motifs

k, t, N, dna_list = read_file("regulatory_motif_week_3_4_5/gibbs_sampler.txt")
result = run_gibbs_sampler(dna_list, k, t, N, runs=20)
print(" ".join(result))