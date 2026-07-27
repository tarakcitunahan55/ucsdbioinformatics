from collections import Counter

def read_file(file):
    with open(file, "r") as f: #k,t,sequences all in one line
        args = f.readline().rstrip().split()
    k = int(args[0])
    t = int(args[1])
    dna_list = args[2:]
    return k, t, dna_list

def build_profile(motifs, k):
    t = len(motifs)
    profile = {"A": [], "C": [], "G": [], "T": []}
    for pos in zip(*motifs):           # iterate pos directly
        counts = Counter(pos)
        for nuc in "ACGT":
            profile[nuc].append(counts.get(nuc, 0) / t) #.get(nuc, 0) means "give me the count of this nucleotide, or 0 if it doesn't exist in the Counter." Without it you'd get a KeyError every time a nucleotide doesn't appear
    return profile

def profile_most_probable(sequence, k, profile):
    best_prob = -1 # start lower than any possible probability
    best_kmer = sequence[:k]  # default to first kmer 
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k] #Slides window across sequence generating each kmer
        prob = 1
        for pos, nuc in enumerate(kmer): #For each kmer, multiplies profile probabilities at each position (no pseudocounts/laplace)
            prob *= profile[nuc][pos] #A single zero anywhere in the profile kills the entire probability to 0, which is why pseudocounts matter — without them, any unseen nucleotide at any position gives probability 0.
        if prob > best_prob: #Keeps track of kmer with highest probability
            best_prob = prob
            best_kmer = kmer
    return best_kmer

def score(motifs, k):
    t = len(motifs)
    total = 0
    for pos in zip(*motifs): # each pos is a tuple of nucleotides
        counts = Counter(pos)
        total += t - max(counts.values()) #at each position gives the number of nucs that don't match the most frequent nucleotide at that position, aka scoring
    return total

def greedy_motif_search(dna_list, k, t):
    """Faster than median string finding for larger k, though not very accurate especially if you do not use pseudocounts"""
    best_motifs = [seq[:k] for seq in dna_list] #Initializes best motifs as first k-mer from each sequence
    
    for i in range(len(dna_list[0]) - k + 1):
        motif1 = dna_list[0][i:i+k] #Tries every k-mer in the first sequence as the starting motif
        motifs = [motif1] #Starts a fresh motif list with just motif1
        
        for j in range(1, t): #For each subsequent sequence, builds profile from current motifs and finds most probable k-mer
            profile = build_profile(motifs, k)
            motifs.append(profile_most_probable(dna_list[j], k, profile)) #j=2, j=3, j=4: same process, profile grows with each motif added
        
        if score(motifs, k) < score(best_motifs, k): #After collecting all t motifs, checks if this collection scores better (lower) than current best
            best_motifs = motifs
    
    return best_motifs #After trying all k-mers in first sequence as starting motif, returns the collection with lowest score

k, t, dna_list = read_file("regulatory_motif_week_3/greedy_motif.txt")
print(" ".join(greedy_motif_search(dna_list, k, t)))