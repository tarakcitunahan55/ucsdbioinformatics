from collections import Counter

def readFile(filename):
    with open(filename, 'r') as f:
        for l in f:
            return l.strip() #genome is only on line 1 of .txt file
genome = readFile("kmer_clump_data.txt")

k_len=9
clump_window_len=500
repeat_no_clump=3

def get_clump(genome, k, L, t):
    clump_kmers = set()
    
    # count kmers in first window
    kmer_counts = Counter(genome[i:i+k] for i in range(L - k + 1))
    
    # add clumps from first window. 
    for kmer, count in kmer_counts.items():
        if count >= t:
            clump_kmers.add(kmer)
    
    # slide window - only check incoming kmer. Slides the window one position at a time across the entire genome, starting from position 1 (since position 0 was handled above).
    for start in range(1, len(genome) - L + 1):
        outgoing = genome[start - 1:start - 1 + k] #The kmer that is leaving the window as it slides right by one position
        incoming = genome[start + L - k:start + L] #The kmer that is entering the window from the right
        
        kmer_counts[outgoing] -= 1 
        kmer_counts[incoming] += 1 #Updates the counter incrementally — instead of recomputing all kmer counts from scratch, just subtract the outgoing kmer and add the incoming one. This is what makes it O(n) fast.
        
        if kmer_counts[incoming] >= t:  # only check incoming-— it's the only count that increased. Previously found clumps are already safely stored in clump_kmers and never removed, so we don't need to check outgoing. When a kmer leaves the window its count decreases — it can only lose clump status, not gain it.
            #If you added a check for kmer_counts[outgoing] >= t, and outgoing is still >= t after decrementing, it was already in clump_kmers from a previous window anyway — so set.add() just silently ignores the duplicates. Not wrong, just unnecessary work.
            clump_kmers.add(incoming)
    
    return clump_kmers

result = get_clump(genome, k_len, clump_window_len, repeat_no_clump)
print(len(result))   # number of distinct clump kmers
print(*result)       # the kmers themselves

"""Alternative, BUT SLOWER:

from collections import Counter

DNAStr = "AAA"
k = 9
l = 500
t = 3

def find_t_frequent_kmers(seq, k, t):
    kmers=[]
    for i in range(len(seq) - k + 1):
        kmers.append(seq[i:i+k])
    kmer_frequencies = Counter(kmers)

    result = []
    for kmer, count in kmer_frequencies.items():
        if count >= t:
            result.append(kmer)
    return result

def clump_finder(k, l):
    clump_kmers = set()
    for i in range(len(DNAStr) - l + 1):
        possible_oric = DNAStr[i:i+l]
        clump_kmers.update(find_t_frequent_kmers(possible_oric, k, t)) #update unpacks and adds the list's elements to the set
    return " ".join(clump_kmers)

print(clump_finder(k, l))"""

""" NOT CORRECT - compares all pairs of positions against each other across the entire genome, so positions that are each within 500 of a neighbor 
but spread across a total range greater than 500 can all end up in store, falsely appearing to form a clump within a single window of length L

from collections import Counter

def readFile(filename):
    with open(filename, 'r') as f:
        for l in f:
            return l.strip() #genome is only on line 1 of .txt file

genome = readFile("kmer_clump_data.txt")
k_len=9
kmer_list = set(genome[i:i+k_len] for i in range(len(genome) - k_len + 1)) #set() ensures same kmer does not get stored twice or more

holder={}

for kmer in kmer_list:
    holder[kmer] = []  #initialize an empty list for every unique kmer -> {kmer: []}
    for position in range(len(genome) - k_len + 1):
        if genome[position:position + k_len] == kmer:
            holder[kmer].append(position) #{"A": [1, 2, 3, 4]}

clump_window_len=500
repeat_no_clump=3

def get_clump():
    clump_kmers = []
    for key, positions in holder.items():
        store=set() #initialize an empty set for every unique kmer
        for i in range(len(positions)):
            for j in range(i+1,len(positions)):
                if i != j:
                    if abs(positions[j]-positions[i])<clump_window_len: #absolute value difference
                        store.add(positions[j])
                        store.add(positions[i])
        if len(store)>=repeat_no_clump: #find kmers with repeat number at least repeat_no_clump
            clump_kmers.append(key)
    return clump_kmers

print(*get_clump())    """
                   

