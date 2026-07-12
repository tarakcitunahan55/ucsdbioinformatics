from collections import Counter
from itertools import product

genome = "ATGTGCGGCGCTGCGGCCTGTGCTGCATGTGCCTGGCATGCTGGCATGCTGATGTGCGCATGCTGTGCATGCTGATGTGCATGGCTGCTGCGGCTGCCTGATGATGGCCTGCTGTGCGGCGCGCTGCATGGCTGCGGCATGTGCGCTGCGCTGCATGTGCGCGCGCGCGGCATGGGCATGGGCTGCTGCGGCATGCTGATGGGC"
#search for DNA A boxes in a known ori

k = 5
d = 3

def all_kmers(k):
    kmers = []
    for p in product('ACGT', repeat=k): #generate all possible 4^k number of kmers
        kmers.append(''.join(p))
    return kmers

def hamming_distance(seq1, seq2):
    count = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1
    return count #one line version -> return len([(n1, n2) for n1, n2 in zip(seq1,seq2) if n1 != n2])


def reverse_complement(DNAStr):
    table = str.maketrans("ATGC", "TACG")
    newStr = DNAStr.translate(table)
    return newStr[::-1]

def most_frequent(kmer_list):
    kmer_frequencies = Counter(kmer_list)
    highest_frequency = max(kmer_frequencies.values())
    max_list = []
    for key, value in kmer_frequencies.items():
        if value == highest_frequency:
            max_list.append(key)
    return max_list

def approximate_patterns():
    results = []
    for kmer in all_kmers(k):
        count = 0
        for i in range(len(genome) - k + 1):
            if hamming_distance(kmer, genome[i:i+k]) <= d:
                count += 1
            if hamming_distance(reverse_complement(kmer), genome[i:i+k]) <= d:
                count += 1 #add kmer and RC counts together to find most freq kmer+its RC
        if count > 0:
            for _ in range(count):
                results.append(kmer)
    return most_frequent(results)

print(" ".join(approximate_patterns())) #same output as unpacking the list

""" NOT CORRECT: finds the most frequent kmer(s) first, then separately finds the most frequent RCs from most freq kmer(s). It never adds kmer + RC counts together

import itertools
from collections import Counter
sequence="ATGTGCGGCGCTGCGGCCTGTGCTGCATGTGCCTGGCATGCTGGCATGCTGATGTGCGCATGCTGTGCATGCTGATGTGCATGGCTGCTGCGGCTGCCTGATGATGGCCTGCTGTGCGGCGCGCTGCATGGCTGCGGCATGTGCGCTGCGCTGCATGTGCGCGCGCGCGGCATGGGCATGGGCTGCTGCGGCATGCTGATGGGC"
d=3 #number of mismatches at most
k=5 #kmer length
kmers = [''.join(kmer) for kmer in itertools.product('ATGC', repeat=k)] #generate all possible 4^k number of kmers


def get_window():
     return [(sequence[pos:pos+k]) for pos in range(len(sequence)-k+1)] #0-based indexing

windows=get_window()

def most_freq_kmer_with_mismatch():
    store=[]
    for kmer in kmers:
        for window in windows:
            if len([(n1, n2) for n1, n2 in zip(window,kmer) if n1 != n2])<=d: #alternative: if sum(n1 != n2 for n1, n2 in zip(window[0], pattern)) <= d:
                store.append(kmer) #appends kmers

    kmer_frequencies = Counter(store)
    highest_frequency = max(kmer_frequencies.values())

    return [kmer for kmer, frequency in kmer_frequencies.items() if frequency == highest_frequency]

print (*most_freq_kmer_with_mismatch())

def reverse_complement():
    holder=[]
    for kmer in most_freq_kmer_with_mismatch():
        mapping = str.maketrans('ATCG', 'TAGC')
        holder.append(kmer.translate(mapping)[::-1])
    return holder

def most_freq_rc_with_mismatch(): #for reverse complements
    store=[]

    for kmer in reverse_complement():
        for window in windows:
            if len([(n1, n2) for n1, n2 in zip(window,kmer) if n1 != n2])<=d: #alternative: if sum(n1 != n2 for n1, n2 in zip(window[0], pattern)) <= d:
                store.append(kmer) #appends positions

    kmer_frequencies = Counter(store)
    highest_frequency = max(kmer_frequencies.values())

    return [kmer for kmer, frequency in kmer_frequencies.items() if frequency == highest_frequency]

print (*most_freq_rc_with_mismatch())"""

