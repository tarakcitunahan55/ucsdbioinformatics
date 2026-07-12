import itertools
from collections import Counter
sequence="TGATGCCGCTTGCCTCCCGCTTGCTGCCTCCCTCCTGCTGATGCTGCCGCTCGCTCCGTCGCTTGATGCTGCCCGTCCGTTGACCGTTGCTGACCGTCTCCCCGTCGCTCTCCCGCTCGCTTGATGCTGCCCGTTGCTGCCTCCCCGTCTCCCGCTCTCCCGCTTGCCGCTTGACCGTTGACCGTCGCTCCGTTGCCTCCTGCTGCTGACCGTCGCTCCGTCTCCTGCTGCCGCTTGATGCCGCTTGACTCCCCGTTGATGCCTCCTGCTGCTGCCCGTCTCCCCGTCCGTCGCTTGACCGTTGCTGACTCCTGCCGCTTGACCGTCTCCCGCTCTCCTGACTCCCGCTCCGTTGACTCCCTCCTGA"
#search for DNA A boxes in a known ori sequence

d=3 #number of mismatches at most
k=7 #kmer length
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

