import math
from itertools import product

def readFile(filename):
    """Reading a file and returning a list of lines"""
    with open(filename, 'r') as f:
        return [l.strip() for l in f]

sequences = readFile("regulatory_motif_week_3_4_5/median_string.txt")

def hamming_distance(seq1, seq2):
    return len([(n1, n2) for n1, n2 in zip(seq1,seq2) if n1 != n2])

def all_kmers(k):
    kmers = []
    for p in product('ACGT', repeat=k): #generate all possible 4^k number of kmers
        kmers.append(''.join(p))
    return kmers

def get_window(sequence,k):
    return [(sequence[pos:pos+k]) for pos in range(len(sequence)-k+1)]

def DistanceBetweenPatternAndStrings(pattern, sequences):
    """For each sequence, finds the closest kmer to the pattern (minimum hamming distance), then sums these minimums across all sequences:
"""
    k=len(pattern)
    distance=0
    for seq in sequences:
        hammingdistance=math.inf # start with infinity
        for kmer in get_window(seq,k):
            if hammingdistance > hamming_distance(pattern, kmer):
                hammingdistance=hamming_distance(pattern, kmer) # keep minimum hamming distance
        distance+=hammingdistance # add minimum to total
    return distance #total hamming distance of all sequences

def medianstring(sequences, k):
    """Tries every possible kmer pattern of length k, finds the one with the smallest total distance to all sequences — that's the median string, the motif that best represents all sequences simultaneously
    Although very accurate, very slow and impractical for large k"""
    distance = math.inf           # start with infinity, any real distance will be smaller
    for pattern in all_kmers(k):  # try every possible pattern/possible consensus strings
        if distance > DistanceBetweenPatternAndStrings(pattern, sequences):
            distance = DistanceBetweenPatternAndStrings(pattern, sequences)
            median = pattern      # save pattern with smallest total distance
    return median

print(medianstring(sequences,6)) #k-mer pattern that minimizes total hamming distance across all sequences -> median string

