""" !!! Use median_string.py since this algorithm is very slow and inefficient. Just coded for an alternative exploration."""

import math
from itertools import product

def readFile(filename):
    """Reading a file and returning a list of lines"""
    with open(filename, 'r') as f:
        return [l.strip() for l in f]

sequences = readFile("regulatory_motif_week_3/median_string.txt")

def hamming_distance(seq1, seq2):
    return len([(n1, n2) for n1, n2 in zip(seq1,seq2) if n1 != n2])

def suffix(pattern):
    return pattern[1:] # simply remove first character

def neighbors(pattern,d):
    """ Compared to generating kmers with itertools product 4^k kmers at once, neighbors() is a lot more inefficient
    """
    if d == 0:
        return {pattern}  # wrap in set  # no mismatches allowed, only pattern itself
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}  # any single nucleotide is valid
    
    neighborhood=set() #sets use curly braces {} but watch out — empty {} creates a dict, not a set

    suffix_neighbors=neighbors(suffix(pattern), d)
    for text in suffix_neighbors:
        if hamming_distance(suffix(pattern), text) < d:
            for nuc in ["A","T","G","C"]:
                neighborhood.add(nuc+text)
        else:
            neighborhood.add(pattern[0]+text)
    return neighborhood

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
    """Tries every possible kmer pattern of length k, finds the one with the smallest total distance to all sequences — that's the median string, the motif that best represents all sequences simultaneously"""
    distance = math.inf          
    for seq in sequences:
         for i in range(len(seq) - k + 1):
            window = seq[i:i+k]
            for neighbor in neighbors(window,k): #generates all posssible d=k kmers for every window -> very slow compared to median_string.py, where 4^k kmers are generated at once
                if distance > DistanceBetweenPatternAndStrings(neighbor, sequences):
                    distance = DistanceBetweenPatternAndStrings(neighbor, sequences)
                    median = neighbor      # save pattern with smallest total distance
    return median

print(medianstring(sequences,6)) #k-mer pattern that minimizes total hamming distance across all sequences -> median string

