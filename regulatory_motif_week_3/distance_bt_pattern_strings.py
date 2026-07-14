import math
from itertools import product


def read_file(file): # data was given on one line including the pattern we search for and DNA sequences
    with open(file, "r") as f:
        first_line = f.readline().rstrip()
        args = first_line.split()
        pattern = args[0]
        dna_list = []
        for i in range(1, len(args)):
            dna_list.append(args[i])

    return pattern, dna_list

pattern, sequences = read_file("regulatory_motif_week_3/distance_bt_.txt")


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
    k=len(pattern)
    distance=0
    for seq in sequences:
        hammingdistance=math.inf
        for kmer in get_window(seq,k):
            if hammingdistance > hamming_distance(pattern, kmer):
                hammingdistance=hamming_distance(pattern, kmer)
        distance+=hammingdistance
    return distance

print(DistanceBetweenPatternAndStrings(pattern,sequences))

