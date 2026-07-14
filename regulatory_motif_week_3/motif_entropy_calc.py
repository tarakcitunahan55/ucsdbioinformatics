from collections import Counter
import math #for logs

def readFile(filename):
    """Reading a file and returning a list of lines"""
    with open(filename, 'r') as f:
        return [l.strip() for l in f]

sequences = readFile("regulatory_motif_week_3/motif_entropy.txt")
t=len(sequences) #number of strings

counts = [Counter(position) for position in zip(*sequences)] #make a Counter dict for every position/column of sequences

def entropy():
    """ Entropy of a motif matrix is defined as the sum of the entropies of its columns
    Find individual entropies of every nucleotide frequency in every position and add all entropies together in all columns
    H=-sigma(pxlog2p) where p=nuc frequency in column
    """
    holder=0
    for counter in counts: #for every column
        holder+=sum(counter[key]/t*math.log2(counter[key]/t) for key in counter if counter[key]>0) #divide by t to find nuc freq in column
        #put "if counter[key]>0", since log2(0) is undefined
    return -holder #change the negative sign to positive by multiplying by negative -> (- x - = +)

print(entropy())
