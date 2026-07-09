""" Core idea: Generate all kmers within hamming distance d of a pattern by solving a smaller version of the same problem 
— recursive function - just like factorial solves 5! by first solving 4!.
"""

def hamming_distance(seq1, seq2):
    return len([(n1, n2) for n1, n2 in zip(seq1,seq2) if n1 != n2])

def suffix(pattern):
    return pattern[1:] # simply remove first character

def neighbors(pattern,d):
    """ Recursive function """
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
        

print(*neighbors("CAA",1)) #unpack the set