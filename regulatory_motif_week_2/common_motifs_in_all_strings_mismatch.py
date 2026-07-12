k=5
d=2
sequence1="CCCGTGCAGCTAGTTTCGTGTATAA"
sequence2="AGCGCGAAGGTCGTGAATCCGGTCG"
sequence3="GCTGGGATAATGGGTACATGTGATA"
sequence4="TTTCCTACAATCTTGCCGCGGACAG"
sequence5="TACAGAGGAAGCGTGATCCAACGGT"
sequence6="TAACGGCCTGAATCACAAATAATAA"

all_sequences=[sequence1,sequence2,sequence3,sequence4,sequence5,sequence6]

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

def approximate_patterns():
    holder=set()
    for sequence in all_sequences:

        for i in range(len(sequence) - k + 1):
            window = sequence[i:i+k]
         

            for neighbor in neighbors(window, d):
                # check if neighbor appears in ALL sequences with at most d mismatches
                if all(
                    any(hamming_distance(neighbor, seq[j:j+k]) <= d 
                        for j in range(len(seq) - k + 1))
                    for seq in all_sequences
                ):
                    holder.add(neighbor)
       
    return holder


print("All (k, d)-motifs in Dna", *approximate_patterns())