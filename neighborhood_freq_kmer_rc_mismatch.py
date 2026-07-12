genome = "ATGTGCGGCGCTGCGGCCTGTGCTGCATGTGCCTGGCATGCTGGCATGCTGATGTGCGCATGCTGTGCATGCTGATGTGCATGGCTGCTGCGGCTGCCTGATGATGGCCTGCTGTGCGGCGCGCTGCATGGCTGCGGCATGTGCGCTGCGCTGCATGTGCGCGCGCGCGGCATGGGCATGGGCTGCTGCGGCATGCTGATGGGC" 
#search for DNA A boxes in a known ori

k = 5
d = 3

def hamming_distance(seq1, seq2):
    return len([(n1, n2) for n1, n2 in zip(seq1,seq2) if n1 != n2])


def reverse_complement(DNAStr):
    table = str.maketrans("ATGC", "TACG")
    newStr = DNAStr.translate(table)
    return newStr[::-1]


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

def most_frequent_from_counts(counts):
    highest_frequency = max(counts.values())
    max_list = []
    for key, value in counts.items():
        if value == highest_frequency:
            max_list.append(key)
    return max_list

def approximate_patterns():
    counts = {}

    for i in range(len(genome) - k + 1):
        window = genome[i:i+k]
        window_rc = reverse_complement(window)

        for neighbor in neighbors(window, d):
            if neighbor in counts:
                counts[neighbor] += 1
            else:
                counts[neighbor] = 1

        for neighbor in neighbors(window_rc, d):
            if neighbor in counts:
                counts[neighbor] += 1
            else:
                counts[neighbor] = 1

    return most_frequent_from_counts(counts)
    
print(*approximate_patterns())