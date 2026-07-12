import requests
k = 9
d = 1 #mismatches

url = f"https://bioinformaticsalgorithms.com/data/Salmonella_enterica.txt"
response = requests.get(url)
FASTADict={}
FASTALabel = ""
for line in response.text.splitlines():
    line = line.strip()
    if '>' in line:
        FASTALabel = "Salmonella enterica subsp. enterica serovar Typhi str"
        FASTADict[FASTALabel] = ""
    else:
        FASTADict[FASTALabel] += line    

sequence=FASTADict[FASTALabel]

def ori_finder(seq):
    """ Find ori (with skew data) that is crucial to enable us search for most freq 9-mer and its RC - DNA A boxes"""
    skew=[0] #initial skew set to zero before we start looking at sequence
    for i in range(len(seq)): 
        if seq[i]=="C":
            skew.append(skew[i]-1)
        elif seq[i]=="G":
            skew.append(skew[i]+1)
        else: 
            skew.append(skew[i])

    min_skew=min(skew) #minimum skew (G-C difference)

    return [pos for pos in range(len(seq)) if skew[pos+1]==min_skew] #return position of minimum skew in 0-based indexing -> reflects transition from reverse half strand to forward half strand (G-C difference decreasing to increasing)(skew slope turning from negative to positive) -> ori

ori=ori_finder(sequence)[0] #there are two min_skew points 2 nuc away from each other in Salmonella genome. Take the first pos.

sequence_possible_1=sequence[ori:ori+501] # try a 500-nuc window either starting, ending, or centered at the position of minimum skew to locate ori's possible locations
sequence_possible_2=sequence[ori-500:ori+1]
sequence_possible_3=sequence[ori-250:ori+251]

print(f"Possible ori site starts at either {ori-500} or {ori-250} or {ori}")

all_sequences=[sequence_possible_1,sequence_possible_2,sequence_possible_3]


def hamming_distance(seq1, seq2):
    return len([(n1, n2) for n1, n2 in zip(seq1,seq2) if n1 != n2])


def reverse_complement(DNAStr):
    table = str.maketrans("ATGC", "TACG")
    newStr = DNAStr.translate(table)
    return newStr[::-1]


def suffix(pattern):
    return pattern[1:] # simply remove first character


def neighbors(pattern,d):
    """ Recursive function - gives kmers of at most d mismatches """
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
    """ 
Find most frequent kmer(s) and its RCs within each candidate ori window independently
Collect all winners into holder
If the same kmer wins in multiple possible ori sites, it appears once in the set (since it's a set)
"""
    holder=set()
    for possible_oric_site in all_sequences:

        counts = {}

        for i in range(len(possible_oric_site) - k + 1):
            window = possible_oric_site[i:i+k]
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

        for i in most_frequent_from_counts(counts):
            holder.add(i)
       
    return holder


print("Possible DNA A box 9-mers:", *approximate_patterns())