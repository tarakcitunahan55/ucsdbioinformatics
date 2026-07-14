sequence="TATATGACGCAGACTAGAGTACGCGCAAGTACTACTGCGTAGGAAGAAGGTTGAGTGACGCCGTGATACGGAGCGAGTGATTGAGATCACATTCAAGGTGGGCAGGATGCACGCGAAAGTTAAGGAGCAGTACGTCTAATCGAGCCATGCCTCCGATGCCGAGGTCAAAATGCTCGCGCCTCATGACCTTGGGATACGATATCGCAAGACCATACCGGAACGACGACTAGATAAGTGGCGCATTTATCACGGGAAGAGCGCAGAGATAAGATGAGCGCCGAGACAAGAGGAATTCAGCCTTAAGCACATACCCTATCTCGAACATACTAACCTATTCAGAAACGGTACGCAGGACCATTTATTGGCCGTATAAGCGTTGTTTGCGATTTCGAGTGGCGCGTCAAAGAGCCTCGTACGGGAGGCAGTAGATTCCGGGCACTGTGAGTCCGGATCGGGCATAAATGGACGTACAAGATAGAGTGCGGCAATACCGTTTATCTGCAGTGCTGGTGCAAGTGCGGAAGCATAACCAGCCCTACGGCCAGCTTTTGTATCAGAGACCAGGATGACTCTAGAGTACCCACGAGGTCAAGCTACCTCAATTGCTAAGCAGTTTAGGCACCACTCCTGTGAGGCCTCTAACCAGTAACACTCATCAACGATTCGCAACGCAGAAGGAGAACACAATGGTTGTTGGGGTGAAACATTGTCTTGGCCCCTCCGCAATGAAGCTATATGTCGTACTACCGCGGTGGGGTTAAGAAATTTGCGCCCATGAGGAAGTCGCGCTAGGATCCTGCAGAGACATGAACAGAGCCCGCCATGCTAGGTCGGGTTGAGATGCATTTAGGAAATAAAGCTACTCGGAGAGACCAGGTACGAACTTTGTAAATGCTGAGGTTCTATGTATAAGCATCGGGAACAGGCAGTAGGGGACTGTCCTCACGGTACGGCTCTCTATGCGACCCAAAATATTGAGTCAATACGACAAGTGGCACGT"
k=14

def read_file(file):
    """Form lists for each nuc in every position from the given profile matrix"""
    with open(file, "r") as f:
        a_list = f.readline().split() #split() returns a list of strings
        c_list = f.readline().split() #readline automatically moves to the next line in the file
        g_list = f.readline().split()
        t_list = f.readline().split()

    return a_list, c_list, g_list, t_list     

a_list, c_list, g_list, t_list = read_file("regulatory_motif_week_3/profile_kmer.txt")

def get_window(sequence,k):
    return [(sequence[pos:pos+k]) for pos in range(len(sequence)-k+1)]

def profile_kmer():
    holder={}
    for seq in get_window(sequence,k):
        product=1
        for i in range(len(seq)):
            if seq[i]=="A":
                product=product*float(a_list[i])
            elif seq[i]=="C":
                product=product*float(c_list[i])
            elif seq[i]=="G":
                product=product*float(g_list[i])
            else: #"T"
                product=product*float(t_list[i])
        holder[seq]=product
    
    #maxseq = max(holder, key=holder.get) -> gets the key corresponding to max value of holder dictionary
    #return maxseq -> if there is a tie, it only returns one 

    highest_product = max(holder.values())
    return [seq for seq, product in holder.items() if product == highest_product] #if there is a tie, returns all corresponding keys

print (*profile_kmer())
