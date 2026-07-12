pattern="TGTCTGT" #kmer
sequence="CAACCTGGGGACGTTCGACGTTAGGTAAAGTTGACTTGTACGCGACAACGTCTGTTCAGTCCCGGGAGCTGCCCCATCCGCATCCCCGCGTACCCTAGTCAACTATTCGTATCTGCATCATCTGTCACATGACAAACGGATTATAAGGCTAAGGTATCGCCCGACTTCGTACAGCCGCTGTAAAGCTGTCTAGTTCCTGTCTTGCTCTTAACTCAAGCTGTTCATGACTTCGCGGGGTGTTATGTGTCTGTACAATCATGCTCAGTGACGGAGTAACCCACAATTGGGTACACCATGTATAGTGTGTGGATCGACTATTCTCTAACAAAAGAGATATAAGCCTTCCCGAATCCCCCCAGAATTGATTCGCAACGTTCCTGCGCGGTGAC"
d=3

def get_window():
     """ Find all starting positions where pattern appears as a substring of sequence with at most d mismatches. """
     return [(sequence[pos:pos+len(pattern)],pos) for pos in range(len(sequence)-len(pattern)+1)] #0-based indexing

windows=get_window()

store=[]
for window in windows:
     if len([(n1, n2) for n1, n2 in zip(window[0],pattern) if n1 != n2])<=d: #alternative: if sum(n1 != n2 for n1, n2 in zip(window[0], pattern)) <= d:
          store.append(window[1]) #appends positions

print(*store) #prints position of pattern with at most d-mismatches
print(len(store)) #prints the total number of patterns with at most d-mismatches 
          
