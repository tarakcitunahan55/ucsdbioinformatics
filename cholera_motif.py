""" Use the whole cholera genome to identify DNAA protein binding specific 9-mer positions 
and origin of replication that is implied by the 9-mer motifs at most 1000 bp apart"""

def readFile(filename):
    with open(filename, 'r') as f:
        for l in f:
            return l.strip() #genome is only on line 1 of .txt file

genome = readFile("cholera_genome.txt")


motif="CTTGATCAT" #a 9-mer that DNA A protein binds to in origin of replication. We want to see this 9-mer's starting position in whole genome.

def get_motif(sequence):
    return [position for position in range(len(sequence) - len(motif) + 1) #return position in 0-based indexing of python
            if sequence[position:position + len(motif)] == motif]

positions=get_motif(genome)

print("Position:", *positions)  #gives the initial position the motifs are found in sequence by unpacking the list
print("Count:", len(positions)) #the number of motifs in sequence

store=set() # use a set to avoid duplicates
for i in range(len(positions)):
        for j in range(len(positions)):
            if i != j:
                if abs(positions[j]-positions[i])<1000: #absolute value difference
                    store.add(positions[i])
                    store.add(positions[j])

print("9-mer positions at most 1000 bp apart:", *store)

#form clumps, i.e., appear close to each other in a small region of the genome. You may check that the same conclusion is reached when searching for ATGATCAAG (its reverse complement, another 9-mer that DNA A binds)
#evidence that ATGATCAAG/CTTGATCAT 9-mers represent the message to DnaA to start replication (DnaA box in origin of replication)                