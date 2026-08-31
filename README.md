# UCSD Bioinformatics I — Finding Hidden Messages in DNA

Python implementations of core algorithms from the Finding Hidden Messages in DNA course (genome replication origins and regulatory motif discovery).

## Contents

### `dnaabox_ori_kmers_week1_2/` — Replication Origin Detection
Algorithms for locating the origin of replication (*oriC*) in bacterial genomes.

| File | Description |
|---|---|
| `skew.py` | Computes GC-skew to narrow down the likely origin region |
| `find_ori_from_skew.py` | Identifies *ori* candidates from skew minimum |
| `pattern_with_hamming_distance.py` | Approximate pattern matching using Hamming distance |
| `neighborhood_mismatch.py` / `neighborhood_freq_kmer_rc_mismatch.py` | Generates k-mer neighborhoods within a mismatch threshold |
| `freq_kmer_with_mismatches.py` / `freq_kmer_and_rc_with_mismatch.py` | Frequent k-mers allowing mismatches and reverse complements |
| `kmer_clump.py` | Detects genomic regions with unusually dense k-mer clumps |
| `find_dnaabox_final.py` | Combines the above to predict the DnaA box / *oriC* |final product of all codes
| `cholera_motif.py` |

### `regulatory_motif_week_3_4_5/` — Motif Finding
Algorithms for discovering conserved regulatory motifs across a set of DNA sequences.

| File | Description |
|---|---|
| `median_string.py` / `median_string_by_neighbors.py` | Median String algorithm for motif finding |
| `profile_most_probable_kmer.py` | Finds the most probable k-mer given a profile matrix |
| `greedy_motif_search.py` / `greedymotifsearch_laplace_pseudocount.py` | Greedy Motif Search, with Laplace smoothing |
| `randomized_motif_search.py` | Randomized Motif Search |
| `gibbs_sampler.py` | Gibbs Sampling for motif discovery |
| `motif_entropy_calc.py` | Computes motif matrix entropy |
| `distance_bt_pattern_strings.py` | Pairwise distance between pattern and string set |
| `common_motifs_in_all_strings_mismatch...py` | Finds motifs common to all sequences under mismatch constraints |

## Topics covered
Genome replication origin detection · GC-skew analysis · approximate pattern matching · greedy/randomized/Gibbs motif search · consensus and entropy scoring

## Notes
Coursework implementations; shared for portfolio purposes only, DO NOT COPY.
