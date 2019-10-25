import sys
import numpy as np

# This script makes a transition probability matrix for the specified mutation rate.
# Assumes that only mutations between adjacent genotypes are possible and all have equal
# probability.

# Adapted from https://github.com/mirzaevinom/fitness_landsace_abm/blob/master/abm_model.py

# Computes the hamming distance between two genotypes.
def hammingDistance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

# This script requires two command-line arguments: the number of genotypes to use and
# the total probability of mutating away from the current genotype.
if len(sys.argv) != 3:
    print("Usage: python make_trans_prob.py [number of genotypes] [mutation rate]")

n_genotypes = int(sys.argv[1])
mut_rate = float(sys.argv[2])

# Initiallize empty transition matrix
probs = np.zeros([n_genotypes, n_genotypes])

# Calculate hamming distance between all pairs of genotypes
for i in range(n_genotypes):
    for j in range(n_genotypes):
        probs[i,j] = hammingDistance(format(i, '010b'), format(j, '010b'))

# It is impossible to mutate to a genotype more than 1 bit-flip away
probs[probs>1] = 0
# Divide each 1 by the number of neighbors each genotype has, to get
# the fraction of the time a mutation will go to that genotype specifically
probs = probs/probs.sum(axis=1)
# Multiply in the actual probability of a mutation happening
probs *= mut_rate

# Handle the diagonal (probability of not mutating)
for i in range(n_genotypes):
    probs[i,i] = 1 - mut_rate

# Print matrix
for row in probs:
    print(",".join([str(i) for i in row]))