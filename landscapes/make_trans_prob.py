import sys
import numpy as np

# Adapted from https://github.com/mirzaevinom/fitness_landsace_abm/blob/master/abm_model.py

# Computes the hamming distance between two genotypes.
def hammingDistance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

if len(sys.argv) != 3:
    print("Usage: python make_trans_prob.py [number of genotypes] [mutation rate]")

n_genotypes = int(sys.argv[1])
mut_rate = float(sys.argv[2])

probs = np.zeros([n_genotypes, n_genotypes])

for i in range(n_genotypes):
    for j in range(n_genotypes):
        probs[i,j] = hammingDistance(format(i, '010b'), format(j, '010b'))

probs[probs>1] = 0
probs = probs/probs.sum(axis=1)
probs *= mut_rate

for i in range(n_genotypes):
    probs[i,i] = 1 - mut_rate

for row in probs:
    print(",".join([str(i) for i in row]))