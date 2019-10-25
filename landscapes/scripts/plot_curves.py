import numpy as np
import matplotlib.pyplot as plt

# This script reproduces the drug response curves in Ogbunugafor et. al, 2016
# to confirm that wer're doing everything right.

# Note that the equation expects log10 of drug concentration and IC 50
# Also note that we're working with molar concentrations while Ogbunugafor et. al
# used micromolar

# This is the value of c used by Ogbunugafor et. al, 2016 (obtained via personal communication)
c = -.6824968
# IC50 values are from supplemental table 1 of Ogbunugafor et. al, 2016/supplemental 
# data table 3 of Brown et al, 2010 (https://academic.oup.com/mbe/article/27/12/2682/1072079)
# (the same data are in both tables). They are log10.
logic50s = [-6.286,-5.812,-4.239,1,-6.046,-5.774,-3.732,-3.55,-5.724,-5.491,-4.015,-4.6,-5.773,-5.624,-3.587,-3.3]
# These are the untransformed drugless growth values from supplemental 
# data table 3 of Brown et al, 2010. Normalized versions of these values were used
# by Ogbunugafor et. al, 2016 
g_druglesses = [0.000969794, 0.000884475, 0.000851618, 0, 0.000950368, 0.000953728, 0.000969172, 0.000845918, 0.000776222, 0.000821543, 0.000906315, 0.000693783, 0.000883164, 0.000889632, 0.001005913, 0.000867504]


xs = []
ys = []
for _ in g_druglesses:
    ys.append([])

for i in np.logspace(-9, 5, num=100000):
    xs.append(i)
    for genotype in range(len(logic50s)):
        ys[genotype].append((g_druglesses[genotype]/(1+np.exp((logic50s[genotype] - np.log10(i))/c)))/0.000693783)
    

for genotype in range(len(logic50s)):
    plt.semilogx(xs, ys[genotype])

plt.show()