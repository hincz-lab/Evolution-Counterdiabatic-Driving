import numpy as np
import matplotlib.pyplot as plt

g_drugless = 1.398
c = -.6824968
ic50 = 0.000005176068319505681
logic50 = -6.286


logic50s = [-6.286,-5.812,-4.239,1,-6.046,-5.774,-3.732,-3.55,-5.724,-5.491,-4.015,-4.6,-5.773,-5.624,-3.587,-3.3]
g_druglesses = [0.000969794, 0.000884475, 0.000851618, 0, 0.000950368, 0.000953728, 0.000969172, 0.000845918, 0.000776222, 0.000821543, 0.000906315, 0.000693783, 0.000883164, 0.000889632, 0.001005913, 0.000867504]

for genotype in range(len(logic50s)):
    print(g_druglesses[genotype]/(1+np.exp((10**logic50s[genotype] - 0)/c)), g_druglesses[genotype]/(1+np.exp((10**logic50s[genotype] - 0)/c))*2)


xs = []
ys = []
for _ in g_druglesses:
    ys.append([])

for i in np.logspace(-9, 5, num=100000):
    print(i)
    xs.append(i)
    for genotype in range(len(logic50s)):
        ys[genotype].append((g_druglesses[genotype]/(1+np.exp((logic50s[genotype] - np.log10(i))/c)))/0.000693783)
    # ys2.append(g_drugless/(1+np.exp((np.exp(logic50) - i)/c)))
    

for genotype in range(len(logic50s)):
    plt.semilogx(xs, ys[genotype])
# plt.semilogx(ys2)

plt.show()