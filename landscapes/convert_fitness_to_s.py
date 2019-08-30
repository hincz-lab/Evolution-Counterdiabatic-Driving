import sys
import numpy as np

data = []

with open(sys.argv[1]) as infile:
    data = [float(i.strip()) for i in infile.readline().split(",")]

max_val = max(data)

# data = [(max_val - i)/max_val for i in data]

data = [np.format_float_positional(data[-1]/i - 1) if i != 0 else 10000000000000 for i in data]

#data = [i + max_val for i in data]

with open(sys.argv[1], "w") as outfile:
    outfile.write(",".join([str(i) for i in data]))