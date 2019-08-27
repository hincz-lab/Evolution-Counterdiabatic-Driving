import sys
import numpy as np

data = []

with open(sys.argv[1]) as infile:
    data = [float(i.strip()) for i in infile.readline().split(",")]

min_val = min([i for i in data if i  > 0])
print(min_val)

# data = [(max_val - i)/max_val for i in data]

data = [np.format_float_positional(i/min_val) for i in data]

#data = [i + max_val for i in data]

with open(sys.argv[1], "w") as outfile:
    outfile.write(",".join([str(i) for i in data]))