import sys
import numpy as np

# This script takes in a file with fitness values separated by commas
# and converts the values to be s values (relative fitness as used in 
# the model) instead. 
# WARNING: Overwrites given file!

if len(sys.argv) < 2:
    print("Usage: python convert_fitness_to_s.py [name of file to convert]")

data = []

# Read in fitness values
with open(sys.argv[1]) as infile:
    data = [float(i.strip()) for i in infile.readline().split(",")]

# Do conversion
data = [np.format_float_positional(data[-1]/i - 1) if i != 0 else 10000000000000 for i in data]

# Write out s values
with open(sys.argv[1], "w") as outfile:
    outfile.write(",".join([str(i) for i in data]))