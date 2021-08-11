import numpy as np
# to convert the grain coordinates to a pdb file

filename = "fort.1401"
with open(filename, 'r') as f:
  data = f.readlines()
  data = [[float(j)/10 for j in i.split()] for i in data]

# for file formatting
header = """ TITLE     SIMULATION CELL
CRYST1  226.300  226.300  226.300  90.00  90.00  90.00 P 1           1
"""
line = "ATOM   {:>4d}    I GRN X {:>3d}     {:>7.3f} {:>7.3f} {:>7.3f}  1.00  1.00\n"

# get info we need
pdbname = filename + ".pdb"
test = np.array(data)
with open(pdbname, 'w') as f:
  f.write(header)
  for i, grain in enumerate(data):
    f.write(line.format(i, i, *grain))

  

