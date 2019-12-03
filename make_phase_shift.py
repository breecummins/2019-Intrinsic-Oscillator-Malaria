import sys
import numpy as np
import random

fname = sys.argv[1]

if fname[-3:] == "csv":
    delimiter = ","
elif fname[-3:] == "tsv":
    delimiter = "\t"
else:
    raise ValueError("File type not recognized. Use .csv or .tsv format.")

A = np.loadtxt(open(sys.argv[1]), delimiter=delimiter, dtype = str)

for j in range(1,A.shape[0],1):
    phi = random.choice(range(len(A[j])-1))
    B = A[j].copy()
    for i in range(len(A[1])-1):
        if 1+i+phi < len(B):
            A[j,i+1] = B[1+ i+phi]
        else:
            A[j,i+1] = B[1+ i+phi-(len(B) -1)]

new_filename = (sys.argv[1])[:-4] + '_shifted.csv'
np.savetxt(new_filename, A, delimiter = ',', fmt = '%s')
