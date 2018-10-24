import numpy as np
from CIS import *

basis = ((0,0), (0,1), (1,0), (1,1), (2,0), (2,1))
Beryllium = CIS(4,basis)
A = Beryllium.organize()

eigvals, eigvecs = np.linalg.eigh(A)

print(A)
print(eigvals)
