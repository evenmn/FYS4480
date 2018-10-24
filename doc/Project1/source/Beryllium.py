import numpy as np
from CIS import *

basis = ((0,0), (0,1), (1,0), (1,1), (2,0), (2,1))
Beryllium = CIS(4,basis)
A = np.zeros((5,5))

# --- <c|H|c> ---
A[0,0] = Beryllium.c_H_c()

# --- <c|H|p_i^a> ---
A[0,1] = Beryllium.c_H_ia(0,2)
A[0,2] = Beryllium.c_H_ia(1,3)
A[0,3] = Beryllium.c_H_ia(0,4)
A[0,4] = Beryllium.c_H_ia(1,5)

# --- <p_i^a|H|c> ---
A[1,0] = Beryllium.c_H_ia(0,2)
A[2,0] = Beryllium.c_H_ia(1,3)
A[3,0] = Beryllium.c_H_ia(0,4)
A[4,0] = Beryllium.c_H_ia(1,5)

# --- <p_i^a|H|p_j^b> ---
# <12|H|21>
A[1,1] = Beryllium.ia_H_jb(0,2,0,2)
A[1,2] = Beryllium.ia_H_jb(0,2,1,3)
A[2,1] = Beryllium.ia_H_jb(1,3,0,2)
A[2,2] = Beryllium.ia_H_jb(1,3,1,3)

# <12|H|31>
A[1,3] = Beryllium.ia_H_jb(0,2,0,4)
A[1,4] = Beryllium.ia_H_jb(0,2,1,5)
A[2,3] = Beryllium.ia_H_jb(1,3,0,4)
A[2,4] = Beryllium.ia_H_jb(1,3,1,5)

# <13|H|21>
A[3,1] = Beryllium.ia_H_jb(0,4,0,2)
A[3,2] = Beryllium.ia_H_jb(1,5,0,2)
A[4,1] = Beryllium.ia_H_jb(0,4,1,3)
A[4,2] = Beryllium.ia_H_jb(1,5,1,3)

# <13|H|31>
A[3,3] = Beryllium.ia_H_jb(0,4,0,4)
A[3,4] = Beryllium.ia_H_jb(0,4,1,5)
A[4,3] = Beryllium.ia_H_jb(1,5,0,4)
A[4,4] = Beryllium.ia_H_jb(1,5,1,5)


eigvals, eigvecs = np.linalg.eigh(A)

print(A)
print(eigvals)
