import numpy as np
from CIS import *

Beryllium = CIS(Z=4)
A = np.zeros((5,5))

# --- <c|H|c> ---
A[0,0] = Beryllium.c_H_c()

# --- <c|H|p_i^a> ---
A[0,1] = Beryllium.c_H_ia(0,1,0,0)
A[0,2] = Beryllium.c_H_ia(0,1,1,1)
A[0,3] = Beryllium.c_H_ia(0,2,0,0)
A[0,4] = Beryllium.c_H_ia(0,2,1,1)

# --- <p_i^a|H|c> ---
A[1,0] = Beryllium.c_H_ia(0,1,0,0)
A[2,0] = Beryllium.c_H_ia(0,1,1,1)
A[3,0] = Beryllium.c_H_ia(0,2,0,0)
A[4,0] = Beryllium.c_H_ia(0,2,1,1)

# --- <p_i^a|H|p_j^b> ---
# <12|H|21>
A[1,1] = Beryllium.ia_H_jb(0,1,0,1,0,0,0,0)
A[1,2] = Beryllium.ia_H_jb(0,1,0,1,0,0,1,1)
A[2,1] = Beryllium.ia_H_jb(0,1,0,1,1,1,0,0)
A[2,2] = Beryllium.ia_H_jb(0,1,0,1,1,1,1,1)

# <12|H|31>
A[1,3] = Beryllium.ia_H_jb(0,1,0,2,0,0,0,0)
A[1,4] = Beryllium.ia_H_jb(0,1,0,2,0,0,1,1)
A[2,3] = Beryllium.ia_H_jb(0,1,0,2,1,1,0,0)
A[2,4] = Beryllium.ia_H_jb(0,1,0,2,1,1,1,1)

# <13|H|21>
A[3,1] = Beryllium.ia_H_jb(0,2,0,1,0,0,0,0)
A[3,2] = Beryllium.ia_H_jb(0,2,0,1,1,1,0,0)
A[4,1] = Beryllium.ia_H_jb(0,2,0,1,0,0,1,1)
A[4,2] = Beryllium.ia_H_jb(0,2,0,1,1,1,1,1)

# <13|H|31>
A[3,3] = Beryllium.ia_H_jb(0,2,0,2,0,0,0,0)
A[3,4] = Beryllium.ia_H_jb(0,2,0,2,1,1,0,0)
A[4,3] = Beryllium.ia_H_jb(0,2,0,2,0,0,1,1)
A[4,4] = Beryllium.ia_H_jb(0,2,0,2,1,1,1,1)


eigvals, eigvecs = np.linalg.eigh(A)

print(A)
print(eigvals)
