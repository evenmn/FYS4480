import numpy as np
from CIS import *

Helium = CIS(Z=2)
A = np.zeros((5,5))

# --- <c|H|c> ---
A[0,0] = Helium.c_H_c()

# --- <c|H|p_i^a> ---
A[0,1] = Helium.c_H_ia(0,1,0,0)
A[0,2] = Helium.c_H_ia(0,1,1,1)
A[0,3] = Helium.c_H_ia(0,2,0,0)
A[0,4] = Helium.c_H_ia(0,2,1,1)

# --- <p_i^a|H|c> ---
A[1,0] = Helium.c_H_ia(0,1,0,0)
A[2,0] = Helium.c_H_ia(0,1,1,1)
A[3,0] = Helium.c_H_ia(0,2,0,0)
A[4,0] = Helium.c_H_ia(0,2,1,1)

# --- <p_i^a|H|p_j^b> ---
# <12|H|21>
A[1,1] = Helium.ia_H_jb(0,1,0,1,0,0,0,0)
A[1,2] = Helium.ia_H_jb(0,1,0,1,0,0,1,1)
A[2,1] = Helium.ia_H_jb(0,1,0,1,1,1,0,0)
A[2,2] = Helium.ia_H_jb(0,1,0,1,1,1,1,1)

# <12|H|31>
A[1,3] = Helium.ia_H_jb(0,1,0,2,0,0,0,0)
A[1,4] = Helium.ia_H_jb(0,1,0,2,0,0,1,1)
A[2,3] = Helium.ia_H_jb(0,1,0,2,1,1,0,0)
A[2,4] = Helium.ia_H_jb(0,1,0,2,1,1,1,1)

# <13|H|21>
A[3,1] = Helium.ia_H_jb(0,2,0,1,0,0,0,0)
A[3,2] = Helium.ia_H_jb(0,2,0,1,1,1,0,0)
A[4,1] = Helium.ia_H_jb(0,2,0,1,0,0,1,1)
A[4,2] = Helium.ia_H_jb(0,2,0,1,1,1,1,1)

# <13|H|31>
A[3,3] = Helium.ia_H_jb(0,2,0,2,0,0,0,0)
A[3,4] = Helium.ia_H_jb(0,2,0,2,1,1,0,0)
A[4,3] = Helium.ia_H_jb(0,2,0,2,0,0,1,1)
A[4,4] = Helium.ia_H_jb(0,2,0,2,1,1,1,1)


eigvals, eigvecs = np.linalg.eigh(A)

print(A)
print(eigvals)
