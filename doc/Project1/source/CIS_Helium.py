import numpy as np
from matrix_elements import TBME, OBME
    
    
u = TBME(Z=2)
v = OBME(Z=2)

A = np.zeros((5,5))
A[0,0] = ground_state(u,v)

# <12|H|21>
A[1,1] = single_excitations(1,1,u,v)
A[1,2] = single_excitations(1,1,u,v)
A[2,1] = single_excitations(1,1,u,v)
A[2,2] = single_excitations(1,1,u,v)

# <12|H|31>
A[1,3] = single_excitations(1,2,u,v)
A[1,4] = single_excitations(1,2,u,v)
A[2,3] = single_excitations(1,2,u,v)
A[2,4] = single_excitations(1,2,u,v)

# <13|H|21>
A[3,1] = single_excitations(2,1,u,v)
A[3,2] = single_excitations(2,1,u,v)
A[4,1] = single_excitations(2,1,u,v)
A[4,2] = single_excitations(2,1,u,v)

# <13|H|31>
A[3,3] = single_excitations(2,2,u,v)
A[3,4] = single_excitations(2,2,u,v)
A[4,3] = single_excitations(2,2,u,v)
A[4,4] = single_excitations(2,2,u,v)


eigvals, eigvecs = np.linalg.eigh(A)

print(A)
print(eigvals)
    
