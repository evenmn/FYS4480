import numpy as np
from matrix_elements import TBME, OBME

def single_excitations(a, b, u, v):
    ''' Calculate <p_1^a|H|p_1^b> '''
    
    return 0.5*(v[a, b] - v[0, 0] + u[a,0,0,b] - u[a,0,b,0])
    
    
def ground_state(u, v):
    ''' Calculate <c|H|c> '''
    
    OBT = 0
    TBT = 0
    
    for i in range(3):
        OBT += v[i,i]
        for j in range(3):
            TBT += u[i,j,i,j] - u[i,j,j,i]
    return OBT + 0.5*TBT
    
    
    
if __name__ == '__main__':
    u = TBME(Z=1)
    v = OBME(Z=1)
    
    A = np.zeros((5,5))
    A[0,0] = ground_state(u, v)
    
    # <12|H|21>
    A[1,1] = single_excitations(1,1,u,v)
    A[1,2] = single_excitations(1,1,u,v)
    A[2,1] = single_excitations(1,1,u,v)
    A[2,2] = single_excitations(1,1,u,v)
    
    # <13|H|21>
    A[1,3] = single_excitations(1,2,u,v)
    A[1,4] = single_excitations(1,2,u,v)
    A[2,3] = single_excitations(1,2,u,v)
    A[2,4] = single_excitations(1,2,u,v)
    
    # <12|H|31>
    A[3,1] = single_excitations(2,1,u,v)
    A[3,2] = single_excitations(2,1,u,v)
    A[4,1] = single_excitations(2,1,u,v)
    A[4,2] = single_excitations(2,1,u,v)
    
    # <13|H|31>
    A[3,3] = single_excitations(2,2,u,v)
    A[3,4] = single_excitations(2,2,u,v)
    A[4,3] = single_excitations(2,2,u,v)
    A[4,4] = single_excitations(2,2,u,v)
    
    print(A)
    
    eigvals, eigvecs = np.linalg.eigh(A)
    print(eigvecs)
    
