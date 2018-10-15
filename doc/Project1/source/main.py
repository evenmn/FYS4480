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
    
    A = np.zeros((3,3))
    
    A[0,0] = ground_state(u, v)
    for i in range(1,3):
        for j in range(1,3):
            A[i,j] = single_excitations(i, j, u, v)
    
    eigvals, eigvecs = np.linalg.eigh(A)
    print(eigvals)
    
