import numpy as np
from Matrix_elements import *

def h_HF_elements(i,j, C, Z, S=2):
    '''Hartree-Fock matrix elements (a,b) given C and Z
    i<3 => i_s=0, j<3 => j_s=0
    '''
    
    if i<3:
        i_s = 0
    else: 
        i = i-3
        i_s = 1
        
    if j<3:
        j_s = 0
    else:
        j = j-3
        j_s = 1
    
    C_conj = np.conj(C)
    v = OBME(Z)
    
    Result = v[i, j]
    
    for p in range(3):
        for p_s in range(S):
            for c in range(3):
                for c_s in range(S):
                    for d in range(3):
                        for d_s in range(S):
                            Result += C_conj[S*p+p_s,S*c+c_s]*C[S*p+p_s,S*d+d_s]*s2r_antisym(Z, i,c,j,d, i_s,c_s,j_s,d_s)
                            
    return Result
    
    
def h_HF_matrix(C, Z):
    '''HF-matrix'''
    
    HF = np.empty((6,6))
    
    for i in range(6):
        for j in range(6):
            HF[i,j] = h_HF_elements(i,j,C,Z)
    
    return HF
   
   
def h_HF_iter(Z):
    '''Solving HF with an iterative scheme'''
    
    C = np.eye(6)
    eigvals = np.ones(6)
    
    print(eigvals)
    
    while abs(eigvals[0]) > 0.001:
        HF = h_HF_matrix(C, Z)
        eigvals, C = np.linalg.eigh(HF)
        print(eigvals)
        print(C)
    
    return C
    

if __name__ == '__main__':
    h_HF_iter(2)
