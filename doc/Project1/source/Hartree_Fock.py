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
    
    v = OBME(Z)
    
    Result = v[i, j]
    
    for p in range(2):
      for p_s in range(S):
        for c in range(2):
          for c_s in range(S):
            for d in range(2):
              for d_s in range(S):
                Result += C[S*p+p_s,S*c+c_s]*C[S*p+p_s,S*d+d_s]*s2r_antisym(Z, i,c,j,d, i_s,c_s,j_s,d_s)
                            
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
    
    #while abs(eigvals[0]) > 0.001:
    for i in range(10):
        HF = h_HF_matrix(C, Z)
        eigvals, C = np.linalg.eigh(HF)
        print(eigvals)
        print(C)
    
    return C
    
    
def calc_E(C, Z, S=2):
    '''Calculate energy'''
    
    v = OBME(Z)
    
    E = 0
    for p in range(2):
      for a in range(3):
        for a_s in range(S):
          for b in range(3):
            for b_s in range(S):
              E += C[p,a]*C[p,b]*v[a,b]
              for q in range(2):
                for c in range(3):
                  for c_s in range(S):
                    for d in range(3):
                      for d_s in range(S):
                        E += C[p,S*a+a_s]*C[q,S*b+b_s]*C[p,S*c+c_s]*C[q,S*d+d_s]*s2r_antisym(Z, a,b,c,d, a_s,b_s,c_s,d_s)
                        
    return E
    

if __name__ == '__main__':
    C = h_HF_iter(2)
    print(calc_E(C, Z=2))
