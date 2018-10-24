import numpy as np
from Matrix_elements import *

class CIS:
    '''Configuration Interction Singles class for atomic structure'''
    
    def __init__(self, Z, S=2):
        '''
        Arguments:
        ----------
        Z:      Int.
                Atomic number (proton number)
        S:      Int.
                Number of spin states, should be fixed to 2
        '''
        
        self.S = S
        self.v = OBME(Z)
        self.Z = Z
        self.Orbits = int(Z/2)          # Assume fully occupied GS orbits


    def c_H_c(self):
        '''Reference energy'''
        
        v = self.v
        S = self.S
        Z = self.Z
        Orbits = self.Orbits
        
        OBT = 0
        for i in range(Orbits):
            for i_s in range(S):
                OBT += v[i,i]
                
        TBT = 0
        for i in range(Orbits):
            for j in range(Orbits):
                for i_s in range(S):
                    for j_s in range(S):
                        TBT += 0.5*s2r_antisym(Z,i,j,i,j, i_s,j_s,i_s,j_s)
        
        return OBT + TBT
        
        
    def c_H_ia(self, i,a,i_s,a_s):
        '''hggh'''
        
        v = self.v
        S = self.S
        Z = self.Z
        Orbits = self.Orbits
        
        OBT = v[i,a]
        
        TBT = 0
        for j in range(Orbits):
            for j_s in range(S):
                TBT += s2r_antisym(Z,a,j,i,j, a_s,j_s,i_s,j_s)
        
        return OBT + TBT
        
        
    def ia_H_jb(self, i,a,j,b,i_s,a_s,j_s,b_s):
        '''hhh'''
        
        v = self.v
        S = self.S
        Z = self.Z
        Orbits = self.Orbits
        
        '''
        Result = s2r_antisym(Z,a,j,i,b, a_s,j_s,i_s,b_s)
        
        if i==j and i_s==j_s:
            Result -= v[a,b]
        if a==b and a_s==b_s:
            Result += v[i,j]
        
        '''
        Result = s2r(Z,a,j,i,b, a_s,j_s,i_s,b_s)
        
        if a==b and a_s==b_s:
            Result -= v[i,j]
            for k in range(Orbits):
                for k_s in range(S):
                    Result -= s2r_antisym(Z,i,k,j,k, i_s,k_s,j_s,k_s)
                    
            if i==j and i_s==j_s:
                for k in range(Orbits):
                    for k_s in range(S):
                        Result += v[k,k]
                        for l in range(Orbits):
                            for l_s in range(S):
                                Result += 0.5*s2r_antisym(Z,k,l,k,l, k_s,l_s,k_s,l_s)
                                
        if i==j and i_s==j_s:
            Result += v[a,b]
            for k in range(Orbits):
                for k_s in range(S):
                    Result += s2r_antisym(Z,a,k,b,k, a_s,k_s,b_s,k_s)
        
        
        return Result
