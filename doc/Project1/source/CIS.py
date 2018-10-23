import numpy as np
from matrix_elements import *

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
        self.u = TBME(Z)
        self.v = OBME(Z)
        self.Orbits = int(Z/2)          # Assume fully occupied GS orbits
        
        
    def s2r(self, p,q,r,s, p_s,q_s,r_s,s_s):
        ''' Returning value of matrix element with given 
        quantum number p,q,r,s and spin number p_s,q_s,r_s,s_s.
        Spin up/down is 0/1'''
        
        u = self.u
        
        if p_s==r_s and q_s==s_s:
            return u[p,q,r,s]
        else:
            return 0
        

    def s2r_antisym(self,p,q,r,s, p_s,q_s,r_s,s_s):
        ''' Returning value of matrix element with given 
        quantum number p,q,r,s and spin number p_s,q_s,r_s,s_s.
        Spin up/down is 0/1.
        Element is now antisymmetric, such that
        <pq|H|rs>AS=<pq|H|rs>-<pq|H|sr>'''
        
        return CIS.s2r(self,p,q,r,s, p_s,q_s,r_s,s_s)-CIS.s2r(self,p,q,s,r, p_s,q_s,s_s,r_s)
    


    def c_H_c(self):
        '''Reference energy'''
        
        u = self.u
        v = self.v
        S = self.S
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
                        TBT += 0.5*CIS.s2r_antisym(self,i,j,i,j, i_s,j_s,i_s,j_s)
        
        return OBT + TBT
        
        
    def c_H_ia(self, i,a,i_s,a_s):
        '''hggh'''
        
        u = self.u
        v = self.v
        S = self.S
        Orbits = self.Orbits
        
        OBT = v[i,a]
        
        TBT = 0
        for j in range(Orbits):
            for j_s in range(S):
                TBT += CIS.s2r_antisym(self,a,j,i,j, a_s,j_s,i_s,j_s)
        
        return OBT + TBT
        
        
    def ia_H_jb(self, i,a,j,b,i_s,a_s,j_s,b_s):
        '''hhh'''
        
        u = self.u
        v = self.v
        S = self.S
        Orbits = self.Orbits
        
        '''
        Result = CIS.s2r_antisym(self,a,j,i,b, a_s,j_s,i_s,b_s)
        
        if i==j and i_s==j_s:
            Result -= v[a,b]
        if a==b and a_s==b_s:
            Result += v[i,j]
        
        '''
        Result = CIS.s2r(self,a,j,i,b, a_s,j_s,i_s,b_s)
        
        if a==b and a_s==b_s:
            Result -= v[i,j]
            for k in range(Orbits):
                for k_s in range(S):
                    Result -= CIS.s2r_antisym(self,i,k,j,k, i_s,k_s,j_s,k_s)
                    
            if i==j and i_s==j_s:
                for k in range(Orbits):
                    for k_s in range(S):
                        Result += v[k,k]
                        for l in range(Orbits):
                            for l_s in range(S):
                                Result += 0.5*CIS.s2r_antisym(self,k,l,k,l, k_s,l_s,k_s,l_s)
                                
        if i==j and i_s==j_s:
            Result += v[a,b]
            for k in range(Orbits):
                for k_s in range(S):
                    Result += CIS.s2r_antisym(self,a,k,b,k, a_s,k_s,b_s,k_s)
                    
        
        return Result
