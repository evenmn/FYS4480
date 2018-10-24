import numpy as np
from Matrix_elements import *

class CIS:
    '''Configuration Interction Singles class for atomic structure'''
    
    def __init__(self, Z, basis):
        '''
        Arguments:
        ----------
        Z:      Int.
                Atomic number (proton number)
        '''
        
        self.Elements = Integrals(Z,basis)  # Matrix elements/integrals
        self.n = Z                          # Assume neutral atom


    def c_H_c(self):
        '''Reference energy'''
        
        OBT = 0
        for i in range(self.n):
            OBT += self.Elements.OBME(i,i)
                
        TBT = 0
        for i in range(self.n):
            for j in range(self.n):
                TBT += 0.5*self.Elements.Antisym(i,j,i,j)
        
        return OBT + TBT
        
        
    def c_H_ia(self, i,a):
        '''Single excited ket'''
        
        OBT = self.Elements.OBME(i,a)
        
        TBT = 0
        for j in range(self.n):
            TBT += self.Elements.Antisym(a,j,i,j)
        
        return OBT + TBT
        
        
    def ia_H_jb(self, i,a,j,b):
        '''Single excited bra and ket'''
        
        Result = self.Elements.Antisym(a,j,i,b)
        
        if a==b:
            Result -= self.Elements.OBME(i,j)
            for k in range(self.n):
                Result -= self.Elements.Antisym(i,k,j,k)
                    
            if i==j:
                for k in range(self.n):
                    Result += self.Elements.OBME(k,k)
                    for l in range(self.n):
                        Result += 0.5*self.Elements.Antisym(k,l,k,l)
                                
        if i==j:
            Result += self.Elements.OBME(a,b)
            for k in range(self.n):
                Result += self.Elements.Antisym(a,k,b,k)
        
        return Result
        
        
    def organize(self):
        '''Set up matrix
        NEEDS TO BE GENERALIZED AND SIMPLIFIED'''
    
        A = np.zeros((5,5))

        # --- <c|H|c> ---
        A[0,0] = CIS.c_H_c(self)

        # --- <c|H|p_i^a> ---
        A[0,1] = CIS.c_H_ia(self,0,2)
        A[0,2] = CIS.c_H_ia(self,1,3)
        A[0,3] = CIS.c_H_ia(self,0,4)
        A[0,4] = CIS.c_H_ia(self,1,5)

        # --- <p_i^a|H|c> ---
        A[1,0] = CIS.c_H_ia(self,0,2)
        A[2,0] = CIS.c_H_ia(self,1,3)
        A[3,0] = CIS.c_H_ia(self,0,4)
        A[4,0] = CIS.c_H_ia(self,1,5)

        # --- <p_i^a|H|p_j^b> ---
        # <12|H|21>
        A[1,1] = CIS.ia_H_jb(self,0,2,0,2)
        A[1,2] = CIS.ia_H_jb(self,0,2,1,3)
        A[2,1] = CIS.ia_H_jb(self,1,3,0,2)
        A[2,2] = CIS.ia_H_jb(self,1,3,1,3)

        # <12|H|31>
        A[1,3] = CIS.ia_H_jb(self,0,2,0,4)
        A[1,4] = CIS.ia_H_jb(self,0,2,1,5)
        A[2,3] = CIS.ia_H_jb(self,1,3,0,4)
        A[2,4] = CIS.ia_H_jb(self,1,3,1,5)

        # <13|H|21>
        A[3,1] = CIS.ia_H_jb(self,0,4,0,2)
        A[3,2] = CIS.ia_H_jb(self,1,5,0,2)
        A[4,1] = CIS.ia_H_jb(self,0,4,1,3)
        A[4,2] = CIS.ia_H_jb(self,1,5,1,3)

        # <13|H|31>
        A[3,3] = CIS.ia_H_jb(self,0,4,0,4)
        A[3,4] = CIS.ia_H_jb(self,0,4,1,5)
        A[4,3] = CIS.ia_H_jb(self,1,5,0,4)
        A[4,4] = CIS.ia_H_jb(self,1,5,1,5)
        
        return A
