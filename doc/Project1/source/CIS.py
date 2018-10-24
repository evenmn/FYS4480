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
