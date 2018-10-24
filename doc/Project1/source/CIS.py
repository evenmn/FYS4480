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
        '''hggh'''
        
        OBT = self.Elements.OBME(i,a)
        
        TBT = 0
        for j in range(self.n):
            TBT += self.Elements.Antisym(a,j,i,j)
        
        return OBT + TBT
        
        
    def ia_H_jb(self, i,a,j,b):
        '''hhh'''
        
        
        '''
        Result = s2r_antisym(Z,a,j,i,b, a_s,j_s,i_s,b_s)
        
        if i==j and i_s==j_s:
            Result -= v[a,b]
        if a==b and a_s==b_s:
            Result += v[i,j]
        
        '''
        Result = self.Elements.TBME(a,j,i,b)
        
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
        
if __name__ == '__main__':
    basis = ((0,0), (0,1), (1,0), (1,1), (2,0), (2,1))
    
    Helium = CIS(2,basis)
    
    print(Helium.c_H_c())
    print(Helium.c_H_ia(0,4))
