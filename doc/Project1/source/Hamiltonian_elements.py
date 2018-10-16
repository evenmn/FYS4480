import numpy as np

def single_excitations(a, b, u, v):
    ''' Calculate <p_1^a|H|p_1^b> '''
    
    return 0.5*(v[a, b] - v[0,0] + u[a,0,0,b] - u[a,0,b,0])
    
    
def ground_state(u, v):
    ''' Calculate <c|H|c> '''
    
    '''
    OBT = 0
    TBT = 0
    
    for i in range(3):
        OBT += v[i,i]
        for j in range(3):
            TBT += u[i,j,i,j] - u[i,j,j,i]
    '''
    
    OBT = v[0,0] + v[1,1] + v[2,2]
    TBT = u[0,1,0,1] - u[0,1,1,0] + u[0,2,0,2] - u[0,2,2,0] + u[1,2,1,2] - u[1,2,2,1]
    
    return OBT + TBT
