import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian_elements import ground_state
from matrix_elements import TBME, OBME

N = 10

GS_energies = np.empty(N)

for Z in range(1,N+1):
    u = TBME(Z)
    v = OBME(Z)
    
    print(v[0,0], v[1,1], v[2,2])

    GS_energies[Z-1] = ground_state(u, v)
    
print(GS_energies)
plt.plot(GS_energies)
plt.show()
    
