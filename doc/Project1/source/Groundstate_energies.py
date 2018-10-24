import numpy as np
import matplotlib.pyplot as plt
from Hamiltonian_elements import ground_state
from matrix_elements import TBME, OBME

N = 4

GS_energies = np.empty(N)

for Z in range(1,N+1):
    u = TBME(Z)
    v = OBME(Z)
    

    GS_energies[Z-1] = ground_state(u, v, Z)
    
print(GS_energies)
#plt.plot(GS_energies)
#plt.show()

plt.axhline(GS_energies[0], color='r')
plt.axhline(-0.4999426, linestyle='--', color='r')
plt.axhline(GS_energies[1], color='b')
plt.axhline(-2.903694, linestyle='--', color='b')
plt.axhline(GS_energies[2], color='g')
plt.axhline(-7.279, linestyle='--', color='g')
plt.axhline(GS_energies[3], color='k')
plt.axhline(-14.6674, linestyle='--', color='k')
plt.show()
    
