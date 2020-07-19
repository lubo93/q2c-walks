import numpy as np
import networkx as nx
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from walkerlib.hamiltonians import CHamiltonian, QHamiltonian
from walkerlib.integration import simulate_QSW
from walkerlib.analytical import stationary_distribution, long_time_average

# linear dimension lattice
L = 10
G = nx.erdos_renyi_graph(L**2, p=0.1, seed = 123)

# quantum Hamiltonian
HQ = QHamiltonian(G)
    
# classical Hamiltonian
HC = CHamiltonian(G)

# time step and simulation time
dt = 0.01
T = 10
N = int(T/dt)

# reset rate
r = 0.3

# initializing wave function
phim1Q0 = np.ones(len(G))
phim1Q0 /= np.linalg.norm(phim1Q0)    
phiQ0 = np.copy(phim1Q0)

# initializing probability distribution
phiC0 = np.ones(len(G))
phiC0 /= sum(phiC0) 

# classical to quantum interpolation parameter
eps = 0.1

start_time = time.time()   

# simulate quantum stochastic walk
degree_arr, QSW_observation_prob = simulate_QSW(G, 
                                                HC, 
                                                HQ, 
                                                phim1Q0, 
                                                phiQ0, 
                                                N, 
                                                dt, 
                                                r, 
                                                samples = 100,
                                                eps = eps)

print("--- %s seconds runtime ---" % (time.time() - start_time))

degree_arr, classical_observation_prob_analytical = \
                        stationary_distribution(G, HC, phiC0, r)
                        
degree_arr, quantum_observation_prob_analytical = \
                        long_time_average(G, HQ, phim1Q0, r)

fig, ax = plt.subplots()
plt.title(r'reset rate $r=%1.1f$'%r)
plt.plot(degree_arr, quantum_observation_prob_analytical, marker = 'd', \
         ls = 'None', color = 'grey')
plt.plot(degree_arr, classical_observation_prob_analytical, marker = 's', \
         ls = 'None', color = 'grey')
plt.plot(degree_arr, QSW_observation_prob, marker = 'o', ls = 'None', \
         markersize = 4, label=r'QSW, $\epsilon=%1.1f$'%eps)
plt.xlim(0,20)
plt.ylim(0,2)  
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
plt.legend(loc = 2, frameon = False, fontsize = 8)
plt.xlabel(r'degree')
plt.ylabel(r'occupation prob. ($10^{-3}$)')
plt.tight_layout()
plt.show()