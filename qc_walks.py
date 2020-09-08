import numpy as np
import networkx as nx
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from walkerlib.hamiltonians import CHamiltonian, QHamiltonian
from walkerlib.integration import simulate_qc
from walkerlib.analytical import stationary_distribution, long_time_average

# linear dimension lattice
L = 40
G = nx.erdos_renyi_graph(L**2, p=0.005, seed = 123)

# quantum Hamiltonian
HQ = QHamiltonian(G)
    
# classical Hamiltonian
HC = CHamiltonian(G)

# time step and simulation time
dt = 0.01
T = 20
N = int(T/dt)

# reset rate
r = 0.3

# initializing wave function
phim1Q0 = np.ones(len(G))
phim1Q0 /= np.linalg.norm(phim1Q0)    
phiQ0 = phim1Q0-HQ.dot(phim1Q0)*dt

# initializing probability distribution
phiC0 = np.ones(len(G))
phiC0 /= sum(phiC0) 

start_time = time.time()   

# simulate classical and quantum walker
degree_arr, classical_observation_prob, \
            quantum_observation_prob = simulate_qc(G, 
                                                   HC, 
                                                   HQ, 
                                                   phim1Q0, 
                                                   phiQ0, 
                                                   phiC0, 
                                                   N, 
                                                   dt, 
                                                   r, 
                                                   average_threshold = 100)

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
plt.plot(degree_arr, quantum_observation_prob, marker = 'o', ls = 'None', \
         markersize = 4, label=r'quantum, $q^\prime_k(r)$')
plt.plot(degree_arr, classical_observation_prob, marker = 'x', ls = 'None', \
         markersize = 4, label=r'classical, $p^\prime_k(r)$')
plt.xlim(0,20)
plt.ylim(0,2)  
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
plt.legend(loc = 2, frameon = False, fontsize = 8)
plt.xlabel(r'degree')
plt.ylabel(r'occupation prob. ($10^{-3}$)')
plt.tight_layout()
plt.show()
