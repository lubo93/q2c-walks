import numpy as np
import networkx as nx
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from walkerlib.hamiltonians import GroverHamiltonianSymmetric, GroverHamiltonianChilds
from walkerlib.integration import Qintegrate

# linear dimension lattice
L = 10
G = nx.complete_graph(100)#nx.erdos_renyi_graph(L**2, p=0.1, seed = 123)##

# quantum Hamiltonian
w = np.zeros(len(G))
w[0] = 1
gamma = 1/len(G)

HG, L = GroverHamiltonianChilds(G, gamma, w)

eigenvalues, eigenvectors = np.linalg.eigh(L.todense())

# initializing wave function
phim1Q0 = np.ravel(eigenvectors[:,0])
phim1Q0 /= np.linalg.norm(phim1Q0)   
 
#phiQ0 = phim1Q0-HG.dot(phim1Q0)*dt
gamma_arr = np.linspace(1e-9,2*gamma,100)
v0_arr = []
v1_arr = []

deltaE = []

for gamma in gamma_arr:
    HG, L = GroverHamiltonianChilds(G, gamma, w)
    eigenvalues, eigenvectors = np.linalg.eigh(HG.todense())
    
    deltaE.append(eigenvalues[1]-eigenvalues[0])
    v0 = np.ravel(eigenvectors[:,0])
    v1 = np.ravel(eigenvectors[:,1])
    
    v0_arr.append(v0)
    v1_arr.append(v1)

plt.figure()
plt.plot(gamma_arr, [np.dot(w,v0)**2 for v0 in v0_arr])
plt.plot(gamma_arr, [np.dot(w,v1)**2 for v1 in v1_arr], '--')
plt.plot(gamma_arr, [np.dot(phim1Q0,v0)**2 for v0 in v0_arr])
plt.plot(gamma_arr, [np.dot(phim1Q0,v1)**2 for v1 in v1_arr], '--')
plt.xlabel(r'$\gamma$')
plt.plot(gamma_arr, deltaE)
plt.show()
    
# time step and simulation time
dt = 0.01
T = 20
N = int(T/dt)

# reset rate
r = 0.0



start_time = time.time()   

# simulate Grover walk

phim1Q = np.copy(phim1Q0)
phiQ = np.copy(phiQ0)        
    
for i in range(N):
    # evolution of quantum walker
    phip1Q = Qintegrate(phiQ, phim1Q, HG, dt)
    phim1Q = np.copy(phiQ)
    phiQ = np.copy(phip1Q)
    
    plt.figure()
    plt.plot(np.conj(phiQ)*phiQ)
    plt.show()
    
print("--- %s seconds runtime ---" % (time.time() - start_time))

fig, ax = plt.subplots()
plt.tight_layout()
plt.show()