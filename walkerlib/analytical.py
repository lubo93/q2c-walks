import numpy as np
import time
from collections import Counter
from scipy.sparse.linalg import spsolve
from scipy.sparse import identity

def stationary_distribution(G, HC, phiC0, r):
    """ 
    Stationary occupation probability of a classical random walker on a graph G. 
    
    Parameters: 
    G (graph): networkx graph object 
    HC (sparse matrix): classical Hamiltonian
    phiC0 (array): initial probability vector
    r (float): reset rate

    Returns: 
    arrays: degrees, stationary occupation probability distribution
  
    """
    
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    degreeCount = Counter(degree_sequence)
    degreeProbDictC_analtyic = {x : 0 for x in degreeCount} 
    deg, cnt = zip(*degreeCount.items())

    stationary_distr_cl = \
    spsolve(1/(r+1e-12)*(r*identity(HC.shape[0])+HC), phiC0)
    
    degreeProbDictC_analtyic = {x : 0 for x in degreeCount} 
    
    for j in range(stationary_distr_cl.size):
        degreeProbDictC_analtyic[G.degree(j)] += stationary_distr_cl[j]   
    
    classical_observation_prob_analytical = \
    np.asarray([x for x in degreeProbDictC_analtyic.values()])/cnt*sum(cnt)
    
    return deg, classical_observation_prob_analytical
        

def long_time_average(G, HQ, phim1Q0, r):
    """ 
    Stationary occupation probability of a quantum random walker on a graph G. 
    
    Parameters: 
    G (graph): networkx graph object 
    HQ (sparse matrix): quantum Hamiltonian
    phim1Q0 (array): initial wave function
    r (float): reset rate

    Returns: 
    arrays: degrees, stationary occupation probability distribution
  
    """
    
    eigenvalues, eigenvectors = np.linalg.eigh(HQ.todense())

    stationary_distr_q = np.zeros(phim1Q0.size, dtype = complex)

    start_time = time.time()   

    print("--- long time average computation started ---")
    
    for k in range(eigenvalues.size):
        print("--- %d/%d eigenvalues completed ---"%(k,eigenvalues.size))
        for kp in range(len(eigenvalues)):
            
            v1 = np.ravel(eigenvectors[:,k])
            v2 = np.ravel(eigenvectors[:,kp])
            
            prod1 = np.dot(phim1Q0, v1)
            prod2 = np.dot(phim1Q0, v2)
            
            stationary_distr_q += \
            v1*v2*r/(r+1.0j*(eigenvalues[k]-eigenvalues[kp]))*prod1*prod2
     
    print("--- long time average computation completed after %s ---" \
          % (time.time() - start_time))
           
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    degreeCount = Counter(degree_sequence)
    degreeProbDictQ_analtyic = {x : 0 for x in degreeCount} 
    deg, cnt = zip(*degreeCount.items())

    for j in range(stationary_distr_q.size):
        degreeProbDictQ_analtyic[G.degree(j)] += stationary_distr_q[j]   
    
    quantum_observation_prob_analytical = \
    np.asarray([float(x) for x in degreeProbDictQ_analtyic.values()])/cnt*sum(cnt)
        
    return deg, quantum_observation_prob_analytical