import numpy as np
from collections import Counter

def Cintegrate(phi,
               HC, 
               dt):
    """ 
    Explicit Euler integration to simulate a classical random walker.
    
    Parameters: 
    phi (array): prob. distr. at current time
    HC (sparse matrix): classical Hamiltonian
    dt (float): time step
  
    Returns: 
    sparse matrix: prob. distr. at next time step
  
    """

    phip1 = phi-HC.dot(phi)*dt
    
    return phip1

def Qintegrate(phi, 
               phim1, 
               HQ, 
               dt,
               eps = 0):
    """ 
    Integration scheme according to: Askar, Attila, and Ahmet S. Cakmak. 
    "Explicit integration method for the time‐dependent Schrodinger equation for collision problems." 
    The Journal of Chemical Physics 68.6 (1978): 2794-2798. 
    
    Parameters: 
    phi (array): wave function at current time
    phim1 (array): wave function at previous time
    HQ (sparse matrix): quantum Hamiltonian
    dt (float): time step
    eps (float): interpolation parameter (eps = 0 --> quantum, eps = 1 --> classical)

    Returns: 
    sparse matrix: wave function at next time step
  
    """
    phip1 = -2*1.0j*(1-eps)*dt*HQ.dot(phi)+phim1 
    
    return phip1

def jump_process(phi, 
                 dt, 
                 n, 
                 m):
    """ 
    Jump process to simulate classical dynamics in quantum-stochastic walks.
    
    Parameters: 
    phi (array):  wave function at current time
    dt (float): time step
    n,m (integer): jump process numbers
  
    Returns: 
    new phi vector
  
    """
    
    if n != m:
        phi_new = np.zeros(len(phi), dtype = complex)    
    
        phi_new[n] = np.copy(phi[m]/abs(phi[m]))
    
    else:
        phi_new = np.zeros(len(phi), dtype = complex)    
    
        phi_new[n] = np.copy(phi[m]/abs(phi[m]))*np.exp( 1.0j * np.pi * 0.5 )

    return phi_new

def select_jump(HC, 
                HC_row, 
                HC_col, 
                phi, 
                phim1, 
                dt, 
                eps, 
                Pjump):
    """ 
    Determine jump process in quantum-stochastic walks.
    
    Parameters: 
    HC (sparse matrix): classical Hamiltonian
    HC_row (array): indices of occupied rows in classical Hamiltonian 
    HC_col (array): indices of occupied columns in classical Hamiltonian
    phi (array): wave function at current time
    phi (array):  wave function at previous time
    dt (float): time step
    eps (float): interpolation parameter (eps = 0 --> quantum, eps = 1 --> classical)
    Pjump (float): jump probability (Pjump = 2*eps*dt)

    Returns: 
    new phi arrays
  
    """
     
    Psigma_cum = 0

    rnd = np.random.rand()
             
    for n,m in zip(HC_row,HC_col):

        if n != m:
            Psigma_cum -= dt * eps * HC[n,m] * abs(phi[m])**2  
        
        else:
            Psigma_cum += dt * eps * HC[n,m] * abs(phi[m])**2 
 
        if rnd < Psigma_cum/Pjump:
            phip1Q = jump_process(phi, dt, n, m)          
            
            phim1Q = np.copy(phip1Q)
            phiQ = np.copy(phip1Q)   
                        
            return phim1Q, phiQ
        
def simulate_qc(G, 
                HC, 
                HQ, 
                phim1Q0, 
                phiQ0, 
                phiC0, 
                N, 
                dt, 
                r, 
                average_threshold = 100):
    """ 
    Run simulation.
    
    Parameters: 
    G (graph): networkx graph object 
    HC (sparse matrix): classical Hamiltonian
    HQ (sparse matrix): quantum Hamiltonian
    phim1Q0 (array): initial wave function (t0)
    phiQ0 (array): initial wave function (t0 + dt)
    phiC0 (array): initial probability vector
    N (int): number of time steps
    dt (float): time step
    r (float): reset rate
    average_threshold (int): number of initial steps before averaging

    Returns: 
    degree array, occupation probabilities
  
    """
    
    assert N > average_threshold, "N should be larger than average_threshold"
    
    phim1Q = np.copy(phim1Q0)            
    phiQ = np.copy(phiQ0)        
    phiC = np.copy(phiC0)        

    nodes = sorted(G.nodes())
    degrees = G.degree(nodes)

    degree_sequence = sorted([d for n, d in degrees], reverse=True)
    degreeCount = Counter(degree_sequence)
    
    # dictionaries to store probabilities of walkers to occupy nodes
    # with certain degrees
    degreeProbDictQ = {x : 0 for x in degreeCount}
    degreeProbDictC = {x : 0 for x in degreeCount}
    
    number_samples = 0
    
    for i in range(N):
        
        eps = np.random.rand()
        
        # reset
        if eps < r*dt:
            
            # reset quantum walker
            phim1Q = np.copy(phim1Q0)  
            phiQ = np.copy(phiQ0)     
            
            # reset classical walker
            phiC = np.copy(phiC0) 
    
        # walker evolution according to Hamiltonian dynamics
        else:
            
            # evolution of quantum walker
            phip1Q = Qintegrate(phiQ, phim1Q, HQ, dt)
            phim1Q = np.copy(phiQ)
            phiQ = np.copy(phip1Q)      
                        
            # evolution of classical walker
            phip1C = Cintegrate(phiC, HC, dt)
            phiC = np.copy(phip1C)  
            
            if i > average_treshold:
                
                number_samples += 1
                
                for j in range(len(nodes)):
                    degreeProbDictQ[G.degree(j)] += abs(phiQ[j])**2  
                    degreeProbDictC[G.degree(j)] += phiC[j]
            
    for x in degreeCount:
        degreeProbDictQ[x] *= 1./number_samples
        degreeProbDictC[x] *= 1./number_samples
    
    deg, cnt = zip(*degreeCount.items())

    cnt = np.asarray(cnt)

    classical_observation_prob = \
    np.asarray([x for x in degreeProbDictC.values()])/cnt*1e3
    
    quantum_observation_prob = \
    np.asarray([x for x in degreeProbDictQ.values()])/cnt*1e3
    
    return deg, classical_observation_prob, quantum_observation_prob

def simulate_QSW(G, 
                 HC, 
                 HQ, 
                 phim10, 
                 phi0, 
                 N, 
                 dt, 
                 r, 
                 samples,
                 eps):
    """ 
    Run simulation.
    
    Parameters: 
    G (graph): networkx graph object 
    HC (sparse matrix): classical Hamiltonian
    HQ (sparse matrix): quantum Hamiltonian
    phim10 (array): initial wave function (t0)
    phi0 (array): initial wave function (t0 + dt)
    N (int): number of time steps
    dt (float): time step
    r (float): reset rate
    samples (int): number of samples
    eps (float): interpolation parameter (eps = 0 --> quantum, eps = 1 --> classical)

    Returns: 
    degree array, occupation probabilities
  
    """
    
    # jump probability
    Pjump = 2 * eps * dt    
    
    # reset probability
    Pr = r * dt
        
    phim1 = np.copy(phim10)            
    phi = np.copy(phi0)        
    
    HC_row = HC.tocoo().row
    HC_col = HC.tocoo().col 
    
    nodes = sorted(G.nodes())
    degrees = G.degree(nodes)

    degree_sequence = sorted([d for n, d in degrees], reverse=True)
    degreeCount = Counter(degree_sequence)
    
    # dictionary to store probabilities of walkers to occupy nodes
    # with certain degrees
    degreeProbDictQSW = {x : 0 for x in degreeCount}
    
    for sample in range(samples):
        
        for i in range(N):
            
            rnd = np.random.rand()
            
            # jump process
            if rnd < Pjump+Pr: 
                
                rnd = np.random.rand()
                
                # reset
                if rnd < Pr/(Pr+Pjump):
                    phim1 = np.copy(phim10)  
                    phi = np.copy(phi0)  
                
                # classical diffusion
                else:
                    phim1, phi = select_jump(HC, HC_row, HC_col, phi, \
                                             phim1, dt, eps, Pjump)
                        
            # no jump process
            else:
                phip1 = Qintegrate(phi, phim1, HQ, dt, eps)
        
                phim1 = np.copy(phi)
                phi = np.copy(phip1)                
                
        for j in range(len(phi)):
            degreeProbDictQSW[G.degree(j)] += phi[j]**2  
        
    for x in degreeCount:
        degreeProbDictQSW[x] *= 1./samples

    degreeCount = Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    cnt = np.asarray(cnt)
    
    QSW_observation_prob = \
    np.asarray([x for x in degreeProbDictQSW.values()])/cnt*1e3
    
    return deg, QSW_observation_prob

