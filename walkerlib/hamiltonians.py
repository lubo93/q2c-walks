from numpy import asarray, newaxis
from scipy.sparse import diags, kron
from networkx import laplacian_matrix, normalized_laplacian_matrix
from networkx.linalg.graphmatrix import adjacency_matrix

def CHamiltonian(G):
    """ 
    Hamiltonian of a classical random walker on a graph G. 
    
    Parameters: 
    G (graph): networkx graph object 
  
    Returns: 
    sparse matrix: Hamiltonian
  
    """
    
    nodes = sorted(G.nodes())
    
    L = laplacian_matrix(G, nodelist=nodes)
    
    degrees = G.degree(nodes)
    degrees_inv = asarray([x[1]**-1.0 for x in degrees])
    
    Dm1 = diags(degrees_inv)
    
    HC = L.dot(Dm1)
    
    return HC

def QHamiltonian(G):
    """ 
    Hamiltonian of a quantum random walker on a graph G. 
    
    Parameters: 
    G (graph): networkx graph object 
  
    Returns: 
    sparse matrix: quantum Hamiltonian
  
    """
    nodes = sorted(G.nodes())
   
    HQ = normalized_laplacian_matrix(G, nodelist=nodes)
    
    return HQ

def GroverHamiltonianSymmetric(G, gamma, w):
    """ 
    Grover Hamiltonian of a quantum random walker on a graph G. 
    
    Parameters: 
    G (graph): networkx graph object 
    gamma (float): Laplacian prefactor 
    w (float): search (or target) state 

    Returns: 
    sparse matrix: quantum walk Hamiltonian for Grover search and
    normalized Laplacian
  
    """
    nodes = sorted(G.nodes())
    
    wT = w[:,newaxis]
    
    HG = gamma*normalized_laplacian_matrix(G, nodelist=nodes) + kron(wT, w) 
    
    return HG, normalized_laplacian_matrix(G, nodelist=nodes)

def GroverHamiltonianChilds(G, gamma, w):
    """ 
    Grover Hamiltonian of a quantum random walker on a graph G. 
    
    Parameters: 
    G (graph): networkx graph object 
    gamma (float): Laplacian prefactor 
    w (float): search (or target) state 

    For further information, see Childs, A. M., & Goldstone, J. (2004). 
    Spatial search by quantum walk. Physical Review A, 70(2), 022314.

    Returns: 
    sparse matrix: quantum walk Hamiltonian for Grover search and Laplacian
  
    """
    nodes = sorted(G.nodes())
    
    wT = w[:,newaxis]
    
    HG = gamma*laplacian_matrix(G, nodelist=nodes) - kron(wT, w) 
    
    return HG, laplacian_matrix(G, nodelist=nodes)
