from numpy import asarray, newaxis, kron
from scipy.sparse import diags
from networkx import laplacian_matrix, normalized_laplacian_matrix

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
