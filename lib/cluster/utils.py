#############################################################################
# utils.py
#
# Utility functions for the 'phylogeny-aware' clustering algorithm
#############################################################################
import os, sys
import numpy as np 
from scipy.stats import binom
import pickle 
import matplotlib.pyplot as plt
from omicsdata.ssm import constants, supervariants
from omicsdata.tree import adj, parents

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "metrics")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "orchard")))

from constants import KEYS_ORCH_NPZ
import neutree 

# CONSTANTS
LLH_KEY = "llhs"
CLUSTERS_KEY = "clusters"
PARENTS_KEY = "parents"
MIN_LINKAGE="min_linkage"
AVG_LINKAGE="avg_linkage"

EPSILON = np.exp(-30)
NULL_NODE = -1
INF = np.inf
BIC = "BIC"
GIC = "GIC"
model_selection_choices = {BIC:lambda x: np.log(x[0]), GIC:lambda x: np.log(x[0])*np.log(x[1])}


join_sort = lambda x: "".join(sorted([x[0],x[1]])) # lambda function for sorting variant ids

class PD_Dict:
    """Pairwise distance dictionary class
    
    Attributes
    ------------
    __dict : dict
        a python dictionary containing all of the pairwise distances between pairs of nodes
    """
    def __init__(self):
        self.__dict = {}

    def put(self, v1, v2, value):
        """Maps two keys v1, v2 to a single value"""
        self.__dict[join_sort((v1,v2))] = value

    def get(self, v1, v2):
        """Gets a value based on two keys v1, v2"""
        return self.__dict.get(join_sort((v1,v2)))

def cluster_llh(variants, cluster):
    """Computes the binomial likelihood a single cluster
    
    Parameters
    ------------
    variants : dict
        a dictionary where the keys are 'id' values and the values are python dictionaries containing the read
        count data for the mutation signified by a particular 'id'
    cluster : list
        a list of lists where each sublist is a particular cluster containing a set of unique 'id' values denoting
        the mutations in that cluster

    Returns
    ---------
    float
        the log-likelihood of all of the clusters under a binomial likelihood model
    """
    V_hat, N_hat, _, = prep_cluster_data(variants, cluster)
    P_c = compute_P(np.sum(V_hat, axis=0), np.sum(N_hat, axis=0))
    return np.sum(binom.logpmf(V_hat, N_hat, P_c))

def min_linkage(pairs_dict, c1, c2):
    """A function to compute the minimum linkage criterion
    
    Parameters
    -----------
    pairs_dict : object
        instance of the PD_Dict class that allows for O(1) access of the distance between an unordered pair of mutations {i,j}
    c1 : list
        a list of 'id' values for the mutations in a particular cluster
    c2 : list
        a list of 'id' values for the mutations in a particular cluster

    Returns
    --------
    float
        the minimum distance between a pair of mutations i \in c1, j \in c2 
    """
    return np.min([pairs_dict.get(v1,v2) for v1 in c1 for v2 in c2]) 

def avg_linkage(pairs_dict, c1, c2):
    """A function to compute the unweighted average linkage criterion
    
    Parameters
    -----------
    pairs_dict : object
        instance of the PD_Dict class that allows for O(1) access of the distance between an unordered pair of mutations {i,j}
    c1 : list
        a list of 'id' values for the mutations in a particular cluster
    c2 : list
        a list of 'id' values for the mutations in a particular cluster

    Returns
    --------
    float
        the average unweighted distance all pairs of mutations (i,j) where i \in c1, j \in c2 
    """
    return np.sum([pd_dict.get(v1,v2) for v1 in c1 for v2 in c2]) / (len(c1)*len(c2))

def compute_P(V, N, epsilon=1e-5):
    """Computes an MLE binomial success probability given V successes of N trials
    
    V : ndarray
        a 2D numpy array where each row is the variant read counts for a particular mutation/cluster
    N : ndarray
        a 2D numpy array where each row is the total read counts for a particular mutation
    epsilon : float, optional
        a minimum cellular prevalence that a mutation can have. This is used to prevent computing log(0))

    Returns
    --------
    ndarray
        a (1D or 2D) numpy array of data-implied cellular prevelance for each mutation
    """
    P = np.divide(V, N)
    P = np.minimum(1 - epsilon, np.maximum(epsilon, P))
    if not isinstance(P, np.ndarray):
        P = np.array([P])
    return P

def load_inputs(input_fn):
    """
    Loads either an omicsdata.npz.archive.Archive zipped archive or a Neutree file

    Parameters
    -----------
    input_fn : str
        A path to either an omicsdata Archive or a Neutree file

    Returns
    --------
    ndarray
        a list of lists where each sub-list is a parents vector
    ndarray
        a list of 2D numpy arrays where each 2D numpy array is the 
        cellular prevalence matrix for a particular parents vector/tree
    """
    structs, Fs = None, None
    exceptions = []
    # try openining as an Archive
    try:
        data = Archive(input_fn)
        structs = data.get(KEYS_ORCH_NPZ.struct_key)
        Fs = data.get(KEYS_ORCH_NPZ.F_key)
    except Exception as e:
        exceptions.append(e)
    
    # try to open as a Neutree if we obtained an exceptionf
    if len(exceptions) > 0:
        try:
            ntree = neutree.load(input_fn)
            structs = ntree.structs
            Fs = ntree.phis
        except Exception as e:
            exceptions.append(e)
            print("Multiple exceptions occurred on loading %s" % input_fn)
            for e in exceptions:
                print(e)
            exit()
    
    return structs, Fs 

def select_clustering(clusterings, llhs, criterion, num_nodes, num_samples, plot=False):
    """
    Parameters
    -----------
    clusterings : list 
        a list of lists, where each sublist of a list of lists describing a particular clustering of mutations
    llhs : list
        a list of floats, where each float is the log-likelihood of the clustering at the same index in clusterings
    criterion : str
        the name of the information criterion to use for model selection 
    num_nodes : int
        the number of nodes in the tree used for clustering
    num_samples : int 
        the number of samples for the dataset
    plot : bool
        a flag to tell the program whether or not to plot the log-likelihood and information criterion scores

    Returns
    ---------
    int 
        the index of the clustering that minimizes the model selection criterion
    """
    n = num_nodes
    m = num_samples
    best_clustering = (0, np.inf) # (index of clone tree, BIC or GIC score)

	# use the model selection criterion to pick a clone tree that minimizes either the BIC or GIC
    scores = []
    for i, llh in enumerate(llhs):
        C = (n-i)*m*model_selection_choices[criterion]((n,m)) - 2*llh
        scores.append(C)
        if C < best_clustering[1]:
            best_clustering = (i, C)

    # plot scores of information criterion compared to cluster log-likelihood
    if plot:
        inv_scores = scores[::-1] # make xaxis go from size 1...n
        inv_llhs = llhs[::-1]
        cluster_size = n - best_clustering[0]
        xaxis = np.arange(1, n+1)
        plt.plot(xaxis, inv_llhs, label="log-likelihood")
        plt.plot(xaxis, inv_scores, label="%s" % criterion)
        plt.scatter(cluster_size, scores[best_clustering[0]], 
                    label="Chosen model (%d nodes)" % (cluster_size), color="green")
        plt.xlabel("Clone tree size")
        plt.legend()
        plt.show()

    return best_clustering[0] 

def prep_cluster_data(variants, cluster):
    """Prepares a cluster's data for computing the binomial likelihood by adjusting their read counts using the 
    variant read probability (similar to a supervariant approximation)
    
    Parameters
    -----------
    variants : dict
        a dictionary where the keys are 'id' values and the values are python dictionaries containing the read
        count data for the mutation signified by a particular 'id'
    cluster : list
        a list of lists where each sublist is a particular cluster containing a set of unique 'id' values denoting
        the mutations in that cluster

    Returns
    --------
    ndarray
        a (1D or 2D) numpy array of the adjusted variant read counts for each mutation
    ndarray
        a (1D or 2D) numpy array of the adjusted total read counts for each mutation
    ndarray
        a (1D or 2D) numpy array of the updated variant read probability for each mutation
    """
    # extract the data for all mutations in each cluster
    cluster_variants = [variants[vid] for vid in cluster]

    # extract the read counts and variant read probabilities for cluster 1
    N = np.array([var[constants.Variants_Keys.TOTAL_READS] for var in cluster_variants])
    V = np.array([var[constants.Variants_Keys.VAR_READS] for var in cluster_variants])
    omega_v = np.array([var[constants.Variants_Keys.OMEGA_V] for var in cluster_variants])
    
    N_hat = 2*N*omega_v 
    V_hat = np.minimum(V, N_hat)
    omega_hat = 0.5*np.ones(len(N_hat))
    
    return V_hat + 1, N_hat + 2, omega_hat

