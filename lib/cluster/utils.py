#############################################################################
# utils.py
#
# Utility functions for the 'phylogeny-aware' clustering algorithm
#############################################################################
import os, sys
import numpy as np 
from scipy.stats import binom
import matplotlib.pyplot as plt
from omicsdata.ssm import constants, supervariants
from omicsdata.tree import adj, parents
from omicsdata.ssm.columns import PARAMS_Columns
from omicsdata.tree import neutree
from omicsdata.tree.columns import NEUTREE_Columns
from omicsdata.npz.archive import Archive

# CONSTANTS
LLH_KEY = "llhs"
CLUSTERS_KEY = "clusters"
PARENTS_KEY = "parents"
SV_CLUSTERS_KEY = "sv_clusters"
F_KEY = "F"
MIN_LINKAGE="min_linkage"
AVG_LINKAGE="avg_linkage"
WARD_LINKAGE="ward_linkage"

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

def cluster_llh(variants, cluster, F):
    """Computes the binomial likelihood a single cluster
    
    Parameters
    ------------
    variants : dict
        a dictionary where the keys are 'id' values and the values are python dictionaries containing the read
        count data for the mutation signified by a particular 'id'
    cluster : list
        a list of lists where each sublist is a particular cluster containing a set of unique 'id' values denoting
        the mutations in that cluster
    F : ndarray
        a matrix of mutation frequencies for the mutations in the cluster

    Returns
    ---------
    float
        the log-likelihood of all of the clusters under a binomial likelihood model
    """
    V_hat, N_hat, omega_v, F_cluster = prep_cluster_data(variants, cluster, F)
    P_c = np.sum(omega_v*F_cluster*N_hat, axis=0)/np.sum(N_hat,axis=0)
    return binom.logpmf(V_hat, N_hat, P_c).sum()


def dfs_find_lineages(parents, ROOT=0, NO_CLUSTER=-1):
    """Computes a list of list, where each sublist corresponds to a lineage found in the tree
    
    Parameters
    -----------
    parents : ndarray
        a numpy array where each index represents a node, and the value at that index represents
        that nodes direct ancestor (parent) in the tree
    ROOT : int, optional
        The integer value that represents the root node. Default = 0.
    NO_CLUSTER : int, optional
        The integer value that represents a node is part of another cluster. Default = -1.

    Returns
    ---------
    list
        a list where each sublist contains the indices for all nodes on the same branch in the tree
    """
    parents_copy = np.copy(parents)
    parents_ = [ROOT]
    queue = [ROOT]
    lineages = []
    new_lineage = []
    while len(queue) > 0:
        u = queue.pop()
        children = np.where(parents_copy == u)[0] + 1
        parents_copy[children-1] = NO_CLUSTER
        queue += children.tolist()
        if u != ROOT:
            new_lineage.append(u)
        if (len(children) > 1 or len(children) == 0) and len(new_lineage) > 0:
            lineages.append(new_lineage.copy())
            for _ in range(len(children)):
                parents_.append(len(lineages))
            new_lineage = [] 
    if len(new_lineage) > 0:
        lineages.append(new_lineage)
    return lineages, parents_

def min_linkage(pairs_dict, variants, c1, c2):
    """A function to compute the minimum linkage criterion
    
    Parameters
    -----------
    pairs_dict : object
        instance of the PD_Dict class that allows for O(1) access of the distance between an unordered pair of mutations {i,j}
    variants : dict
        a dictionary where the keys are 'id' values and the values are python dictionaries containing the read
        count data for the mutation signified by a particular 'id'
    c1 : list
        a list of 'id' values for the mutations in a particular cluster
    c2 : list
        a list of 'id' values for the mutations in a particular cluster

    Returns
    --------
    float
        the minimum distance between a pair of mutations i in c1, j in c2 
    """
    return np.min([pairs_dict.get(v1,v2) for v1 in c1 for v2 in c2]) 

def avg_linkage(pairs_dict, variants, c1, c2):
    """A function to compute the unweighted average linkage criterion
    
    Parameters
    -----------
    pairs_dict : object
        instance of the PD_Dict class that allows for O(1) access of the distance between an unordered pair of mutations {i,j}
    variants : dict
        a dictionary where the keys are 'id' values and the values are python dictionaries containing the read
        count data for the mutation signified by a particular 'id'
    c1 : list
        a list of 'id' values for the mutations in a particular cluster
    c2 : list
        a list of 'id' values for the mutations in a particular cluster

    Returns
    --------
    float
        the average unweighted distance all pairs of mutations (i,j) where i in c1, j in c2 
    """
    return np.sum([pairs_dict.get(v1,v2) for v1 in c1 for v2 in c2]) / (len(c1)*len(c2))

def ward_linkage(F, variants, c1, c2):
    """A function to compute the ward linkage criterion. 
    
    Parameters
    -----------
    F : ndarray
        a matrix of mutation frequencies for all mutations
    variants : dict
        a dictionary where the keys are 'id' values and the values are python dictionaries containing the read
        count data for the mutation signified by a particular 'id'
    c1 : list
        a list of 'id' values for the mutations in a particular cluster
    c2 : list
        a list of 'id' values for the mutations in a particular cluster

    Returns
    --------
    float
        the distance for the Ward linkage criterion
    """
    V_hat1, N_hat1, omega_v1, F_cluster1 = prep_cluster_data(variants, c1, F)
    F_c1 = np.mean(F_cluster1,axis=0)

    V_hat2, N_hat2, omega_v2, F_cluster2 = prep_cluster_data(variants, c2, F)
    F_c2 = np.mean(F_cluster2,axis=0)

    n1 = F_cluster1.shape[0]
    n2 = F_cluster2.shape[0]
    return ((n1*n2)/(n1+n2)) * np.power(F_c1 - F_c2, 2).sum() # ward linkage ((n1*n2)/(n1+n2))*sum((m1 - m2)**2)

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
    ndarray
        a list of log-likelihoods for each of the parents vectors
    """
    structs, Fs = None, None
    exceptions = []
    # try openining as an Archive
    try:
        data = Archive(input_fn)
        structs = data.get("struct")
        Fs = data.get("phi")
        llhs = data.get("llh")

    except Exception as e:
        exceptions.append(e)
    
    # try to open as a Neutree if we obtained an exception
    if len(exceptions) > 0:
        try:
            ntree = neutree.load(input_fn)
            structs = ntree.structs
            Fs = ntree.phis
            llhs = ntree.logscores
        except Exception as e:
            exceptions.append(e)
            print("Multiple exceptions occurred on loading %s" % input_fn)
            for e in exceptions:
                print(e)
            exit()
    
    return structs, Fs, llhs

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
    clustering_sizes = []

	# use the model selection criterion to pick a clone tree that minimizes either the BIC or GIC
    scores = []
    for i, (llh, cl) in enumerate(zip(llhs, clusterings)):
        C = 2*(len(cl))*m*model_selection_choices[criterion]((n,m)) - 2*llh
        scores.append(C)
        if C < best_clustering[1]:
            best_clustering = (i, C)

    # plot scores of information criterion compared to cluster log-likelihood
    if plot:
        inv_llhs = llhs[::-1]
        xaxis = [len(cl) for cl in clusterings[::-1]]
        cluster_size = len(clusterings[best_clustering[0]])
        plt.plot(xaxis, inv_llhs, label="log-likelihood")
        plt.plot(xaxis, scores[::-1], label="%s" % criterion)
        plt.scatter(cluster_size, scores[best_clustering[0]], 
                    label="Chosen model (%d nodes)" % cluster_size, color="green")
        plt.xlabel("Clone tree size")
        plt.legend()
        plt.show()

    return best_clustering[0] 

def prep_cluster_data(variants, cluster, F, eps=1e-4):
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
    F : ndarray
        a matrix of mutation frequencies for all mutations
    eps : float
        the minimum value (and 1-eps is maximum) for entries in the mutation frequency matrix

    Returns
    --------
    ndarray
        a (1D or 2D) numpy array of the adjusted variant read counts for each mutation
    ndarray
        a (1D or 2D) numpy array of the adjusted total read counts for each mutation
    ndarray
        a (1D or 2D) numpy array of the updated variant read probability for each mutation
    ndarray
        a (1D or 2D) numpy array of the mutation frequencies for each mutation
    """
    # extract the data for all mutations in each cluster
    cluster_variants = [variants[vid] for vid in cluster]

    # extract the read counts and variant read probabilities for cluster 1
    N = np.array([var[constants.Variants_Keys.TOTAL_READS] for var in cluster_variants])
    V = np.array([var[constants.Variants_Keys.VAR_READS] for var in cluster_variants])
    omega_v = np.array([var[constants.Variants_Keys.OMEGA_V] for var in cluster_variants])
    
    N_hat = 2*N*omega_v 
    V_hat = np.minimum(V, N_hat)
    omega_hat = 0.5*np.ones(N_hat.shape)
    F_cluster = np.minimum(np.maximum(np.array([F[int(vid[1:])-1] for vid in cluster]), eps), 1-eps)
    
    return V_hat + 1, N_hat + 2, omega_hat, F_cluster