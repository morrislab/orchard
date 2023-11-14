########################################################################
# pa_clustering.py
#
# The source code for the 'phylogeny-aware' clustering algorithm
########################################################################

import numpy as np
from copy import deepcopy 
from utils import cluster_llh, PD_Dict, NULL_NODE, INF
from itertools import chain
from tqdm import tqdm

def _compute_Q(parents, Q, pairs_dict, clusters, linkage_criterion):
    """Computes the distance matrix Q for all adjacent nodes
    
    Parameters
    -----------
    parents : ndarray
        a numpy array where each index represents a node, and the value at that index represents
        that nodes direct ancestor (parent) in the tree
    Q : ndarray
        a 2D numpy array where Q[i,j] represents the distance between the cellular prevalence's of node i and node j
    pairs_dict : object
        instance of the PD_Dict class that allows for O(1) access of the distance between an unordered pair of mutations {i,j}
    clusters : ndarray
        a vector of vectors, where each sub-vector is the 'id' values of the mutations that are in the same cluster
    linkage_criterion : function
        the linkage criterion function

    Returns 
    --------
    ndarray
        the matrix Q that contains the updated distances between all pairs of nodes left in the tree
    """
    # go through each node that's still in the tree
    for node1 in np.where(parents>=0)[0]:
 
        # compute error for joining node1 with each of its
        # children (node2) if this value needs to be recomputed
        possible_nj = np.where(parents==(node1+1))[0] + 1
        for node2 in np.array(possible_nj, dtype=np.int64) - 1:
            if (node1 == node2) or (Q[node1, node2] != INF):
                continue
            else:
                # compute distance between node1 and node2 using a linkage criterion
                dist = linkage_criterion(pairs_dict, clusters[node1], clusters[node2])
                Q[node1, node2] = dist
    return Q

def _join(parents, Q, clusters):
    """Joins two nodes together that minimize the cluster linkage criterion
    
    Parameters
    -----------
    parents : ndarray
        a numpy array where each index represents a node, and the value at that index represents
        that nodes direct ancestor (parent) in the tree
    Q : ndarray
        a 2D numpy array where Q[i,j] represents the distance between the cellular prevalence's of node i and node j
    clusters : ndarray
        a vector of vectors, where each sub-vector is the 'id' values of the mutations/clusters that are in the same cluster

    Returns
    --------
    ndarray
        an updated parents vector containing one fewer nodes in its equivalent tree representation
    ndarray
        an updated Q matrix that is modified such that only the distance for the updated node needs to be computed
    ndarray
        an updated clusters vector, where all of the mutations in node2 are now in node1's cluster
    """
    # return existing clustering if there are no further ways to condense
    if Q.min() == INF: 
        return parents, Q, clusters
    else:
        node1, node2 = np.unravel_index(np.argmin(Q), Q.shape)

    # make the children of node2 be parented by node2's parent
    parents[np.where(parents==node2+1)[0]] = parents[node2]
    parents[node2] = NULL_NODE

    # reset Q matrix so we recompute the distance between our new node and its neighbors
    Q[:,node1] = INF
    Q[node1] = INF
    Q[:,node2] = INF
    Q[node2] = INF

    # update clusters
    clusters[node1] = clusters[node1].union(clusters[node2])
    clusters[node2].clear()

    return parents, Q, clusters

def cluster(supervars, parents, pd_matrix, cluster_ids, original_clusters, F, linkage_func):
    """Main function for performing agglomerative clustering
    
    Parameters
    -----------
    supervars : dict
        a dictionary where the keys are 'id' values and the values are python dictionaries containing the read
        count data for the mutation/cluster signified by a particular 'id'
    parents : ndarray
        a numpy array where each index represents a node, and the value at that index represents
        that nodes direct ancestor (parent) in the tree. This initial parents vector represents a complete tree on n nodes
    pd_matrix : ndarray
        a 2D numpy array containing the pairwise distances between each pair of mutations/clusters i and j
    cluster_id : list
        a list of 'id' values for all clusters/mutations in the tree 
    original_clusters : list 
        a list of lists, where each sublist contains the mutations for each cluster in the original data
    linkage_func : function
        a linkage function for computing a particular type of linkage criterion (min linkage, average linkage, etc.)

    Returns
    --------
    list 
        a list of log-likelihoods under a binomial likehood model for each of the N clone trees
    list
        a list of lists where each sublist is the mutation clusters for a particular clone tree
    list
        a list of lists where each sublist is the parents vector for a particular clone tree
    """
    # number of iterations
    N = len(parents)

    # compute pairwise dictionary containing all coclustering probabilities
    pd_dict = PD_Dict()
    for i, vid1 in enumerate(cluster_ids[:-1]):
        for j, vid2 in enumerate(cluster_ids[i+1:]):
            pd_dict.put(vid1, vid2, pd_matrix[i,i+j+1])
        
    Q = np.full((N,N), INF) # matrix of dissimilarities between adjacent nodes
    
    # initialize variables
    clusters = [set({vid}) for vid in cluster_ids]
    llhs_ = [np.sum(cluster_llh(supervars, c, F) for c in clusters)]
    clusters_ = [[list(c) for c in clusters]]
    Q_list = []
    parents_ = [deepcopy(parents)]

    with tqdm(total=N-1) as pbar:
        pbar.set_description("Running Phylogeny-aware clustering")
        for _ in range(N-1):
            # agglomerative clustering
            Q = _compute_Q(parents, Q, pd_dict, clusters, linkage_criterion=linkage_func)
            parents, Q, clusters = _join(parents, Q, clusters)

            # compute metrics for clustering and save data
            llh = np.sum(cluster_llh(supervars, c, F) for c in clusters if len(c) > 0)
            Q_list.append(deepcopy(Q))
            parents_.append(deepcopy(parents))
            llhs_.append(llh)
            clusters_.append([list(chain.from_iterable([original_clusters[int(sid[1:])-1] for sid in c])) for c in clusters])
            pbar.update()
    return llhs_, clusters_, parents_

