###############################################################################################################
# utils.py
# 
# Contains utility functions for computing loss metrics used to compare cancer phylogeny reconstruction
# methods
###############################################################################################################
import numpy as np 

def compute_membership_matrix(clusters):
	"""Computes the membership matrix for a set of clusters
	
	Parameters
	----------
	clusters : list
		a list of lists where each sublist are the 'id' values for the variants in the same cluster

	Returns
	-------
	list
		the sorted list of 'id' values for the variants, sorted by the numeric value of the 'id'
	ndarray
		a 2D array that is the membership matrix where the rows j=1,...,n are mutations and the
		columns i = 1,...,c are clones where if (j,i) = 1, then mutation j is in clone i
	"""
	vids = sorted([vid for C in clusters for vid in C], key=lambda x: int(x[1:]))
	membership = [np.array([vids.index(vid) for vid in C]) for C in clusters]
	membership_matrix = np.zeros((len(vids), len(clusters)))
	for i, members in enumerate(membership):
		if len(members) == 0:
			continue
		membership_matrix[members,i] = 1
	return (vids, membership_matrix)
