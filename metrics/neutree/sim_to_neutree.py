###############################################################################################################
# sim_neutree.py
# 
# Source code for enumerating all trees that adhere to the ISA for a simulated dataset
###############################################################################################################
import os, sys, argparse
import pickle
import numpy as np
from omicsdata.tree.adj import adj_to_parents

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "lib", "orchard"))

from constants import KEYS_ORCH_NPZ, KEYS_PARAMS
import neutree

# CONSTANTS
DFS = "dfs"
BFS = "bfs"
STRUCTURE = "structure"
VIDS_GARBAGE = "vids_garbage"

def init_order(F):
	"""Initializes the order in which nodes are sampled for checking placements that incur no error
	
	Parameters
	----------
	F : ndarray
		the data-implied cellular prevalence matrix

	Returns
	-------
	ndarray
		an array that indicates the order in which nodes should be sorted such that they are decreasing in 
		terms of sum of their cellular prevalence's (the same idea as the F sum node ordering)
	"""
	Fsum = np.sum(F, axis=1)
	order = np.argsort(-Fsum)
	return order

def init_tau(F, order):
	"""Initializes an ancestry matrix where there is no violations of the ISA (i.e., tau)
	
	Parameters
	----------
	F : ndarray
		the data-implied cellular prevalence matrix
	order : ndarray
		an array that indicates the order in which nodes should be sorted such that they are decreasing in 
		terms of sum of their cellular prevalence's (the same idea as the F sum node ordering)

	Returns
	-------
	ndarray
		a 2D array that is an ancestry matrix that does not incur any violations of the ISA in the data-implied cellular prevalence matrix
	"""
	K, S = F.shape
	tau = np.eye(K, dtype=np.int8)

	for i in range(K):
		for j in range(i+1, K):
			u = order[i]
			v = order[j]
			if np.all(F[u] >= F[v]):
				tau[u,v] = 1

	return tau

def find_trees(tau, 
			   F, 
			   order, 
			   traversal=DFS, 
			   store_trees=True, 
			   eps=1e-10, 
			   max_count=np.iinfo(np.uint64).max):
	"""
	Enumerates all trees for a dataset that do not have ISA violations 

	Parameters
	----------
	tau : ndarray
		a 2D array that is essentially a adjacency matrix where u adjacent to only if this relationship does not violate the ISA
	F : ndarray
		the data-implied cellular prevalence matrix 
	order : ndarray
		an array that indicates the order in which nodes should be sorted such that they are decreasing in 
		terms of sum of their cellular prevalence's (the same idea as the F sum node ordering)
	traversal : str, optional
		the traversal type to search for trees that do not violate the ISA (either 'dfs' or 'bfs')
	store_trees : bool, optional
		a flag where if True the trees that do not violate the ISA will be returned, otherwise the program will only printout the number of valid trees
	eps : float
		the amount of error that a cellular prevalence can be off by to still adhere to the ISA
	max_count : int
		the number of complete trees that can be returned

	Returns
	-------
	int 
		the number of trees that were found that adhere to the ISA
	list
		a list of lists where each sublist is a parents vector that describes a completed tree
	"""
	assert (traversal == DFS) or (traversal == BFS), "traversal must either be %s or %s" % (DFS, BFS)
	K = len(tau)
	expected_Fsum = np.ones(K)
	expected_Fsum[0] = 0

	tau_ = np.copy(tau)
	np.fill_diagonal(tau_, 0)
	child_Fsum = np.zeros(F.shape)
	partial_trees = [(1, tau_.astype(np.int8), child_Fsum)]
	completed_trees = []
	num_trees = np.uint64(0)

	# until our queue is empty
	while len(partial_trees) > 0:
		if traversal == DFS:
			next_node, adj, child_Fsum = partial_trees.pop()
		else:
			next_node, adj, child_Fsum = partial_trees.pop(0)
		if next_node == K:
			diff = F - child_Fsum
			assert np.all(expected_Fsum == np.sum(adj, axis=0)), "Sum of cellular prevalence's do not add up to 1"
			assert np.all(0 <= child_Fsum + eps) and np.all(child_Fsum <= 1 + eps), "Sum of the descendant cellular prevalence's are not bounded by (0, 1+eps)"
			assert np.all(child_Fsum <= F + eps), "Sum of the descendant cellular prevalence's for one or more nodes violates the ISA"

			num_trees += np.uint8(1)

			# keep track of tree's if we're not just counting the number of trees
			if num_trees == max_count:
				break
			elif store_trees:
				struct = adj_to_parents(adj)
				completed_trees.append(struct)
				continue

		# select the next node to check
		node = order[next_node]

		# find all of the ancestors of 'node'
		ancestors = np.nonzero(adj[:,node])[0]
		
		# for each ancestor, if placing 'node' as a direct descendant does not violates the ISA then 'resolve' that node by making that ancestor
		# the direct ancestor of 'node'
		for ancestor in ancestors:
			child_Fsum_a = np.copy(child_Fsum[ancestor])
			child_Fsum_a += F[node]
			if np.any(child_Fsum_a > F[ancestor] + eps):
				continue
			child_Fsum_ = np.copy(child_Fsum)
			child_Fsum_[ancestor] = child_Fsum_a
			adj_ = np.copy(adj)
			adj_[:,node] = 0
			adj_[ancestor,node] = 1
			partial_trees.append((next_node + 1, adj_, child_Fsum_))

	return (num_trees, completed_trees)


def main():
	parser = argparse.ArgumentParser(
	description="Python script that generates a 'Neutree' zipped archive either by (1) enumerating all trees that do not violate the \
					ISA for a simulated dataset, or (2) simply placing the ground-truth simulated tree in a 'Neutree' zipped archive.",
	formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('--only-count', action='store_true')
	parser.add_argument('--only-write-truth', action='store_true', help='Write only the single true tree structure')
	parser.add_argument('pickle_fn', type=str, 
						help="Name of the simulated data 'pickle' file that contains all of the results for a file simulated by Pearsim")
	parser.add_argument('neutree_fn', type=str,
						help="Name of 'Neutree' file to write the ground-truth trees to")
	args = parser.parse_args()

	# open pickle file
	with open(args.pickle_fn, 'rb') as f:
		pickle_data = pickle.load(f)

	F = pickle_data[KEYS_ORCH_NPZ.F_key]
	order = init_order(F)
	tau = init_tau(F, order)

	# count trees and return
	if args.only_count:
		num_trees, _ = find_trees(tau, F, order, store_trees=False)
		print("Total number of trees: %d" % num_trees)
		return
	elif args.only_write_truth: # write only the ground-truth tree to a 'Neutree' zipped archive
		num_trees = 1
		structs = np.array(pickle_data[STRUCTURE])
	else: # enumerate all trees that do not violate the ISA and write them to a 'Neutree' zipped archive
		num_trees, structs = find_trees(tau, F, order)

	# write trees to 'Neutree' file
	N = len(structs)
	llhs = np.zeros(N)
	probs = np.ones(N) / N
	F_ = np.array([pickle_data[KEYS_ORCH_NPZ.F_key] for _ in range(N)])
	counts = np.ones(N)
	ntree = neutree.Neutree(
	structs = np.array(structs),
	phis = F_,
	counts = counts,
	logscores = llhs,
	clusterings = [pickle_data[KEYS_PARAMS.clusters_key] for idx in range(N)],
	garbage = pickle_data[VIDS_GARBAGE],
	)
	neutree.save(ntree, args.neutree_fn)

if __name__ == '__main__':
  main()