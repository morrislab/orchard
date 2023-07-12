###############################################################################################################
# neutree.py
# 
# Contains the source code for reading/writing 'Neutree' namedtuples 
###############################################################################################################
import pickle
from collections import namedtuple
import numpy as np

# for compatibility with Pairtree, we reuse the Neutree named tuple
Neutree = namedtuple('Neutree', ('structs', 'phis', 'counts', 'logscores', 'clusterings', 'garbage'))

def save(ntree, neutree_fn):
	"""Saves the data for a bulk DNA cancer phylongeny reconstruction in a generalized format
	that's simply a zipped archive containing a namedtuple
	
	Parameters
	----------
	ntree : namedtuple
		the name tuple that will be written to a zipped archive
	neutree_fn : str
		the file name that the ntree namedtuple will be written to

	Returns
	-------
	None
	"""
	N = len(ntree.structs)
	for K in ('structs', 'phis', 'counts', 'logscores', 'clusterings'):
		assert len(getattr(ntree, K)) == N, '%s has length %s instead of %s' % (K, len(getattr(ntree, K)), N)

	# we always expect data in the Neutree archive to be ndarray's
	arr_vals = {K: np.array(getattr(ntree, K)) for K in ('counts', 'logscores')}
	ntree = ntree._replace(**arr_vals)

	with open(neutree_fn, 'wb') as F:
		pickle.dump(ntree, F)

def load(neutree_fn):
	"""Loads the Neutree namedtuple from a zipped archive
	
	Parameters
	----------
	neutree_fn : str
		the file name that the ntree namedtuple will be loaded from

	Returns
	-------
	pickle
		a pickle file loaded into memory
	"""
	with open(neutree_fn, 'rb') as F:
		return pickle.load(F)
