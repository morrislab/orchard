########################################################################
# results.py
#
# Contains the source code for exporting the results for a run by
# Orchard.
########################################################################
import sys, os
import numpy as np
from scipy.special import softmax

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "projection"))

from constants import KEYS_PARAMS, KEYS_ORCH_NPZ
from projection import calc_llh, fit_F
from omicsdata.npz.archive import Archive
from omicsdata.tree.newick import adj_to_newick
from omicsdata.tree.parents import parents_to_adj

def results_to_npz(results_data, F_data):
    """
    Prepares the output of Orchard and uses the Archive class from omicsdata to write numpy zipped file (npz)

    Parameters
    ----------
    results_data : dataclass
        an instance of the ResultsData dataclass containing all of the information needed to compute the 
        final results
    F_data : dataclass
        an instance of the FData dataclass containing all of the read count and supervariant information
        needed to compute the cellular prevalence matrix F

    Returns
    -------
    None
    """   

    # find unique trees via matching parents
    unique_trees, unique_parents = [], {}
    for t in results_data.best_trees:
        tree_hash = hash(t)
        if tree_hash in unique_parents:
            unique_parents[tree_hash] += 1
            continue
        else:
            unique_trees.append(t)
            unique_parents[tree_hash] = 1

    # summarize the space of the trees using the tree's found during search
    unique_trees, unique_parents = [], {}
    for t in results_data.best_trees:
        tree_hash = hash(t)
        if tree_hash in unique_parents:
            unique_parents[tree_hash] += 1
            continue
        else:
            unique_trees.append(t)
            unique_parents[tree_hash] = 1

    # mimic the .npz outputs from pairtree
    npz_data = {
        KEYS_ORCH_NPZ.action_info: [],
        KEYS_ORCH_NPZ.branches_explored: results_data.branches_explored,
        KEYS_ORCH_NPZ.branches_cut: results_data.branches_cut,
        KEYS_ORCH_NPZ.clusters_key: results_data.params[KEYS_PARAMS.clusters_key],
        KEYS_ORCH_NPZ.count_key: [],
        KEYS_ORCH_NPZ.llh_key: [],
        KEYS_ORCH_NPZ.F_key: [],
        KEYS_ORCH_NPZ.eta_key: [],
        KEYS_ORCH_NPZ.prob_key: [],
        KEYS_ORCH_NPZ.sampnames: results_data.params[KEYS_PARAMS.samples_key],
        KEYS_ORCH_NPZ.seed: results_data.seed,
        KEYS_ORCH_NPZ.struct_key: [],
        KEYS_ORCH_NPZ.newick_key: [],

        # just for compatibility with pairtree
        KEYS_ORCH_NPZ.garbage: [], 
        KEYS_ORCH_NPZ.clustrel_posterior_vids: [],
        KEYS_ORCH_NPZ.clustrel_posterior_rels: []
    }

    # fit F's to trees
    for t in unique_trees:
        F_llh, _, _ = calc_llh(t.F(), F_data.V, F_data.N, F_data.omega)
        newick_fmt = adj_to_newick(parents_to_adj(t.parents()))
        
        npz_data[KEYS_ORCH_NPZ.newick_key].append(newick_fmt)
        npz_data[KEYS_ORCH_NPZ.llh_key].append(np.sum(F_llh))
        npz_data[KEYS_ORCH_NPZ.F_key].append(t.F())
        npz_data[KEYS_ORCH_NPZ.eta_key].append(t.eta())
        npz_data[KEYS_ORCH_NPZ.struct_key].append(t.parents())
        npz_data[KEYS_ORCH_NPZ.count_key].append(unique_parents[hash(t)])
        npz_data[KEYS_ORCH_NPZ.action_info].append(t.action_info())

    # compute posterior probability using llh data and sort by posterior
    if len(npz_data[KEYS_ORCH_NPZ.llh_key]) > 0:
        # sort by log likelihood and parsimony
        nllhs = -np.array(npz_data[KEYS_ORCH_NPZ.llh_key])
        order = nllhs.argsort()
        parsimony_scores = [(np.unique(p, return_counts=True)[1] - 1).sum() for p in npz_data[KEYS_ORCH_NPZ.struct_key]]
        npz_data[KEYS_ORCH_NPZ.eta_key] = np.array(npz_data[KEYS_ORCH_NPZ.eta_key])[order]
        npz_data[KEYS_ORCH_NPZ.F_key] = np.array(npz_data[KEYS_ORCH_NPZ.F_key])[order]
        npz_data[KEYS_ORCH_NPZ.prob_key] = np.array(softmax(npz_data[KEYS_ORCH_NPZ.llh_key]))[order]
        npz_data[KEYS_ORCH_NPZ.struct_key] = np.array(npz_data[KEYS_ORCH_NPZ.struct_key])[order] # needs to be a numpy array
        npz_data[KEYS_ORCH_NPZ.count_key] = np.array(npz_data[KEYS_ORCH_NPZ.count_key])[order]
        npz_data[KEYS_ORCH_NPZ.llh_key] = np.array(npz_data[KEYS_ORCH_NPZ.llh_key])[order]
        npz_data[KEYS_ORCH_NPZ.action_info] = np.array(npz_data[KEYS_ORCH_NPZ.action_info],dtype=object)[order]
        npz_data[KEYS_ORCH_NPZ.newick_key] = np.array(npz_data[KEYS_ORCH_NPZ.newick_key])[order]

 
    # save npz data via the Archive class
    results = Archive(results_data.results_fn)

    for k, v in npz_data.items():
        results.add(k, v)

    results.save()