#############################################################################
# generate_clonetree
#
# A Python program for choosing a clone tree from the output of the
# phylogeny-aware clustering algorithm.
############################################################################# 
import os, sys
import numpy as np 
import argparse, json
from omicsdata.ssm import parse, supervariants
from omicsdata.npz.archive import Archive
from scipy.special import softmax

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', "orchard")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', "orchard", "projection")))

from utils import LLH_KEY, CLUSTERS_KEY, PARENTS_KEY, GIC, model_selection_choices, select_clustering
from projection import fit_F
from constants import KEYS_PARAMS, KEYS_ORCH_NPZ

def main():
    parser = argparse.ArgumentParser(
        description="Generates a zipped archive that contains a selected tree from the zipped archive output \
                        by the phylogeny-aware clustering program",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("ssm_fn", type=str, help="Simple somatic mutation file (.ssm)")

    parser.add_argument("params_fn", type=str, help="Parameters file (.params.json)")

    parser.add_argument("cluster_npz", type=str, help="Phylogeny-aware clustering zipped archive")

    parser.add_argument("cluster_params", type=str, help="New params file to write new clustering data to")

    parser.add_argument("results_fn", type=str, help="Zipped archive for a particular tree")

    parser.add_argument("-m", "--model-selection", type=str, default=GIC, choices=list(model_selection_choices.keys()),
                        help="Allows for user to pick the model selection criterion. This is only used if --num-clones is less than \
                            or equal to zero.")

    parser.add_argument("-n", "--num-clones", type=int, default=0, 
                        help="Tree size to choose from the phylgeny-aware clustering zipped archive. \
                                By default (-n=0), the clone tree size is selected using the information criterion")

    parser.add_argument("-p", "--plot", action="store_true", default=False, 
                        help="Flag to plot log-likelihood and GIC/BIC for each size clustering.")

    args = parser.parse_args()

    # extract ssm/params daata
    variants = parse.load_ssm(args.ssm_fn)
    params = parse.load_params(args.params_fn)
    n = len(list(variants.keys()))
    m = len(params[KEYS_PARAMS.samples_key])

    # extract data from cluster archive
    results = Archive(args.cluster_npz)
    clusterings = results.get(CLUSTERS_KEY)
    all_parents = results.get(PARENTS_KEY)
    llhs = results.get(LLH_KEY)

    # get index of clone tree to use
    if args.num_clones > 0:
        index = len(all_parents) - args.num_clones 
    else:
        index = select_clustering(clusterings, llhs, args.model_selection, n, m, args.plot)
    parents = np.array(all_parents[index])
    clusters = clusterings[index]

    # overwrite clusters in inputted file, and save as new params file
    params[KEYS_PARAMS.clusters_key] = clusters
    with open(args.cluster_params, 'w') as F:
        json.dump(params, F)

    # extract read count data
    supervars = supervariants.clusters_to_supervars(clusters, variants)
    V, R, omega_v = supervariants.supervars_to_binom_params(supervars)

    # mimic the .npz outputs from pairtree
    npz_data = {
        KEYS_ORCH_NPZ.action_info: [],
        KEYS_PARAMS.clusters_key: clusters,
        KEYS_ORCH_NPZ.count_key: [],
        KEYS_ORCH_NPZ.llh_key: [],
        KEYS_ORCH_NPZ.F_key: [],
        KEYS_ORCH_NPZ.eta_key: [],
        KEYS_ORCH_NPZ.prob_key: [],
        KEYS_ORCH_NPZ.sampnames: params[KEYS_PARAMS.samples_key],
        KEYS_ORCH_NPZ.struct_key: [],

        # just for compatibility with pairtree
        KEYS_ORCH_NPZ.garbage: [], 
        KEYS_ORCH_NPZ.clustrel_posterior_vids: [],
        KEYS_ORCH_NPZ.clustrel_posterior_rels: []
    }

    # fit the cellular prevalence matrix for the selected clone tree
    phi, eta, phi_llh = fit_F(parents, 
                                    V, 
                                    R,
                                    omega_v)

    npz_data[KEYS_ORCH_NPZ.llh_key].append(np.sum(phi_llh))
    npz_data[KEYS_ORCH_NPZ.F_key].append(phi)
    npz_data[KEYS_ORCH_NPZ.eta_key].append(eta)
    npz_data[KEYS_ORCH_NPZ.struct_key].append(parents)
    npz_data[KEYS_ORCH_NPZ.count_key].append(1)

    # compute posterior probability using llh data and sort by posterior
    if len(npz_data[KEYS_ORCH_NPZ.llh_key]) > 0:
        order = np.array(npz_data[KEYS_ORCH_NPZ.llh_key]).argsort()[::-1]
        npz_data[KEYS_ORCH_NPZ.eta_key] = np.array(npz_data[KEYS_ORCH_NPZ.eta_key])[order]
        npz_data[KEYS_ORCH_NPZ.F_key] = np.array(npz_data[KEYS_ORCH_NPZ.F_key])[order]
        npz_data[KEYS_ORCH_NPZ.prob_key] = np.array(softmax(npz_data[KEYS_ORCH_NPZ.llh_key]))[order]
        npz_data[KEYS_ORCH_NPZ.struct_key] = np.array(npz_data[KEYS_ORCH_NPZ.struct_key])[order] # needs to be a numpy array
        npz_data[KEYS_ORCH_NPZ.count_key] = np.array(npz_data[KEYS_ORCH_NPZ.count_key])[order]
        npz_data[KEYS_ORCH_NPZ.llh_key] = np.array(npz_data[KEYS_ORCH_NPZ.llh_key])[order]

    # save npz data via results archive
    tree_results = Archive(args.results_fn)

    for k, v in npz_data.items():
        tree_results.add(k, v)

    tree_results.save()

if __name__ == '__main__':
    main()