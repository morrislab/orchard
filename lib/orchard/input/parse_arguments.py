###########################################################################
# parse_arguments.py
#
# This contains the source code parsing the arguments provided to
# the orchard python binary (i.e., /path/to/orchard/bin/orchard)
###########################################################################
import argparse, os, sys
from multiprocessing import cpu_count
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "models"))

from model_types import sbs, MODEL_TYPES
from constants import F_SUM_NODE_ORDER, RANDOM_NODE_ORDER, DIVERSE_NODE_ORDER, node_order_options


def parse_arguments():
    """Parses command line arguments for running orchard"""
    parser = argparse.ArgumentParser(
        description="Orchard: building large mutation trees using stochastic combinatorial search",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("ssm_fn", type=str, 
                        help="Path to a simple somatic mutation file (.ssm).")

    parser.add_argument("params_fn", type=str, 
                        help="Path to a parameter file (.params.json).")

    parser.add_argument("results_fn", help="Zipped archive (.npz) file name to write results to.")

    parser.add_argument("-d", "--debug", default=False, action="store_true",
                         help="Flag for debug printouts.")      

    parser.add_argument("-e", "--expansion-factor", default=1, type=int,
                         help="Number of partial solutions to expand at each iteration of the Orchard algorithm. \
                               If the expansion_factor=1, this amounts to a best-first-search. If the \
                                expansion_factor=beam_width, this amounts to a breadth-first-search.")

    parser.add_argument("-f", "--branching-factor", type=int, default=20,
                         help="Branching factor controls the number of possible parents to consider when choosing an unparented mutation to add to the tree. \
                               The time complexity of Orchard scales linearly with the branching factor.")

    parser.add_argument("-i", "--num-instances", type=int, default=None,
                         help="Number of parallel instances of Orchard to run.")

    parser.add_argument("-k", "--beam-width", type=int, default=1,
                         help="Beam width is used to limit the number of trees sampled for each chain. \
                              The time complexity of Orchard scales linearly with the beam width.")     

    parser.add_argument("-m", "--model", type=str, default=sbs, choices=list(MODEL_TYPES.keys()),
                         help="The type of model to run")

    parser.add_argument("-n", "--num-cpu-cores", default=cpu_count(), type=int,
                         help="Number of CPUs to use for processing candidate trees.")

    parser.add_argument("-p", "--force-monoprimary", default=False, action="store_true",
                         help="Flag that when provided will force Orchard to only search for trees that are monoprimary.")

    parser.add_argument("-r", "--rescale-depth", default=False, action="store_true",
                         help="Rescales the read depth for each sample using the average read depth across the sample. This should be done if the \
                               data has variable coverage in one or more samples. This should not be done if you believe the low coverage for some samples \
                               is due to errors in sequencing.")

    parser.add_argument("-s", "--seed", type=int, default=None,
                         help="Seed for duplicating results.")

    parser.add_argument("-w", "--node-order", type=str, default=F_SUM_NODE_ORDER, choices=node_order_options,
                         help="Option for what order the nodes will be added to the tree in each parallel instance of Orchard. The '%s' will initialize the \
                         node order for each parallel instance of Orchard by the sum of the data-implied cellular prevalence for each node in descending order (F sum node order). \
                         The '%s' will initialize the node order for each parallel instance of Orchard randomly.  \
                         The '%s' will initialize one of the parallel instances to have the F sum node order, and all of the remaining instances to have a randomized node order" \
                          %(F_SUM_NODE_ORDER, RANDOM_NODE_ORDER, DIVERSE_NODE_ORDER))

    parser.add_argument("-x", "--max-placements", type=int, default=10,
                         help="If we are placing some node u as a direct descendant of a node v that is already in the tree \
                               and v has a set of X direct descendants, then max_placements is used as an upper bound on the number \
                               of combinations of children that u could parent in X such that we will try a total of \
                               total_placements_under_v = \sum_{i=1}^{minimum(max_placements, |X|)}{\binom{|X|}{i}}). \
                               This means that max_placements fixes an exponential upper bound on the number of placements for the \
                                u-parents-v placement into the tree.")

    parser.add_argument("-z", "--ignore-zero-probs", default=False, action="store_true",
                         help="Flag that when provided will have Orchard ignore node placements that have a softmax probability of zero.")
                         
    return parser.parse_args()

def process_args(args):
    """Processes command line arguments
    
    Parameters
    ----------
    args : object
        argparse object containing all of the parsed command line arguments

    Returns
    -------
    int
        the number of workers to use in the ProcessPool for parallelizing Orchard
    int
        the number of individual instances of the Orchard algorithm to run in parallel
    int
        the seed used to initialize the numpy.random.default_rng() object with
    """
    # determine how many CPU cores to use
    poolsize = min(max(1, args.num_cpu_cores), cpu_count())

    # determine number of chains to run
    if args.num_instances is None:
        args.num_instances = max(1, poolsize) 

    # select random seed if it's not defined
    if args.seed is None:
        args.seed = np.random.randint(2**32)
    
    return poolsize, args.num_instances, args.seed
 
