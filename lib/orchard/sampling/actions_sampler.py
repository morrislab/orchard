########################################################################
# action_sampler.py
#
# This contains the source code used to propose nodes to place into 
# a partial tree. This is unfortunately one of the most convoluted
# source files.
########################################################################
import numpy as np 
from scipy.special import betainc, softmax
from itertools import combinations 

# CONSTANTS
U_REPLACES_V = 1
U_BRANCHED_FROM_V = 2
EPSILON = np.exp(-30)
MIN_VAR_READS = 3
MIN_OMEGA=1e-3
CONF_HIGH = "high"
CONF_LOW = "low"
ROOT = 0

def max_log_sum(x, min_value=EPSILON):
    """Returns the max-log-sum of the input x with the comparison for the maximum
    of min_value
    
    Parameters
    -----------
    x : ndarray
        vector of mutation frequencies 
    min_value : float
        the minimum value to have as a mutation frequency

    Returns
    -------
    float
        the sum of the logs of the input x where x >= EPSILON
    """
    return np.sum(np.log(np.maximum(x, min_value)))

def find_valid_samples(V1_var_reads, V1_omega, V2_var_reads, V2_omega):
    """Determine which samples we can compare the read counts for variant 1 (V1) and variant 2 (V2)
    
    Parameters
    ----------
    V1_var_reads : ndarray
        an array of variant allele reads for the first variant (V1)
    V1_omega : ndarray
        an array of variant read probabilities for the first variant (V1)
    V2_var_reads : ndarray
        an array of variant allele reads for the second variant (V2)
    V2_omega : ndarray
        an array of variant read probabilities for the second variant (V2)

    Return
    ------
    ndarray
        a boolean array of which samples we can assume that V1 is ancestral to V2
    ndarray
        a boolean array of samples in which we can compare the two variants
    """
    uninformative_samples = V1_var_reads <= MIN_VAR_READS
    u_parents_v_samples = np.logical_and(V1_var_reads > MIN_VAR_READS, V2_var_reads <= MIN_VAR_READS)
    invalid_omega = np.logical_or(V1_omega < MIN_OMEGA, V2_omega < MIN_OMEGA)
    valid_samples = np.logical_not(np.logical_or(np.logical_or(uninformative_samples, invalid_omega), u_parents_v_samples))
    return u_parents_v_samples, valid_samples

def compute_F_sum(nodes, node_to_index, F_t):
    """Computes a sum over the values in the mutation frequency matrix (F) for one or more nodes
    
    Parameters
    ----------
    nodes : list
        a list of nodes
    node_to_index : dictionary
        a dictionary that maps a node to its index in the mutation frequency matrix (F)
    F_t : ndarray
        the mutation frequency matrix (F) fit only to the nodes in the (partial) tree

    Returns
    -------
    ndarray
        an array that contains the sum of the mutation frequencies in each sample for all of the nodes from the input
    """
    F_t_hat = np.sum([F_t[node_to_index[v]] for v in nodes], axis=0)
    return F_t_hat

def wald_interval(p, n, z=2.576):
    """Computes the Wald confidence interval for the binomial distribution
    
    Parameters
    ----------
    p : ndarray
        an array of data implied mutation frequencies in each sample for a node
    n : ndarray
        an array of total read counts in each sample for a node
    z : float
        target error rate, z=2.576 is the 99% confidence interval

    Returns
    -------
    ndarray
        the lower limit of the confidence interval for the mutation frequency in each sample
    ndarray
        the upper limit of the confidence interval for the mutation frequency in each sample
    """
    q = z*np.sqrt(np.maximum(EPSILON, np.divide(p*(1-p),n,where=n>0)))
    return np.maximum(EPSILON, p-q), np.minimum(1-EPSILON, p+q)

def wilson_ci(p, n, z=2.576):
    """Computes the Wilson confidence interval for the binomial distribution
    
    Parameters
    ----------
    p : ndarray
        an array of data implied mutation frequencies in each sample for a node
    n : ndarray
        an array of total read counts in each sample for a node
    z : float
        target error rate, z=2.576 is the 99% confidence interval

    Returns
    -------
    ndarray
        the lower limit of the confidence interval for the mutation frequency in each sample
    ndarray
        the upper limit of the confidence interval for the mutation frequency in each sample
    """
    q = (z/(1+(z**2/n)))*np.sqrt(((p*(1-p))/n) + (z**2/(4*n**2)))
    p_alt = (1/(1+(z**2/n)))*(p+(z**2/(2*n)))
    return np.maximum(EPSILON, p_alt - q), np.minimum(1-EPSILON, p_alt + q)

def compute_betainc(F_data, u, v, x=np.array([]), confidence_adjustment=None):
    """Computes the incomplete beta function on valid samples and properly sets other values for non-valid samples
    
    Parameters
    ----------
    F_data : dataclass
        an instance of the FData dataclass containing all of the read count and supervariant information
        needed to compute the mutation frequency matrix F
    u : int
        the node value (its index is u-1) for the node being placed into the tree
    v : int
        the node(s) that are used to define the mutation frequency values that u must fall below in each sample
    x : ndarray, optional
        the mutation frequency values that u must fall under
    confidence_adjustment : str, optional
        which confidence interval to use when compute the incomplete beta function

    Returns
    -------
    ndarray
        an array of probabilities for the mutation frequency of u falling below the mutation frequency of v in each sample
    """
    # b_u: variant read counts for u, T_u: total read counts for u, o_u: variant reads probability of u
    b_u, T_u, o_u = F_data.V[u-1], F_data.N[u-1], F_data.omega[u-1] 
    b_v, T_v, o_v = None, None, None

    # compute 'supervariant' approximation if v is more than one value and x isn't defined
    if isinstance(v, np.ndarray):
        b_v, T_v = np.zeros(len(b_u)), np.zeros(len(b_u))
        for c_i in v:
            b_v += F_data.V[c_i-1]
            T_v += 2*F_data.N[c_i-1]*F_data.omega[c_i-1]
        T_v = T_v.astype(int)
        b_v = np.minimum(b_v, T_v).astype(int)
        o_v = np.full(len(b_v), 0.5)
        v_F_hat = np.maximum(EPSILON,np.minimum(1, np.divide(b_v, T_v*o_v, where=T_v*o_v>0)))
    elif v == 0: # create dummy values if v=0
        b_v, o_v = np.full(len(b_u), MIN_VAR_READS+1), np.full(len(b_u), MIN_OMEGA)
        v_F_hat = np.full(len(b_u), 1 - EPSILON)
    else: # grab data for v and compute \hat{\F}_v
        b_v, T_v, o_v = F_data.V[v-1], F_data.N[v-1], F_data.omega[v-1]
        v_F_hat = np.minimum(1 - EPSILON, b_v / (np.maximum(1 ,np.array(T_v*o_v).astype(int))))

    # set x=v_F_hat if we aren't provided with an upper limit for the incomplete beta function
    if len(x) == 0:
        x = v_F_hat
 
    # find samples that will provide valid results
    u_parents_v_samples, valid_samples = find_valid_samples(b_u, o_u, b_v, o_v) 

    # adjust F's using binomial confidence interval  
    if (confidence_adjustment != None) and isinstance(T_v, np.ndarray) and (sum(valid_samples) > 0):
        low, high = wald_interval(x, T_v)
        if confidence_adjustment == CONF_HIGH:
            x[valid_samples] = high[valid_samples]
        elif confidence_adjustment == CONF_LOW:
            x[valid_samples] = low[valid_samples]

    # compute incomplete beta function if we have at least one valid sample
    if sum(valid_samples) > 0:
        results = np.full(len(b_u), 1 - EPSILON)
        results[u_parents_v_samples] = EPSILON
        results[valid_samples] = betainc(
                                b_u[valid_samples]+1, 
                                T_u[valid_samples] - b_u[valid_samples]+1, 
                                np.maximum(0, np.minimum(1, o_u[valid_samples]*x[valid_samples]))
                            )
    else: # otherwise return near-zero probabilities across all samples
        results = np.full(len(b_u), EPSILON)
    return results

class Actions_Sampler:
    """Class for sampling placements of node u in a partial tree
    
    Attributes
    ----------
    __current_node_idx : int
        the index in the __node_order for the next node that will be placed into the partial tree
    __max_placements : int
        the maximum number of placements that can be checked for placing node u in the partial tree
    __num_nodes : int
        the total number of nodes that will be in the complete tree
    __num_samples : int
        the number of samples for each node
    __node_order : list
        the order in which the nodes will be placed into the tree 
    __node_weights : list 
        the sum of the mutation frequencies for each node across all samples (used to compute the F sum node order)
    __ignore_zero_probs : bool
        a flag that prevents the fitting the mutation frequency matrix and scoring partial trees for placements of node u
        that have a probability of zero computed by the node placement function
    __force_monoprimary : bool
        a flag that when true the Action_Sampler does not consider node placements that would result in a poly-primary tree
    """
    def __init__(self, 
                 max_placements,
                 num_nodes, 
                 num_samples, 
                 ignore_zero_probs,
                 force_monoprimary,
                 current_node_idx=1, 
                 node_order=[],
                 node_weights=[]):

        self.__current_node_idx = current_node_idx
        self.__max_placements = max_placements
        self.__num_nodes = num_nodes
        self.__num_samples = num_samples
        self.__node_order = node_order
        self.__node_weights = node_weights
        self.__ignore_zero_probs = ignore_zero_probs
        self.__force_monoprimary = force_monoprimary

    def _compute_F_sum(self, F_data):
        """Computes the F-sum node weights
        
        Parameters
        ----------
        F_data : dataclass
            a dataclass containing the read count data for all supervariants
        """
        self.__node_weights = []
        for u in range(self.__num_nodes):
            # compute the data-implied VAFs, and sum those as the weights
            u_F_hat = np.divide(F_data.V[u], (F_data.N[u] * F_data.omega[u]),where=(F_data.N[u] * F_data.omega[u])>0)
            self.__node_weights.append(np.sum(u_F_hat))

    def _sample_node_order(self, generator, randomize=False):
        """Samples a node order -- either a randomized node order or a F-sum node order
        
        Parameters
        ----------
        generator : object
            a numpy default_rng object used for reproducible random sampling
        randomized : bool, optional
            a flag where if True a randomized node order is generated, and otherwise the F-sum node order is generated
        """
        if randomize:
            node_order = list(range(1,self.__num_nodes+1))
            generator.shuffle(node_order)
        else:
            node_order = np.argsort(self.__node_weights)[::-1] + 1
        self.__node_order = np.insert(node_order, 0, 0)

    def copy(self):
        """Returns a copy of the Action_Sampler object"""
        A = Actions_Sampler(num_nodes=self.__num_nodes, 
                            num_samples=self.__num_samples,
                            current_node_idx=self.__current_node_idx, 
                            max_placements=self.__max_placements,
                            ignore_zero_probs=self.__ignore_zero_probs,
                            force_monoprimary=self.__force_monoprimary,
                            node_order=self.__node_order.copy(),
                            node_weights=self.__node_weights.copy())
        return A

    def propose_actions(self, parents, F_data, F_t, eta_t, branching_factor, generator):
        """Proposes possible placements for the node u
        
        Parameters
        ----------
        parents : ndarray
            a numpy array where each index represents a node, and the value at that index represents
            that nodes direct ancestor (parent) in the tree
        F_data : dataclass
            an instance of the FData dataclass containing all of the read count and supervariant information
            needed to compute the mutation frequency matrix F
        F_t : ndarray
            the mutation frequency matrix (F) fit only to the nodes in the (partial) tree
        eta_t : ndarray
            the subpopulation frequency matrix (eta) fit only to the nodes in the (partial) tree
        branching_factor : int
            the number of possible extensions of the branch to fit the mutation frequency matrix for F and
            score under a binomial likelihood
        generator : object
            a numpy default_rng object used for reproducible random sampling
        
        Returns
        -------
        list
            a list of tuples where each tuple (u,v,c,ACTION) defines the following:
                u : node being placed
                v : parent node of u
                c : nodes that were direct descendants of v, that will now be direct descendents of u (aka q from the paper)
                ACTION: the action being taken either (1) U_REPLACES_V or (2) U_BRANCHED_FROM_V
        """
        actions, probs = [], [] # possible placements for u and the related probabilities
        u = self.__node_order[self.__current_node_idx] # node to place into tree

        nodes_in_tree = self.__node_order[:self.__current_node_idx] 
        argsort_nodes_in_tree = np.argsort(nodes_in_tree)
        node_to_index = {v:argsort_nodes_in_tree.tolist().index(i) for i,v in enumerate(nodes_in_tree)} # map a node to its index in the F/eta matrix

        # iterate through the nodes in the tree and see which placement is best for u
        for i, v_idx in enumerate(argsort_nodes_in_tree):
            v = nodes_in_tree[v_idx] # we do this so we can properly index into the F_t, eta_t matrix

            # (1) u replaces v
            # sum of log probabilities for u being a descendant of v, and u parenting
            # the combinations of children of u
            children = np.where(parents == v)[0] + 1 # needs a plus one

            if len(children) > 0:# and (v in u_replaces_v_considerations):
                n_children = np.minimum(len(children), self.__max_placements) # ensure we don't degrade to exponential runtime
                for n in range(1, n_children+1):
                    for c in combinations(children, n): # powerset of children
                        c = np.array(c) # children considered for u being a parent of
                        rc = np.array(list(set(children)-set(c))) # remaining children u is not a parent of
                        c_F_t = compute_F_sum(c, node_to_index, F_t) # compute sum of F for c
                        rc_F_t = compute_F_sum(rc, node_to_index, F_t) # compute remaining sum of F for children \ c
                        v_betainc = compute_betainc(F_data, u, v, F_t[i]-rc_F_t, confidence_adjustment=CONF_HIGH)
                        c_betainc = compute_betainc(F_data, u, c, c_F_t, confidence_adjustment=CONF_LOW)
                        u_replaces_v = max_log_sum(v_betainc-c_betainc)

                        # u : node being placed 
                        # v : parent node of u 
                        # c : nodes that were direct descendants of v, that will now be direct descendents of u (aka q from the paper)
                        # U_REPLACES_V : the action where u replaces v as the parent of all nodes in c
                        actions.append((u,v,c,U_REPLACES_V))
                        probs.append(u_replaces_v)

            # (2) Special case: v parents u, with u parenting the empty set {}  
            b_betainc = compute_betainc(F_data, u, v, eta_t[i], confidence_adjustment=CONF_HIGH)
            u_branched_from_v = max_log_sum(b_betainc)

            # u : node being placed 
            # v : parent node of u 
            # [] : u is the parent of the empty set (i.e., has no children)
            # U_BRANCHED_FROM_V : the action where u is on its own branch
            if (v == ROOT) and self.__force_monoprimary:
                pass # do nothing if we're forcing a monoprimary tree
            else:
                actions.append((u,v,[],U_BRANCHED_FROM_V))
                probs.append(u_branched_from_v)

        # if we don't have any actions, then we will equally consider placing u as a child of any node v that's already in the tree
        if len(actions) == 0:
            actions = [(u,v,[],U_BRANCHED_FROM_V) for v in nodes_in_tree]
            probs = np.full(len(actions), 1/len(actions))
        else:
            probs = softmax(np.array(probs) - np.max(probs)) # rescale probabilities by largest to prevent precision problems

        # sort actions and add index that allows us to keep track of how the tree was formed
        sorted_indices = np.argsort(probs)[::-1]
        sorted_actions = np.array([(i,u,v,c,model) for i, (u,v,c,model) in enumerate(np.array(actions, dtype=object)[sorted_indices])],dtype=object)
        sorted_probs = np.array(probs)[sorted_indices]

        if self.__ignore_zero_probs: # if true, we'll ignore placements with a probability of zero
            non_zero_probs = sorted_probs > 0
            sorted_actions = sorted_actions[non_zero_probs]
            sorted_probs = sorted_probs[non_zero_probs]

        # take the top f
        chosen_actions = sorted_actions[:np.minimum(branching_factor, len(sorted_probs))]

        self.__current_node_idx += 1 # used to keep track of which nodes we've already placed into the tree
        return chosen_actions