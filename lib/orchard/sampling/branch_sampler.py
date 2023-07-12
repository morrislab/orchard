########################################################################
# branch_sampler.py
#
# Contains the source code for sampling unique extensions of a 
# Branch object.
########################################################################
import sys, os
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "projection"))

from projection import fit_F
from constants import NO_RELATIONSHIP
from gumbel import gumbel_with_maximum
from omicsdata.tree.parents import compute_partial_parents
from actions_sampler import U_REPLACES_V

def F_llh(partial_parents, F_data, bool_mask):
    """Computes the log-likelihood of the cellular prevalence matrix F
    
    Parameters
    ----------
    partial_parents : ndarray
        an array of size l such that 1 <= l <= n containing only the parents information
        for the nodes currently in the (partial) tree
    F_data : dataclass
        an instance of the FData dataclass containing all of the read count and supervariant information
        needed to compute the cellular prevalence matrix F
    bool_mask : ndarray
        a boolean mask that is used to extract the read count data only for the nodes current in the
        (partial) tree

    Returns
    -------
    ndarray
        the cellular prevalence matrix F fit to the (partial) tree
    ndarray
        the subpopulation frequency matrix eta fit to the (partial) tree
    float
        the log-likelihood of the cellular prevalence matrix fit to the (partial) tree
    """
    F, eta, F_llh = fit_F(partial_parents, 
                          F_data.V[bool_mask], 
                          F_data.N[bool_mask] - F_data.V[bool_mask],
                          F_data.omega[bool_mask])
    return F, eta, F_llh

def propose_branches(branch, F_data, branching_factor, generator):
    """Proposes a set of actions based on the current state of the tree 
    and the read count data
    
    Parameters
    ----------
    branch : object
        an instance of the Branch class for which extensions will be proposed
    F_data : dataclass
        an instance of the FData dataclass containing all of the read count and supervariant information
        needed to compute the cellular prevalence matrix F
    branching_factor : int
        the number of possible extensions of the branch to fit the cellular prevalence matrix for F and
        score under a binomial likelihood
    generator : object
        a numpy default_rng object used for reproducible random sampling

    Returns
    -------
    list 
        a list of Branch objects that are unique extensions of the branch passed in
    """
    chosen_actions = branch.propose_actions(F_data, branching_factor, generator)

    new_branches = []
    for (i,u,v,c,model) in chosen_actions:
        b = branch.copy()
        if model == U_REPLACES_V:
            for c_i in c:
                b.add_parent(c_i-1, u, overwrite=True)
            b.add_parent(u-1, v, overwrite=True)
        else:
            b.add_parent(u-1, v, overwrite=True)

        b.add_action_info((i,u,v,c,model))
        new_branches.append(b)
    return new_branches

def sample_branches(branch, F_data, model_data, generator):
        """A general sampling function that will propose extensions of the branch passed in
        
        Parameters
        ----------
        branch : object
            an instance of the Branch class for which extensions will be proposed
        F_data : dataclass
            an instance of the FData dataclass containing all of the read count and supervariant information
            needed to compute the cellular prevalence matrix F
        model_data : dataclass
            a dataclass containing all of the information needed to setup the model
        generator : object
            a numpy default_rng object used for reproducible random sampling

        Returns
        -------
        list 
            a list of Branch objects that are unique extensions of the branch passed in that
            are sorted by their scores (log-likelihoods)
        int
            the number of unique extensions that are being returned in the possible_branches list
        """
        # initialize array of possible candidates to return
        possible_branches = []

        # propose extensions of the current branch
        proposed_branches = propose_branches(branch, F_data, model_data.branching_factor, generator)

        for new_branch in proposed_branches:

            # compute data for partial tree
            bool_mask = new_branch.parents() != NO_RELATIONSHIP
            partial_parents = compute_partial_parents(new_branch.parents())

            # compute the log-likelihood of the partial tree 
            F, eta, llh = F_llh(partial_parents, F_data, bool_mask)

            # set updated values (F matrix, eta matrix, negative log-likelihood) for partial tree
            new_branch.set_score(-llh)
            new_branch.set_eta(eta)
            new_branch.set_F(F)
            possible_branches.append(new_branch)

        return sorted(possible_branches), len(possible_branches)

def sample_branches_bs(branch, F_data, model_data, generator):
    """Samples possible nodes to add to branch for beam search
    
    Parameters
    ----------
    branch : object
        an instance of the Branch class for which extensions will be proposed
    F_data : dataclass
        an instance of the FData dataclass containing all of the read count and supervariant information
        needed to compute the cellular prevalence matrix F
    model_data : dataclass
        a dataclass containing all of the information needed to setup the model
    generator : object
        a numpy default_rng object used for reproducible random sampling

    Returns
    -------
    list 
        a list that is at most of length beam_width containing the top Branch objects sorted by their 
        scores (log-likelihoods)
    int
        the number of branches that are not being explored (i.e., the number of branches discarded)
    """
    possible_branches, cut_branches = sample_branches(branch, F_data, model_data, generator)
    return possible_branches[:model_data.beam_width], cut_branches + max(len(possible_branches) - model_data.beam_width, 0)

def sample_branches_sbs(branch, F_data, model_data, generator):
    """Samples possible nodes to add to branch for stochastic beam search
    
    Parameters
    ----------
    branch : object
        an instance of the Branch class for which extensions will be proposed
    F_data : dataclass
        an instance of the FData dataclass containing all of the read count and supervariant information
        needed to compute the cellular prevalence matrix F
    model_data : dataclass
        a dataclass containing all of the information needed to setup the model
    generator : object
        a numpy default_rng object used for reproducible random sampling

    Returns
    -------
    list 
        a list that is at most of length beam_width containing the top Branch objects sorted by their 
        scores (total perturbed log-likelihoods)
    int
        the number of branches that are not being explored (i.e., the number of branches discarded)
    """
    possible_branches, cut_branches = sample_branches(branch, F_data, model_data, generator)

    # only process branches branches if there are some, otherwise we'll obtain an error
    if len(possible_branches) == 0:
        return [], 0

    else:
        # perturb the log-likelihood using Gumbel noise and shift the scores
        gumbel_shifted_scores = gumbel_with_maximum(
                                    np.array([-new_branch.score() for new_branch in possible_branches]), 
                                    -branch.score(), 
                                    generator=generator
                                )
                            
        # overwrite scores
        for new_branch, score in zip(possible_branches, gumbel_shifted_scores):
            new_branch.set_score(-score)

        # sort and return based on beam width
        possible_branches = sorted(possible_branches)

        return possible_branches[:model_data.beam_width], cut_branches + max(len(possible_branches) - model_data.beam_width, 0)
