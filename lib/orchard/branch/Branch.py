########################################################################
# Branch.py
#
# This contains the source code for the Branch class that's used to 
# represent individual partial solutions.
########################################################################
import numpy as np
import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from constants import NO_RELATIONSHIP

class Branch:
    """The Branch class is used to keep track of a single partial solution.
    It's name is derived from the generalization that it's the 'branch' part of 'branch-and-bound'.
    
    Attributes
    ----------
    __parents : ndarray
        1d numpy array of size n that contains the direct parent for each node
    __score : float
        likelihood of the partial solution
    __action_sampler : object
        an instance of the Action_Sampler class which is used to sample possible node placements
    __action_info : list
        a list that is used to keep track of the steps the lead to the current state of the partial solution
    __nsamples : int
        the number of samples (features) for the dataset
    __eta : ndarray
        the subpopulation frequency matrix that is used to propose node placements
    __F : ndarray
        the cellular prevalence matrix that is used to propose node placements
    """

    def __init__(self, size, samples, actions_sampler, action_info=[]):
        """Initialize the Branch class
        
        Parameters
        ----------
        size : int
            the number of nodes (i.e., mutations or subclones) contained in the dataset 
        samples : int
            the number of samples (features) for the dataset
        action_sampler : object
            an instance of the Action_Sampler class which is used to sample possible node placements
        action_info : list
            a list that is used to keep track of the steps the lead to the current state of the partial solution
        """
        self._parents = np.full(size, NO_RELATIONSHIP)
        self._score = 0
        self.__actions_sampler = actions_sampler
        self.__action_info = action_info
        self.__nsamples = samples
        self.__eta = np.array([np.ones(samples)])
        self.__F = np.array([np.ones(samples)])

    def __hash__(self):
        """Bytes magic method"""
        return hash(self._parents.tobytes())

    def __gt__(self, branch):
        """Greater than override"""
        return self.score() > branch.score()

    def __ge__(self, branch):
        """Greater than or equal to override"""
        return self.score() >= branch.score()

    def __lt__(self, branch):
        """Less than override"""
        return self.score() < branch.score()

    def __le__(self, branch):
        """Less than or equal to override"""
        return self.score() <= branch.score()

    def __eq__(self, branch):
        """Equal to override"""
        return self.score() == branch.score()

    def __ne__(self, branch):
        """Not equal to override"""
        return self.score() != branch.score()

    def __len__(self):
        """Returns length of parented nodes"""
        return np.sum(self._parents >= 0)

    def action_info(self):
        """Returns all action information recorded"""
        return self.__action_info.copy()

    def add_action_info(self, action_info):
        """
        Allows us to keep track of what actions we take input should be where action_info = (i,u,v,model) where
        these are defined as follows:
            i: index of action selected
            u: node to add
            v: node with respect to take the action
            model: the pairwise action we take 

        The __action_info class member should end up as a vector of size N, where N is the number of nodes.
        """
        self.__action_info.append(action_info)

    def add_parent(self, node, parent, overwrite=False):
        """Adds a parent for a particular node
        
        i.e., 
        parents[node] = parent

        Parameters
        ----------
        node : int 
            the index of the node which is being 'parented'
        parent : int 
            the node number for the parent node (i.e., the index of the parent node + 1)

        Returns
        -------
        self : object
            a reference to the Branch object
        """
        if not overwrite:
            assert self._parents[node] == NO_RELATIONSHIP, "Parent already exists [%d: %d]" % (node, self._parents[node])
        assert (node) >= 0 and (node) < len(self._parents), "node %d is not a valid node for parents vector of size %d" % (node, len(self._parents))
        assert parent >= 0 and parent <= len(self._parents), "parent %d is not a valid parent for parents vector of size %d" % (parent, len(self._parents))
        
        self._parents[node] = parent

        return self

    def copy(self):
        """Returns a copy of the Branch by recursively copying objects it contains"""
        C = Branch(size=self.size(), 
                   samples=self.nsamples(), 
                   actions_sampler=self.__actions_sampler.copy(),
                   action_info=self.action_info())
        C._parents = self.parents().copy()
        C.set_eta(self.eta())
        C.set_F(self.F())
        C.set_score(self.score())
        return C
        
    def equals(self, branch):
        """Checks to see if the parents vectors are the same"""
        return all(self.parents() == branch.parents())

    def eta(self):
        """Getter for the eta"""
        return self.__eta.copy()

    def is_tree(self):
        """Returns whether or the parents vector now represents a tree"""
        return len(self) == self.size()

    def nsamples(self):
        """Returns the number of samples for the underlying data we're working with"""
        return self.__nsamples 

    def nbytes(self):
        """Returns the number of bytes for the most space-complex data members of the Branch class"""
        return self.__eta.nbytes + self.__F.nbytes + self._parents.nbytes + sys.getsizeof(self.__action_info)
        
    def parents(self):
        """Getter for the parents vector"""
        return self._parents.copy()

    def F(self):
        """Getter for the cellular prevalence matrix (F)"""
        return self.__F.copy()

    def propose_actions(self, F_data, branching_factor, generator, use_choice=False):
        """Proposes an action using the Action object held in each Branch
        
        Parameters
        ---------
        F_data : object
            an instance of the F_Data class which contains all of the read count data (var_reads, ref_reads, var_read_prob)
            for all of the mutations in the dataset
        branching_factor : int
            also known as f, the number of partial solutions to fit a cellular prevalence matrix for using the projection algorithm
        generator : object
            a numpy default_rng object used for reproducible random sampling
        use_choice : bool 
            a flag for whether or not to sample from the possible node placements. 
            If use_choice=False, we take the f placements with the highest likelihood (determined by branching_factor)
        
        
        """
        return self.__actions_sampler.propose_actions(self.parents(), F_data, self.F(), self.eta(), branching_factor, generator)

    def score(self):
        """Getter for the partial solution score"""
        return self._score

    def set_eta(self, eta):
        """Setter for the subpopulation frequency matrix (eta) fit to the partial solution
        
        Parameters
        ----------
        eta : ndarray
            matrix of subpopulation frequencies limited to the nodes in the tree
        """
        self.__eta = eta

    def set_F(self, F):
        """Setter for the cellular prevalence matrix (F) fit to the partial solution
        
        Parameters
        ----------
        eta : ndarray
            matrix of cellular prevalences limited to the nodes in the tree
        """
        self.__F = F
        
    def set_score(self, score):
        """Setter for the score of this partial solution
        
        Parameters
        -----------
        score : float
            the likelihood (log-likelihood, total perturbed log-likelihood, etc.) of the partial tree
        """
        assert np.isscalar(score), "score must be a scalar value, type %s was provided" % type(score)
        self._score = score
    
    def size(self):
        """Returns the size of the parents vector"""
        return len(self._parents)
