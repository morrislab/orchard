########################################################################
# BNBModel.py
#
# A base class for branch-and-bound models.
########################################################################
import sys, os
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "branch"))
from BranchTracker import BranchTracker

class BNBModel:
    """Simple base class for a branch and bound style model
    
    
    Attributes
    ----------
    __branch_tracker : object
        an instance of the BranchTracker class which is for debug purposes
    __data_queue : object
        a python queue class or a class that inherits from the python queue class 
    __generator : object
        a numpy default_rng object used for reproducible random sampling
    __model_data : dataclass
        a dataclass containing the information for setting up the BNBModel
    __phi_data : dataclass
        a dataclass containing the read count data for all supervariants
    __progress_queue : object
        a thread safe python queue 
    """

    def __init__(self, init_branches, model_data, phi_data, generator, data_queue, progress_queue=None):
        """Initializes member variables
        
        Parameters
        ----------
        init_branches : list
            a list of Branch objects that can be used to provide a BNB class (or Orchard) with arbitrary starting branches
        model_data : dataclass
            a dataclass containing the information for setting up the BNBModel
        phi_data : dataclass
            a dataclass containing the read count data for all supervariants
        generator : object
            a numpy default_rng object used for reproducible random sampling
        data_queue : object
            a python queue class or a class that inherits from the python queue class 
        progress_queue : object
            a thread safe python queue 
        """
        self._branch_tracker = BranchTracker(debug=model_data.debug)
        self._data_queue = data_queue
        self._generator = generator
        self._model_data = model_data
        self._phi_data = phi_data
        self._progress_queue = progress_queue

        self._best_trees = []
        self._unique_best_trees = set()
        self._lower_bound = np.inf

        self.add_branches(init_branches)

    def add_branch(self, branch):
        """Adds a single branch to the priority queue
        
        Parameters
        ----------
        branch : object
            an instance of the Branch class
        """
        self._put(branch)

    def add_branches(self, branches):
        """Adds a list of branches to the priority queue
        
        Parameters
        ----------
        branches : list
            a list of Branch objects
        """
        for branch in branches:
            self.add_branch(branch)

    def best_trees(self):
        """Returns a list of the best trees found during search"""
        return sorted(self._best_trees.copy())

    def branches_cut(self):
        """Returns the number of branches cut"""
        return self._branch_tracker.branches_cut()

    def branches_explored(self):
        """Returns the number of branches explored"""
        return self._branch_tracker.branches_explored()

    def lower_bound(self):
        """Returns the lower bound"""
        return self._lower_bound

    def model_data(self):
        """Returns the dataclass used to initial the model type"""
        return self._model_data

    def search(self):
        """Search method must be overrided in subclass"""
        raise NotImplementedError

    def _add_all_trees(self, trees):
        """Adds all trees to the queue
        
        Parameters
        ----------
        trees : list
            a list of Branch objects each of which is a complete tree
        """
        n_unique_trees = len(self._unique_best_trees)
        for tree in trees:
            assert tree.is_tree(), "Object of type %s is not a valid tree" % type(tree)

            self._best_trees.append(tree)
            # add hash of parents vector to list of uniques so we can determine
            # if we've found a new unique tree and can decrease our heap size
            self._unique_best_trees.add(hash(tree))

        if len(self._unique_best_trees) > n_unique_trees:
            self._data_queue.decrease_heapsize()
            # update progress queue if it exists
            if self._progress_queue is not None:
                self._progress_queue.put(1)

    def _add_best_tree(self, tree):
        """Adds a tree to the list of our best trees found so far
        
        Parameters
        ----------
        tree : object
            an instance of a Branch object that is a complete tree
        """
        assert tree.is_tree(), "Object of type %s is not a valid tree" % type(tree)
        self._best_trees.append(tree)
        self._data_queue.decrease_heapsize()
        # update progress queue if it exists
        if self._progress_queue is not None:
            self._progress_queue.put(1)

    def _evaluate_tree(self, tree):
        """Evaluates a tree to see if it matches or improves upon our lower bound
        
        Parameters
        ----------
        tree : object
            an instance of a Branch object that is a complete tree
        """
        if self._within_bound(tree):
            self._set_lower_bound(tree.score())
            self._add_best_tree(tree)

    def _evaluate_branch(self, branch):
        """Evaluates a branch to see if it improves upon our lower bound, or 
        is a partial solution and may be able to improve upon our lower bound
        
        Parameters
        ----------
        branch : object
            an instance of a Branch object 
        """
        if branch.is_tree():
            self._evaluate_tree(branch)
        else:
            if branch.nllh() <= self.lower_bound():
                self._put(branch)
            else:
                self._branch_tracker.cutting_branches()

    def _get(self):
        """Wrapper around whatever data structure used to get branches"""
        return self._data_queue.get()

    def _put(self, branch):
        """Wrapper around whatever data structure used to store branches
        
        Parameters
        ----------
        branch : object
            an instance of a Branch object 
        """
        self._data_queue.put(branch)

    def _set_lower_bound(self, lower_bound):
        """Setter for the lower bound
        
        Parameters
        ----------
        lower_bound : float
            a lower bound on the likelihood for the optimal solution
        """
        self._lower_bound = lower_bound

    def _within_bound(self, branch):
        """Checks if the negative log-likelihood of the branch is still within
        our current lower bound
        
        Parameters
        ----------
        branch : object
            an instance of a Branch object 
        """
        return branch.nllh() <= self.lower_bound()