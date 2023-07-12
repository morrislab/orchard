########################################################################
# BeamSearch.py
#
# A class that implements a classical beam search.
########################################################################
import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "queue"))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "branch"))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "sampling"))

from BNBModel import BNBModel
from Branch import Branch
from UniquePriorityQueue import UniquePriorityQueue
from branch_sampler import sample_branches_bs


class BeamSearch(BNBModel):
    """Beam search class -- note that it inherits most of its functionality from the BNBModel class
    
    Attributes
    ----------
    _beam_width : int 
        the number of Branch objects (partial solutions) to keep track of at each iteration
    _branching_factor : int 
        the number of node placements to fit the cellular prevalence matrix for (F) and score
    """
    def __init__(self, init_branches, bs_data, F_data, generator, progress_queue=None):
        """Initializes member variables for beam search
        
        Parameters
        ----------
        init_branches : list
            a list of Branch objects that can be used to provide a BNB class (or Orchard) with arbitrary starting branches
        bs_data : dataclass
            a dataclass containing the information for setting up beam search
        F_data : dataclass
            a dataclass containing read count data used to compute the cellular prevalence matrix (F)
        generator : object
            a numpy default_rng object used for reproducible random sampling
        progress_queue : object
            a thread safe python queue 
        """
        self._beam_width = bs_data.beam_width
        self._branching_factor = bs_data.branching_factor
        
        super().__init__(self._init_branches(init_branches), 
                         bs_data, 
                         F_data,
                         generator,
                         data_queue=UniquePriorityQueue(maxsize=0, heapsize=bs_data.beam_width, expansion_factor=bs_data.expansion_factor),
                         progress_queue=progress_queue)


    def beam_width(self):
        """Getter for the beam width"""
        return self._beam_width

    def best_trees(self):
        """Returns a list of the best trees found during search"""
        return sorted(self._best_trees.copy())[:self.beam_width()]

    def branching_factor(self):
        """Getter for the branching factor"""
        return self._branching_factor

    def search(self):
        """Performs beam search"""
        while not self._data_queue.empty():
            # when we call _data_queue.get_m(), we cut all but the number of branches 
            # defined by our beam_width from the priority queue
            self._branch_tracker.cutting_branches(max(0, self._data_queue.qsize() - self.beam_width()))

            # used for tracking search space exploration
            self._branch_tracker.exploring_branch(self._data_queue.expansion_factor(), self._data_queue.qsize())

            # explore m branches, where m is our expansion factor
            for branch in self._data_queue.get_m(): 

                # Conditional for determining if:
                #   (1) branch is invalid
                #   (2) branch should be explored
                if not isinstance(branch, Branch):
                    continue
                else:

                    # compute valid branches and sample k of them
                    sampled_branches, cut_branches = self._sample_branches(branch)

                    # we add all configurations if the sampled_branches are complete trees, 
                    # and we'll decide later which are the top-k
                    if all([branch.is_tree() for branch in sampled_branches]):
                        self._add_all_trees(sampled_branches)
                    else:
                        # otherwise, we'll add the branches back to the queue
                        for new_branch in sampled_branches:
                            self._put(new_branch)

                    self._branch_tracker.cutting_branches(cut_branches)
                    self._branch_tracker.set_branch_info(branch.nbytes(), len(branch))


    def _choose_branches(self, branches):
        """Returns the top k branches
        
        Parameters
        ----------
        branches : list
            a list of Branch objects
        """
        return branches[:self.beam_width()], max(len(branches)-self.beam_width(),0)

    def _init_branches(self, branches):
        """Given a list of initial branches to search, sample those to search on.
        Only issue with this is we do not account for these branches as being cut.
        
        Parameters
        ----------
        branches : list
            a list of Branch objects
        """
        init_branches, _ = self._choose_branches(sorted(branches))
        return init_branches

    def _sample_branches(self, branch):
        """Wrapper around the sample_branches_bs method which samples branches using beam search
        
        Parameters
        ----------
        branch : object
            an instance of a Branch class
        """
        return sample_branches_bs(branch, 
                                  self._F_data, 
                                  self._model_data,
                                  self._generator)
