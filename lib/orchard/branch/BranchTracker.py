
########################################################################
# BranchTracker.py
#
# This contains the source code for the BranchTracker class.
# This class is primarily for debug purposes and characterizing 
# the Orchard algorithm.
#
########################################################################
import time

class BranchTracker:
    """Class for keeping track of exploration and cutting during search
    
    Attributes
    -----------
    __branches_explored : int 
        the number of branches (partial solutions) that have been explored by the algorithm
    __branches_in_queue : int 
        the number of branches in the priority queue 
    __branches_cut : int 
        the number of branches that have been discarded because exploring them cannot lead to an optimal solution
    __branch_nbyes : int
        the number of bytes that are contained within an instance of the Branch class
    __branch_size : int 
        the current partial solution size (i.e., the number of nodes in the partial tree)

    Parameters
    -----------
    debug : bool 
        a flag for whether or not to printout debug information
    """

    def __init__(self, debug):

        self.__debug = debug

        # variables keeping track of current branching
        self.__branches_explored = 0
        self.__branches_in_queue = 0
        self.__branches_cut = 0
        self.__branch_nbytes = 0
        self.__branch_size = 1

        self.__start_time = time.time()

    def branches_cut(self):
        """Getter for the number of branches cut"""
        return self.__branches_cut

    def branches_explored(self):
        """Getter for the number of branches explored"""
        return self.__branches_explored

    def set_branch_info(self, nbytes, size):
        """Setter for branch information
        
        Parameters
        ----------
        nbytes : int 
            the number of bytes in the Branch object
        size : int 
            the number of nodes in the partial tree
        """
        self.__branch_nbytes = nbytes 
        self.__branch_size = size + 1

    def cutting_branches(self, num_cut_branches=1):
        """Keeps track of the number of branch cut during sampling
        
        Parameters
        -----------
        num_cut_branches : int 
            the number of branches that have been discarded
        """
        self.__branches_cut += num_cut_branches

    def exploring_branch(self, expansion_factor=1, branches_in_queue=0):
        """Each time we look for new nodes to parent, we'll increment the number of branches explored
        
        Parameters
        -----------
        expansion_factor : int
            the number of branches that are being explored in parallel 
        branches_in_queue : int 
            the number of branches in the priority queue
        """
        self.__branches_explored += expansion_factor
        self.__branches_in_queue = branches_in_queue
        if self.__debug:
            self.print_exploration_info()

    def print_exploration_info(self):
        """Printout helpful information.
        The values for each of the following should be:
            Num Branches In Queue  = expansion_factor*beam_width
            Num Branches Explored  = l*expansion_factor
            Num Branches Cut       = l*(branching_factor*expansion_factor - beam_width) 

        where l is the number of nodes that have been added to the tree. If multiple cores are being used,
        these printouts will occur simultaneously. It's best to use debug mode when --num-cpu-cores=1.
        """
        elapsed_seconds = time.time() - self.__start_time
        print("----------------------------------------------------------------------------------------------------------------------")
        print("Time (seconds) | Num Branches In Queue | Num Branches Explored | Num Branches Cut | Estimated KB/branch| Branch size")
        print("    %d                 %d                      %d                     %d                   %.4f             %d" 
                %  (elapsed_seconds, self.__branches_in_queue, self.__branches_explored, self.__branches_cut, self.__branch_nbytes/1e3, self.__branch_size))
        print("----------------------------------------------------------------------------------------------------------------------")
        print()
