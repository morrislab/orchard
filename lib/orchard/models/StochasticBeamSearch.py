########################################################################
# StochasticBeamSearch.py
#
# A class that implements a stochastic beam search
########################################################################
import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "sampling"))

from branch_sampler import sample_branches_sbs
from BeamSearch import BeamSearch

class StochasticBeamSearch(BeamSearch):
    """Implements the stochastic beam search algorithm
    Please note that this class inherits most of its functionality from the BeamSearch/BNBModel class"""

    def _sample_branches(self, branch):
        """Wrapper around the sample_branches_bs method which samples branches using the stochastic verison of beam search
        
        Parameters
        ----------
        branch : object
            an instance of the Branch class to explore extensions of

        Returns
        -------
        list 
            a list that is at most of length beam_width containing the top Branch objects sorted by their 
            scores (total perturbed log-likelihoods)
        int
            the number of branches that are not being explored (i.e., the number of branches discarded)
        """
        return sample_branches_sbs(branch, 
                                   self._F_data,
                                   self._model_data,
                                   self._generator)


