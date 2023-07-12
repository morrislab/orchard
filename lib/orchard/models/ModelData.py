########################################################################
# ModelData.py
#
# A dataclass that contains the information necessary to run either
# beam search or stochastic beam search.
########################################################################
from dataclasses import dataclass
from numpy import bool_, int16

@dataclass 
class ModelData:
    """Class for storing model parameters"""
    beam_width: int16 
    branching_factor: int16
    debug: bool_
    ignore_zero_probs: bool_
    expansion_factor: int16
    force_monoprimary: bool_
    max_placements: int16
