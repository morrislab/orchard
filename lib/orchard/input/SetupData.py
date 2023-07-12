########################################################################
# SetupData.py
#
# This contains the source code for the SetupData dataclass.
# This is simply a convenience for passing general setup data
# during the initialization of Orchard.
########################################################################
from dataclasses import dataclass
import numpy as np

@dataclass 
class SetupData:
    """Class for storing setup parameters to run Orchard"""
    n_chains: np.int16
    poolsize: np.int16 
    randomize_nodes: bool
    seed: np.int32
