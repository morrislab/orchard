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
    num_instances: np.int16
    poolsize: np.int16 
    node_order: str
    seed: np.int32
