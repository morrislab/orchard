
########################################################################
# PathData.py
#
# This contains the source code for the PathData dataclass.
# This is simply a convenience for passing path data 
# during the initialization of Orchard.
########################################################################
from dataclasses import dataclass

@dataclass 
class PathData:
    """Class for storing path information"""
    ssm_fn: str
    params_fn: str
    results_fn: str
