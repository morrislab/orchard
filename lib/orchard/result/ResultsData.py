
########################################################################
# ResultsData.py
#
# This contains the source code for the ResultsData dataclass.
# This is simply a convenience for packaging the data needed
# to compute the final results for Orchard.
########################################################################
from dataclasses import dataclass, field

@dataclass 
class ResultsData:
    """Class for storing and transporting the results from running SuperMT"""
    best_trees: list 
    branches_explored: int
    branches_cut: int
    seed: int 
    results_fn: str 

    params: dict = field(default_factory=dict)
