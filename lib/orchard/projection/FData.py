########################################################################
# FData.py
#
# This contains the source code for the FData dataclass.
# This is simply a convenience for passing all of the read count
# and supervariant information.
########################################################################
from dataclasses import dataclass, field

@dataclass 
class FData:
    """Class for storing data needed to calculate cellular prevalences for a given tree"""
    superclusters: list
    V: list
    N: list 
    omega: list

    supervars: dict = field(default_factory=dict)
