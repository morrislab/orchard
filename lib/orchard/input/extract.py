########################################################################
# extract.py
#
# This contains the source code for preparing the read count data 
# for Orchard.
########################################################################
import os, sys
import numpy as np 
from omicsdata.ssm import parse, supervariants, columns

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "projection"))

# orchard imports
from FData import FData

def extract_F_and_params(ssm_fn, params_fn):
    """Extracts the data from the simple somatic mutation file and parameters file 
    needed to run Orchard
    
    Parameters
    -----------
    ssm_fn : str
        the path for the simple somatic mutation file
    params_fn : str
        the path of the parameters file

    Returns
    --------
    F_data : object
        an instance of the FData class which packages together all of the read count data used
        by Orchard
    params : dictionary
        a dictionary containing the key/value pairs directly from the parameters file 
    """

    # obtain data for entire dataset
    variants = parse.load_ssm(ssm_fn)
    params = parse.load_params(params_fn)
    clusters = params[columns.PARAMS_Columns.CLUSTERS]

    # use pairtree utilities to extract mutation data 
    supervars = supervariants.clusters_to_supervars(clusters, variants)
    superclusters = supervariants.make_superclusters(supervars)
    V, R, omega_v = supervariants.supervars_to_binom_params(supervars)

    # package mutation data for later use
    F_data = FData(supervars=supervars,
                   superclusters=superclusters,
                   V=V,
                   N=V+R,
                   omega=omega_v)

    return F_data, params