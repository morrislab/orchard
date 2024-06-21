###############################################################################################################
# constants.py
# 
# Constants that are held inside of python dataclasses for convenience
###############################################################################################################

from dataclasses import dataclass

##############################
# Orchard Keys
##############################

# .npz keys
@dataclass 
class KEYS_ORCH_NPZ:
  action_info: str = "action_info"
  branches_explored: str = "branches_explored"
  branches_cut: str = "branches_cut"
  count_key: str = "count"
  eta_key: str = "eta"
  llh_key: str = "llh"
  prob_key: str = "prob"
  F_key: str = "phi" # keep the actual string key that same as Pairtree for compatibility
  sampnames: str = "sampnames"
  seed: str = "seed"
  struct_key: str = "struct"
  clusters_key: str = "clusters"
  newick_key: str = "newick"

  # just for compatibility with pairtree  
  clustrel_posterior_rels: str = "clustrel_posterior_rels"
  clustrel_posterior_vids: str = "clustrel_posterior_vids"
  garbage: str = "garbage" 

# params.json keys
@dataclass
class KEYS_PARAMS:
  clusters_key: str = "clusters"
  garbage_key: str = "garbage"
  samples_key: str = "samples"

# .ssm keys
@dataclass 
class KEYS_SSM:
  id_key: str = "id"
  name_key: str = "name"
  var_reads_key: str = "var_reads"
  total_reads_key: str = "total_reads"
  var_read_prob_key: str = "var_read_prob"

##############################
# Common definitions
##############################
# these can't be placed in a dataclass because Numba doesn't recognize dataclasses
MAX_SEED = 2**32 - 1
NO_RELATIONSHIP = -1
MISSING_REL_PROB = 1e-14
MIN_VARIANCE = 1e-4

##############################
# Initialization schemes
##############################
F_SUM_NODE_ORDER = "fsum"
RANDOM_NODE_ORDER = "random"
DIVERSE_NODE_ORDER = "diverse"
node_order_options = (F_SUM_NODE_ORDER, RANDOM_NODE_ORDER, DIVERSE_NODE_ORDER)