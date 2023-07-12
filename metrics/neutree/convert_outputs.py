##########################################################################################################
# convert_outputs.py
#
# Python script for converting a zipped archive (.results.npz or .orchard.npz) to the 'Neutree' file
# format
##########################################################################################################
import os, sys, argparse
import numpy as np
from omicsdata.npz.archive import Archive

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "lib", "orchard"))

import neutree
from constants import KEYS_ORCH_NPZ, KEYS_PARAMS

def main():
  parser = argparse.ArgumentParser(
    description="Converts a zipped archive (.npz) from either Pairtree or Orchard into the 'Neutree' file format",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('npz_fn', type=str, help="Zipped archive file to convert to a 'Neutree' format")
  parser.add_argument('neutree_fn', type=str, help="Neutree file name")
  args = parser.parse_args()

  results = Archive(args.npz_fn)
  N = len(results.get(KEYS_ORCH_NPZ.struct_key))
  clusters = list(results.get(KEYS_PARAMS.clusters_key))
  ntree = neutree.Neutree(
    structs = results.get(KEYS_ORCH_NPZ.struct_key),
    phis = results.get(KEYS_ORCH_NPZ.F_key),
    counts = results.get(KEYS_ORCH_NPZ.count_key),
    logscores = results.get(KEYS_ORCH_NPZ.llh_key),
    clusterings = [clusters for idx in range(N)],
    garbage = np.array([]),
  )
  neutree.save(ntree, args.neutree_fn)

if __name__ == '__main__':
  main()