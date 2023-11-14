########################################################################
# params_to_clusters.py
#
# Translates the 'clusters' key in a .params.json file to a CSV file
# that contains key/value pairs of {'variant id':'cluster assignment'}
########################################################################
import argparse, os, sys  
import pandas as pd
from omicsdata.ssm import parse

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', "orchard")))
from constants import KEYS_ORCH_NPZ

def main():
    """Main function to translate params 'cluster' key to a csv file"""

    parser = argparse.ArgumentParser(
        description="Translates outputted params file to a csv file containing id,cluster#",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("params_fn", type=str, help="Parameters file (.params.json)")
    parser.add_argument("output_csv", type=str, help="Csv file to output clusters to")

    args = parser.parse_args()

    # extract 'cluster' key from params file
    params = parse.load_params(args.params_fn)
    clusters = params[KEYS_ORCH_NPZ.clusters_key]

    ids_ = [] 
    clusters_ = []

    # collect necessary cluster info for each id
    for i in range(len(clusters)):
        for vid in clusters[i]:
            ids_.append(vid)
            clusters_.append(i)

    # write it to a csv file
    df = pd.DataFrame({"id":ids_, "clusters":clusters_})
    df["sort_col"] = df["id"].apply(lambda x: int(x[1:]))
    df = df.sort_values(by=["sort_col"])
    df = df.drop("sort_col", axis=1)
    df.to_csv(args.output_csv, header=df.columns, index=False)

if __name__ == '__main__':
    main()

