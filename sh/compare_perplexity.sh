#!/bin/bash

###############################################################################################################
# compare_perplexity.sh
#
# Computes the perplexity of the cellular prevalence matrix for one or more neutree files.
#
###############################################################################################################

# initialize data directory
data_dir=$DATA_DIR

if [ ! -z "$1" ]
  then
    data_dir=$1
fi

####################################
# (1) Evaluate perplexity
####################################
for dir in $data_dir/* ; do

    runid=$(basename $dir | cut -d. -f1)
    truthdir=$dir/truth

    args=""
    for foldername in "${@:2}"
    do
        args="${args} ${foldername}=${dir}/${foldername}/${runid}.neutree.npz"
    done

    python3 $ORCH_DIR/metrics/perplexity.py --ssm-fn $truthdir/$runid.ssm \
        truth=$truthdir/$runid.neutree.npz $args \
        > $dir/$runid.logperplexity.txt
        
done
