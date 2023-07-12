#!/bin/bash

###############################################################################################################
# compare_rrl.sh
#
# Compares the tree structures produced by a method to the ground truth tree structure.
#
###############################################################################################################

# initialize data directory
data_dir=$DATA_DIR

if [ ! -z "$1" ]
  then
    data_dir=$1
fi

####################################
# (1) Evaluate rels
####################################
for dir in $data_dir/* ; do

    runid=$(basename $dir | cut -d. -f1)
    truthdir=$dir/truth
    args=""
    for foldername in "${@:2}"
    do
        args="${args} ${foldername}=${dir}/${foldername}/${runid}.neutree.npz"
    done

    python3 $ORCH_DIR/metrics/relationship_reconstruction_loss.py truth=$truthdir/$runid.neutree.npz $args > $dir/$runid.rrl.txt
        
done
