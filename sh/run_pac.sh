#!/bin/bash

###############################################################################################################
# run_pac.sh
#
# Expects the following command line arguments:
#
# $1 = data_dir, directory that the data is in, which has a required format
# $2 = orchard_dir_name, the name of the folder that contains Orchard's reconstruction results
#
# FLAGS
# intruthdir, flag that means that the .ssm and .params.json files are in a /$dir/truth folder
# convertoutput, flag that means the .orchard.npz file needs to be converted to a neutree format
#
# Performs the following:
#   (1) Converts each $dir/$orchard_dir_name/*.orchard.npz file to a neutree file (optional)
#   (2) Runs phylogeny-aware clustering on each $dir/$orchard_dir_name/*.neutree.npz file
#
###############################################################################################################


##################################
# Process command line arguments #
##################################
# set data directory
# we expect each folder in the data_dir to have a truth folder containing a .ssm and .params.json file
data_dir=$DATA_DIR

if [ ! -z "$1" ]
  then
    data_dir=$1
fi

# set folder to store orchard results in
orchard_dir_name=orchard

if [ ! -z "$2" ]
  then
    orchard_dir_name=$2
fi

# flag to signal that the .ssm and .params.json files are in a folder named 'truth'
intruthdir=false
if [[ "$*" == *"intruthdir"* ]]
then
    intruthdir=true
fi

# flag to signal whether or not the tree outputs (npz) need to be convert to the neutree format
convertoutput=false
if [[ "$*" == *"convertoutput"* ]]
then
    convertoutput=true
fi

##############################
# (1) Convert Outputs
##############################

if [[ $convertoutput = true ]]
then
  # go through each folder and convert the .orchard.npz file to a neutree format
  for dir in $data_dir/* ; 
  do
      
      runid=$(basename $dir | cut -d. -f1)
      
      orcharddir=$dir/$orchard_dir_name    
      resultsfn=$runid.orchard.npz
      
      python3 $ORCH_DIR/metrics/neutree/convert_outputs.py \
                          $orcharddir/$resultsfn \
                          $orcharddir/$runid.neutree.npz 


  done
fi

##############################
# (2) Run clustering
##############################
for dir in $data_dir/* ; 
do
    
    runid=$(basename $dir | cut -d. -f1)
    
    orcharddir=$dir/$orchard_dir_name    

    datapath=$dir
    if [[ $intruthdir = true ]]
    then
      datapath=$dir/truth
    fi 
    
    python3 $ORCH_DIR/bin/pac \
                        $datapath/$runid.ssm \
                        $datapath/$runid.params.json \
                        $orcharddir/$runid.neutree.npz \
                        $orcharddir/cluster
                        
    python3 $ORCH_DIR/lib/cluster/generate_clonetree.py $datapath/$runid.ssm $datapath/$runid.params.json $orcharddir/cluster/clusters.npz $orcharddir/cluster/cluster.params.json $orcharddir/cluster/cluster.results.npz

    python3 $ORCH_DIR/lib/cluster/params_to_clusters.py $orcharddir/cluster/cluster.params.json $orcharddir/clusters.csv

done