#!/bin/bash

###############################################################################################################
# run_orchard.sh
#
# Expects the following command line arguments:
#
# $1 = data_dir, directory that the data is in, which has a required format
# $2 = orchard_dir_name, the name of the folder to store the results from Orchard in
# $3 = num_cores, then number of parallel instances of Orchard to run
# $4 = beam_width, the priority queue size
# $5 = branching_factor, the number of placements to try for puting node u into a partial tree
#
# FLAGS
# random, randomizes the node order
# beam_search, performs a beam search (deterministic search instead of stochastic)
# monoprimary, force Orchard to only search for monoprimary trees
# ignorezeroprobs, tells orchard to ignore node placements with zero probability during search
# intruthdir, signals that the .ssm and .params.json files are in a /$dir/truth folder
#
# Performs the following:
#   (1) Runs orchard on each folder in the $DATA_DIR 
#   (2) Converts the output from orchard to a format that can be compared to other reconstructions
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

# set number of parallel instances (cores) to use
num_cores=1

if [ ! -z "$3" ]
  then
    num_cores=$3
fi

# set the beam width, i.e., the number of unique trees to return per parallel instance
beam_width=1

if [ ! -z "$4" ]
  then
    beam_width=$4
fi

# set the branching factor, i.e., the number of placements to try when adding a node to the tree
branching_factor=20

if [ ! -z "$5" ]
  then
    branching_factor=$5
fi

# flag to randomize the node order
random=""

if [[ "$*" == *"random"* ]]
then
    random="-r"
fi

# flag to tell orchard to run as beam search
beam_search=""

if [[ "$*" == *"beamsearch"* ]]
then
    beam_search="--model=bs"
fi

# flag to tell orchard to force monoprimary trees
monoprimary=""
if [[ "$*" == *"monoprimary"* ]]
then
    monoprimary="-p"
fi

# flag to tell orchard to not consider node placements with a zero probability during search
ignorezeroprobs=""
if [[ "$*" == *"ignorezeroprobs"* ]]
then
    ignorezeroprobs="-z"
fi

# flag to signal that the .ssm and .params.json files are in a folder named 'truth'
intruthdir=false
if [[ "$*" == *"intruthdir"* ]]
then
    intruthdir=true
fi

##############################
# (1) Run Orchard
##############################
for dir in $data_dir/* ; 
do
    
    runid=$(basename $dir | cut -d. -f1)
    
    orcharddir=$dir/$orchard_dir_name

    mkdir $orcharddir
    
    datapath=$dir
    if [[ $intruthdir = true ]]
    then
      datapath=$dir/truth
    fi 

    ssmfn=$datapath/$runid.ssm
    paramsfn=$datapath/$runid.params.json
    resultsfn=$orcharddir/$runid.orchard.npz
    
    time ( python3 $ORCH_DIR/bin/orchard $ssmfn $paramsfn $resultsfn -n $num_cores -f $branching_factor -k $beam_width $random $beam_search $ignorezeroprobs $monoprimary ) \
           2>> $orcharddir/$runid.orchard.out 

done


##############################
# (2) Convert Outputs
##############################
for dir in $data_dir/* ; 
do
    
    runid=$(basename $dir | cut -d. -f1)
    
    orcharddir=$dir/$orchard_dir_name    
    resultsfn=$runid.orchard.npz
    
    python3 $ORCH_DIR/metrics/neutree/convert_outputs.py \
                        $orcharddir/$runid.orchard.npz \
                        $orcharddir/$runid.neutree.npz 


done

