Orchard
========

Orchard uses bulk DNA sequencing data to quickly infer accurate cancer phylogenies.

Table of Contents
=================
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Examples](#examples)
4. [Orchard Parameters](#orchard_parameters)
5. [Orchard Inputs](#orchard_inputs)
6. [Orchard Output](#orchard_output)


Requirements
================
Python 3+


Installation
===============

Clone repository
```
git clone https://github.com/morrislab/orchard && cd orchard
```

Setup conda environment with python 3.8 or higher 

```
conda create -y --name orchard python=3.10


CONDA_ACTIVATE_ENV_VARS=$HOME/.conda/envs/orchard/etc/conda/activate.d/env_vars.sh
CONDA_DEACTIVATE_ENV_VARS=$HOME/.conda/envs/orchard/etc/conda/deactivate.d/env_vars.sh

mkdir -p $HOME/.conda/envs/orchard/etc/conda/activate.d
mkdir -p $HOME/.conda/envs/orchard/etc/conda/deactivate.d

touch $CONDA_ACTIVATE_ENV_VARS
touch $CONDA_DEACTIVATE_ENV_VARS

echo "" >> $CONDA_ACTIVATE_ENV_VARS
echo "# Set Environment Variables" >> $CONDA_ACTIVATE_ENV_VARS
echo "export ORCH_DIR=$PWD" >> $CONDA_ACTIVATE_ENV_VARS

echo "" >> $CONDA_DEACTIVATE_ENV_VARS
echo "# Deactivate Environment Variables" >> $CONDA_DEACTIVATE_ENV_VARS
echo "unset ORCH_DIR" >> $CONDA_DEACTIVATE_ENV_VARS

conda activate orchard
```

OR Setup virtual environment

```
python3 -m venv env
echo "" >> env/bin/activate
echo "# Environment Variables" >> env/bin/activate
echo "export ORCH_DIR=$PWD" >> env/bin/activate
source env/bin/activate
```

Install dependencies

```
python -m pip install -r requirements.txt

cd $ORCH_DIR/lib
git clone https://github.com/ethanumn/projectppm
cd projectppm
bash make.sh
cd $ORCH_DIR
```

Download Pairtree (optional)

```
cd $ORCH_DIR/lib
git clone https://github.com/jwintersinger/pairtree && cd pairtree
pip3 install -r requirements.txt
cd lib
git clone https://github.com/jwintersinger/projectppm
cd projectppm
bash make.sh
cd $ORCH_DIR
```

Install libraries for simulating mutation trees (optional)
```
cd $ORCH_DIR/lib
git clone https://github.com/ethanumn/mutsim.git
```

Download simulated datasets from Pairtree paper (optional)
```
cd $ORCH_DIR
git clone https://github.com/morrislab/pairtree-experiments.git
$SH_DIR/reorganize_folder.sh $ORCH_DIR/pairtree-experiments/inputs/sims.smallalpha.pairtree
```


Examples
===========

Run Orchard on example data

```
example1_dir=$ORCH_DIR/examples/example1
python3 $ORCH_DIR/bin/orchard $example1_dir/example1.ssm $example1_dir/example1.params.json $example1_dir/example1.orchard.npz
```

Run Orchard on a directory, where each subdirectory contains an SSM file and parameter file
```
$ORCH_DIR/sh/run_orchard.sh $ORCH_DIR/examples
``` 


Orchard Parameters
===================

There are quite a few different parameters that Orchard has. We provide a brief description of all of the parameters.

- `ssm_fn` : str, required
    - the file path to the simple somatic mutation file (see [Orchard Inputs](#orchard-inputs) for file format details)
- `params_fn` : str, required
    - the file path to a parameter file (see [Orchard Inputs](#orchard-inputs) for file format details)
- `results_fn` : str, required
    - the file name to write a zipped archive to (see the [omicsdata package](https://omicsdata.readthedocs.io/en/latest/brief.html#writing-data-archives))
- `-e`, `--expansion-factor` : int, default=1
    - the number of partial solutions to expand at each iteration of the algorithm. The default `expansion_factor` is 1, which results in a **best first search**. If the `expansion_factor=beam_width`, this results in a **breadth first search**. Increasing the `expansion_factor` can result in large run times.
- `-f`, `--branching-factor`: int, default=20
    - the maximum number of node placements to consider when placing a new node into a partial tree
-  `-k`, `--beam-width` : int, default=1
    - the maximum of the priority queue. This limits the number of partial trees to keep track of during search
- `-c`, `--n-chains` : int, default=None
    - the number of parallel instances of Orchard's sampling routine to run. If None, it will
    default to the value of `num_cpu_cores`
- `-d`, `--debug` : bool, default=False
    - a flag to printout useful information for debugging
- `-m`, `--model` : str, default=sbs
    - denotes which type of model to run, sbs means that Orchard will run a **stochastic beam search**, bs means that Orchard will run a **deterministic beam search**
- `-n`, `--num-cpu-cores` : int, default=cpu_count()
    - the number of cores to provide Orchard with. If `num_cpu_cores=1`, then orchard will
    run as a single thread.
- `-p`, `--force-monoprimary` : bool, default=False
    - flag to tell orchard to only search for trees that are monoprimary
- `-r`, `--randomize-nodes` : bool, default=False
    - flag to tell orchard to run using a randomized node order
- `-s`, `--seed` : int, default=None
    - seed value for duplicating orchard's results across multiple runs
- `-x`, `--max-placements` : int, default=10
    - the maximum number of placements that can be tried for some node `u`, with respect to a single node `v`. This is not the same as the `branching_factor`. This parameter should not be changed as it's only applicable when the phylogeny may be a star tree
- `-z`, `--ignore-zero-probs` : bool, default=False
    - a flag to ignore node placements of `u` that have a zero probability. Setting this flag may increase runtimes, however, for most reconstruction problems this will only 
    provide a small decrease in runtime

Orchard Inputs
================

Orchard expects two input files to perform its reconstruction. A **simple somatic mutation file**, and a **parameters file**.

Simple Somatic Mutation file (SSM)
-----------------------------------

The SSM file format has been thoroughly detailed on the [omicsdata readthedocs](https://omicsdata.readthedocs.io/en/latest/brief.html#simple-somatic-mutation-file). This file format has also been detailed in the [Pairtree STAR protocol](https://www.sciencedirect.com/science/article/pii/S266616672200586X#sectitle0040) and the [Pairtree GitHub](https://github.com/morrislab/pairtree#input-files).


Here are the first few lines of the `example1` simple somatic mutation file
```
id	name	var_reads	total_reads	var_read_prob
s0	S_0	30,20,22	50,50,50	0.5,0.5,0.5
s1	S_1	0,9,0	50,50,50	0.5,0.5,0.5
s2	S_2	28,1,19	50,50,50	0.5,0.5,0.5
```

Parameters File
-----------------

The parameters file format has also been thoroughly detailed on the [omicsdata readthedocs](https://omicsdata.readthedocs.io/en/latest/brief.html#parameter-file). This file format has also been detailed in the [Pairtree STAR protocol](https://www.sciencedirect.com/science/article/pii/S266616672200586X#p0100) and the [Pairtree GitHub](https://github.com/morrislab/pairtree#input-files).

Here are the most relevant key/value pairs in the `example1` parameters file 
```
{"samples": ["Sample 1", "Sample 2", "Sample 3"], "clusters": [["s9"], ["s0"], ["s1"], ["s2"], ["s4"], ["s3"], ["s7"], ["s5"], ["s8"], ["s6"]]}
```


Orchard Outputs
=================

Orchard outputs a zipped archive containing `numpy` data types. It uses the [Archive class](https://omicsdata.readthedocs.io/en/latest/brief.html#writing-data-archives) in the `omicsdata package` to compress its data outputs. This format is the same as that which is outputted by [Pairtree](https://github.com/morrislab/pairtree/).

The only output by Orchard is a single zipped archive file `*.npz`. The contents of this file can be viewed using the `numpy.load` function (see the [numpy.load docs](https://numpy.org/doc/stable/reference/generated/numpy.load.html)) as follows 

```
# assumes we ran orchard on the example data
import numpy as np 
example_npz_fn = "/path/to/example1.orchard.npz" 
data = np.load(example_npz_fn)
```

Alternatively, the zipped archive output by Orchard can be opened using the `Archive` class 

```
from omicsdata.npz.archive import Archive
example_npz_fn = "/path/to/example1.orchard.npz" 
archive = Archive(example_npz_fn)
parents = archive.get("struct") # get parents list from archive
```

The archive output by Orchard contains `numpy` data types stored under specific keys values. The values for the keys `struct`, `count`, `phi`, `llh`, `prob` should be sorted lists of the same length, and they are sorted such that the same index in each of the lists corresponds to the data for a specific tree. The lists are sorted in descending order according to the log-likelihood of the tree, i.e., the data at index 0 corresponds to the data for the tree with the lowest log-likelihood (the best tree). Here, we describe the key/value pairs contained within the archive output by Orchard: 

* `struct`: a list of `parents` vectors sorted by the likelihood of their corresponding cellular prevalence matrix under a binomial likelihood model.

* `count`: a list denoting the number of times each `parents` vector was found during search

* `phi`: a list of 2D `numpy` arrays, where each 2D array is the cellular prevlaence matrix fit to the tree.

* `llh`: a list of log-likelihoods, where each value is the log of the likelihood of the cellular prevalence matrix in `phi` under a binomial likelihood model.

* `prob`: a list containing the probability that the tree is the best of those that were found during search.

