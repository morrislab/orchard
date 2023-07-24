Orchard
========

Orchard uses bulk DNA sequencing data to quickly infer accurate cancer phylogenies.

Table of Contents
=================
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Examples](#examples)
4. [FAQ](#faq)
5. [Orchard Parameters](#orchard-parameters)
6. [Running Orchard](#running-orchard)
7. [Orchard Inputs](#orchard-inputs)
8. [Orchard Outputs](#orchard-outputs)
9. [Phylogeny-aware clustering](#phylogeny-aware-clustering)



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

Download repository to replicate the results from the Orchard manuscript  (optional)
```
cd $ORCH_DIR
git clone https://github.com/ethanumn/orchard_experiments.git
cd $ORCH_DIR/orchard_experiments
```

Install libraries for simulating mutation trees (optional)
```
cd $ORCH_DIR/lib
git clone https://github.com/ethanumn/mutsim.git
```

Examples
===========

Run Orchard on example data (`--seed` is used to replicate results)

```
example1_dir=$ORCH_DIR/examples/example1
python3 $ORCH_DIR/bin/orchard $example1_dir/example1.ssm $example1_dir/example1.params.json $example1_dir/example1.orchard.npz --seed=123
```

Run Orchard on a directory, where each subdirectory contains an SSM file and parameter file
```
$ORCH_DIR/sh/run_orchard.sh $ORCH_DIR/examples
``` 

Visualize the best tree output by Orchard
```
example1_dir=$ORCH_DIR/examples/example1
python3 $ORCH_DIR/bin/plot $example1_dir/example1.orchard.npz $example1_dir/example1.plot.html --ssm-fn $example1_dir/example1.ssm
open $example1_dir/example1.plot.html
```

Run phylogeny-aware clustering to infer a clone tree from a mutation tree (run Orchard on example data prior to running this)
```
example1_dir=$ORCH_DIR/examples/example1

# convert example1 data to neutree format
python3 $ORCH_DIR/metrics/neutree/convert_outputs.py $example1_dir/example1.orchard.npz $example1_dir/example1.neutree.npz

# run phylogeny-aware clustering on our example data
python3 bin/phylogeny_aware_clustering $example1_dir/example1.ssm $example1_dir/example1.params.json $example1_dir/example1.neutree.npz examples/example1/cluster

# use the generalized information criterion to select a clone tree
python3 lib/cluster/generate_clonetree.py $example1_dir/example1.ssm $example1_dir/example1.params.json $example1_dir/cluster/clusters.npz $example1_dir/cluster/cluster.params.json $example1_dir/cluster/cluster.results.npz -p

# visualize the best clone tree
# the clone tree should have 9 nodes if Orchard was given the seed 123
python bin/plot $example1_dir/cluster/cluster.results.npz $example1_dir/cluster/cluster.plot.html --ssm-fn $example1_dir/example1.ssm

open $example1_dir/cluster/cluster.plot.html 
```

Compare the trees output by Orchard to a simulated ground truth
```
# make sure we've run Orchard on the example1 data prior to running this

# convert simulated data .pickle file to .neutree.npz
python3 $ORCH_DIR/metrics/neutree/sim_to_neutree.py $example1_dir/example1.truth.pickle $example1_dir/example1.truth.neutree.npz

# convert .orchard.npz to neutree.npz for comparison
python3 $ORCH_DIR/metrics/neutree/convert_outputs.py $example1_dir/example1.orchard.npz $example1_dir/example1.orchard.neutree.npz

# compute relationship reconstruction loss, and save it to text file
python3 $ORCH_DIR/metrics/relationship_reconstruction_loss.py truth=$example1_dir/example1.truth.neutree.npz orchard=$example1_dir/example1.orchard.neutree.npz > $example1_dir/relationship_reconstruction_loss.txt

# compute the log perplexity, and save it to a text file
python3 $ORCH_DIR/metrics/perplexity.py truth=$example1_dir/example1.truth.neutree.npz orchard=$example1_dir/example1.orchard.neutree.npz  --ssm-fn $example1_dir/example1.ssm > $example1_dir/logperplexity.txt
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
- `-i`, `--num-instances` : int, default=None
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
- `-s`, `--seed` : int, default=None
    - seed value for duplicating orchard's results across multiple runs
- `-w`, `--node-order` : str, default="f_sum_node_order", choices=("f_sum_node_order", "randomized_node_order", "diverse_node_order")
    - Option for what order the nodes will be added to the tree in each parallel instance of Orchard
- `-x`, `--max-placements` : int, default=10
    - the maximum number of placements that can be tried for some node `u`, with respect to a single node `v`. This is not the same as the `branching_factor`. This parameter should not be changed as it's only applicable when the phylogeny may be a star tree
- `-z`, `--ignore-zero-probs` : bool, default=False
    - a flag to ignore node placements of `u` that have a zero probability. Setting this flag may increase runtimes, however, for most reconstruction problems this will only 
    provide a small decrease in runtime


Running Orchard
==================
**What are the most important parameters for Orchard?**

- The *beam width* (`-k`, `--beam-width`) is used to control the size of the priority queue that Orchard stores partially constructed trees in. The default is `-k=1`, which means Orchard will sample only one tree for each parallel instance it is run. The larger this priority queue, the more complete trees that Orchard will construct. This parameter helps improve reconstructions by exploring more tree structures. The major drawback is that the run time of Orchard increases linearly with the size of the beam width. With `-k=1`, we'll return a single tree from each parallel instance of Orchard.

- The *branching factor* (`-f`, `--branching-factor`) is used to control the number of places Orchard considers putting a new node in an existing partial tree. The default is `-f=20`, which means Orchard will consider 20 possible placements when adding a new node to a partial tree. Every placement of a new node in a partial tree is scored using a fast approximation, and then the *f* placements with the highest scores are formally used to construct trees. A maximum of *k* (i.e., the beam width) of these *f* partial trees can be kept according to Orchard's sampling routine. Therefore, it's very important to choose a large enough value of *f* such that a sufficient number of tree structures are considered during search. At the same time, constructing *f* trees involves fitting a cellular prevalence matrix to each tree which can be computationally expensive with large values of *f*. Generally, the run time of Orchard increases linearly with the size of the branching factor.

- The *expansion factor* (`'-e`, `--expansion-factor`) is used to control how many partial trees are extracted from the priority queue at each iteration. The default is `-e=1`, which results in a *best first search*. This is an efficient style of search that leads us to accurate tree structures quickly. An alternative approach is to set `-e=k`, where *k* is the beam width. When `-e=k`, Orchard will perform a *breadth first search* that will explore extensions of all of the partial trees in its priority at each iteration. Increasing the expansion factor can result in better reconstructions, however, run times will increase especially for large reconstruction problem and when using a large beam width.

- The *number of parallel instances of Orchard* (`-i`, `--num-instances` is used to control how many instances of the Orchard sampling routine to run. Each parallel instance of Orchard starts with its own node ordering, empty priority queue, etc. A related parameter is the *number of CPU core* to provide to Orchard (`-n`, `--num-cpu-cores`), which controls how many CPU cores Orchard is provided with. The default number of CPU cores provided to Orchard is the number of CPU cores on the machine it's being run. If `--num-instances` is not defined, then by default the number of parallel instances of Orchard run is the number of CPU cores provided to Orchard. In general, we recommend running no more than one instance of the Orchard sampling routine per CPU core that is provided to it. We recommend only using the `--num-cpu-cores` parameter to change how many parallel instances of Orchard to run. Running multiple instances of Orchard allows Orchard to explore more tree structures, and allows for the use of different initialization schemes.

- The *node order* (`-w`, `--node-order`) is used to control the order in which nodes are added to the tree during Orchard's sampling routine. We recommend using the *F sum node order* (`-w=fsum`) which orders the nodes in descending order according to the sum of the nodes data-implied cellular prevalences across all samples. If multiple CPU cores are provided to Orchard, it can also be interesting to use a *diverse node ordering* (`-w=diverse`) which initializes one instance of Orchard's sampling routine with the *F sum node order* and then randomizes the node order for all remaining parallel instances of Orchard.

**What is the estimate run time of Orchard?**

There are many possible parameter setups for Orchard. We provide here a table of parameter ranges with their estimates run times. Please not that there are other parameters that can modify the expected run times.


| # clones/mutations  | # samples| Beam width (k) | Branching Factor | Estimated Run Time |
|---------------------| ---------|----------------|------------------|--------------------|
| 1-29                | 1-100    | 1-50           | 20-100 		     | 0-2 minutes        |
| 30-49	              | 1-100    | 1-50           | 20-100           | 1-100 minutes      |
| 50-99               | 1-100    | 1-50           | 20-100  		 | 2-300 minutes      |
| 100-200 			  | 1-100    | 1-50           | 20-100           | 2-600 minutes      | 
| 200-1000            | 1-100    | 1-5            | 20-100           | 2-1800 minutes     |

Although some of these run times ranges may be large (e.g., for 100-200 clones/mutations), Orchard generally finds very accurate trees with parameter setups that fall towards the lower end of the run time ranges. **PLEASE NOTE** that generally Orchard's run time scales linearly with the size of the beam width (*-k*), branching factor (*-f*) and expansion factor (*-e*).

Recommended parameter setups for Orchard
-------------------------


**Datasets with < 50 clones/mutations**

Run time will not be an issue. 

**DO**

- Use larger beam widths (*-k=10+*), larger branching factors (*-f=100*), multiple parallel search instances (*-n=8*)
- Consider using a *diverse search* (*-w=diverse*), *breadth first search* (set *-e=k* where *k* is the beam width)

**DO NOT**
- Limit the search by ignoring low probability node placements (*-z*)

**Example command**
```
python3 $ORCH_DIR/bin/orchard $example1_dir/example1.ssm \
		$example1_dir/example1.params.json \
        $example1_dir/example1.orchard.npz \
        -k 20 -f 100 -n 8
```

**Datasets with 50-100 clones/mutations**

Run times generally should not be an issue. If the data has closer to 100 clones/mutations, expect Orchard to run for 10-15 minutes with the example command provided below.  

**DO**

- Use beam widths up to 20 (*k <= 20*), branching factors up to 100 (-f <= 100*), multiple parallel search instances 
- Consider using a *diverse search* (*-w=diverse*), *breadth first search* should only be used for relatively small choices of *k* (set *-e=k* where *k* is the beam width)

**DO NOT**
- Limit the search by ignoring low probability node placements (*-z*)

**Example command**
```
python3 $ORCH_DIR/bin/orchard $example1_dir/example1.ssm \
		$example1_dir/example1.params.json \
        $example1_dir/example1.orchard.npz \
        -k 10 -f 100 -n 8
```

**Datasets with 100-200 clones/mutations**

Run times can become quite large, especially with beam widths greater than 10 (*k > 10*).

**DO**

- Use beam widths up to 10 (*k <= 10*), branching factors up to 100 (-f <= 100*), multiple parallel search instances 
- Consider using a *diverse search* (*-w=diverse*), *breadth first search* should only be used for relatively small choices of *k* (set *-e=k* where *k* is the beam width)
- Consider limiting the search by ignoring low probability placements  by adding the flag `-z`

**DO NOT**
- Expect Orchard to complete in a timely manner if using large beam widths or large expansion factors

**Example command**
```
python3 $ORCH_DIR/bin/orchard $example1_dir/example1.ssm \
		$example1_dir/example1.params.json \
        $example1_dir/example1.orchard.npz \
        -k 10 -f 100 -n 8 -z
```

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

Here are a few high level things to note about the SSM file format:

- The SSM file is a *tab-separated* file
- The `id` field must be unique for each entry, and an `id` must follow the following regular expression format `s\d+`
- The `name` field can be any string, and does not need to be unique
- The `var_reads`, `total_reads`, and `var_read_prob` fields must be comma separated strings, and must all have the same length
- The values in the comma separated string for `var_reads`, `total_reads`, `var_read_prob` are ordered such that each index in the comma separated string is associated with a specific sample, and the values must appear in the same order that the `sample` names appear in the parameters file.
- The values in the comma separated string for `var_reads` and `total_reads` must be non-negative integers
- The values in the comma separated string for `var_read_prob` must be floats in the range (0,1]


Parameters File
-----------------

The parameters file format has also been thoroughly detailed on the [omicsdata readthedocs](https://omicsdata.readthedocs.io/en/latest/brief.html#parameter-file). This file format has also been detailed in the [Pairtree STAR protocol](https://www.sciencedirect.com/science/article/pii/S266616672200586X#p0100) and the [Pairtree GitHub](https://github.com/morrislab/pairtree#input-files).

Here are the most relevant key/value pairs in the `example1` parameters file 
```
{"samples": ["Sample 1", "Sample 2", "Sample 3"], "clusters": [["s9"], ["s0"], ["s1"], ["s2"], ["s4"], ["s3"], ["s7"], ["s5"], ["s8"], ["s6"]]}
```

Here are a few high level things to note about the parameters file format:

- The parameters file is a structured object format similar to [Javascript Object Notation](https://en.wikipedia.org/wiki/JSON)
- The values in the `samples` list must be unique strings
- The values in the `clusters` list must be lists, where each sublist contains one or more `id` values defined in the corresponding SSM file

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

Phylogeny aware clustering
===========================
The *phylogeny-aware clustering* method is an agglomerative clustering algorithm that infers a clone tree from a clone/mutation tree. It can also be used to find plausible mutation clusters after Orchard has been used to construct a clone/mutation tree on a dataset. 

**When should I use phylogeny-aware clustering?**

1. You do not think other mutation clustering algorithms are constructing realistic clusters
2. You have a priori knowledge about the number of mutation clusters in a cancer
3. You want to analyze a cancer's data using clones trees with varying numbers of clusters

**What are the most important parameters for phylogeny-aware clustering?**

The phylogeny-aware clustering algorithm joins pairs of nodes based on the distance between the cellular prevalence of the mutations in each node. As a result, we can change the criterion for the distance between in each node (`-l`, `--linkage-criterion`), and we can change the actual distance metric that is used (`-m`, `--metric`).

- The *linkage criterion*  (`-l`, `--linkage-criterion`) is used to select how the type of distance to compute between pairs of nodes (u,v). By default, the linkage criterion is the minimum distance between all pairs of mutations (i,j) (`-l=min_linkage`) where i is a mutation in node u, and j is a mutation in node v.

- The *distance metric* (`-m`, `--metric`) is used to compute the distance between all pairs of mutation (i,j). The default distance metric is the *euclidean distance* (`-m=euclidean`). Any metric that can be used with the [sklearn.metrics.pairwise_distances](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances.html) function can be use.

**What are the outputs of the phylogeny-aware clustering algorithm?**

The phylogeny-aware clustering algorithm outputs zipped archive using [omicsdata.npz.archive.Archive](https://omicsdata.readthedocs.io/en/latest/omicsdata.npz.html#omicsdata.npz.archive.Archive). This archive contains all of the clones trees of size *2,...,n+1* and their related data. By default, the name of the zipped archive will be *cluster.npz*.

FAQ
=======

Here are some frequently asked questions/errors:

### `AssertionError: Could not find projectppm library.`

The cause of this error is that the projectppm repository could not be found in its expected path. Please ensure that the [installation directions](#installation) have been followed properly.

### `AssertionError: var_reads for s0 is of length 1, but 3 expected'`

The cause of this error is that mutation `s0` listed in the SSM file provided to Orchard does not have variant read counts listed for all 3 samples that were defined in the corresponding parameters file provided to Orchard. Orchard expects that there are values provided in `var_reads`, `total_reads`, and `var_read_prob` for all of the samples listed in the parameters file. If for whatever reason you have data missing for a sample, please see the [Pairtree STAR Method section on imputing read count data](https://www.sciencedirect.com/science/article/pii/S266616672200586X#sectitle0055).
