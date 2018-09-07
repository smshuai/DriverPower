Tutorial with an working example
===================================

DriverPower aims to identify cancer driver candidates from **somatic** variants of a tumour cohort. The minimal number of
samples required in the PCAWG study is 15. However, more samples will have more power to detect rare driver events.
Both whole-genome sequencing and panel sequencing variants can be used.

The basic test unit accepted by DriverPower is called a **genomic element**, which is a set of disjoint (or continuous)
genomic regions (see the figure below). In our manuscript, we also removed all blacklist regions to reduce the effect of
variant artifacts.
In the real world, a genomic element can be all exons of a gene, or all TF binding sites among a promoter etc.

.. image:: /pics/genomic_element.png

Typically 1K to 100K genomic elements are tested together such as ~20K genes in the human genome.
For each element, given its functional impact, mutational burden as well as thousands of features, DriverPower will
report a p-value and q-value associated with that element.


0: Download example data
------------------------

Our example data are hosted on `figshare
<https://figshare.com/projects/DriverPower_Dataset/36065>`_, on which you can find the following files:

1. random_mutations.tsv.gz [14 MB]
2. train_feature.hdf5.part1 [4 GB]
3. train_feature.hdf5.part2 [4 GB]
4. train_feature.hdf5.part3 [900 MB]
5. test_feature.hdf5 [1.4 GB]
6. train_elements.tsv.gz [7 MB]
7. test_elements.tsv []
8. whitelist.bed.gz []

.. important:: You can run DriverPower for your own data by simply replacing
    **random_mutations.tsv** with your mutations.

The training features has been divided into three parts,
so you will need to merge them:

.. code-block:: bash

    cat train_feature.hdf5.part1 train_feature.hdf5.part2 train_feature.hdf5.part3 > train_feature.hdf5
    md5sum train_feature.hdf5  # should be cb4af6bc7979efa087955f94909dd065

1: Install Python and DriverPower
--------------------------------------

For this tutorial, we used a brand new virtual machine with 15 CPUs and 125GB RAM.
The OS we used is Ubuntu 18.04.1 LTS (GNU/Linux 4.15.0-33-generic x86_64).

To install Python3 using Anaconda3 if you don't already have:

.. code-block:: bash

    # This is the most recent version as of 2018-09-07
    wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh
    # You may need sudo to run this
    ./Anaconda3-5.2.0-Linux-x86_64.sh
    . .bashrc
    # [Optional] Create a new environment for DriverPower
    conda create --name driverpower
    # [Optional] Activate the environment
    source activate driverpower

Then we can install required packages and DriverPower:

.. code-block:: bash

    conda install -c conda-forge xgboost
    conda install -c bioconda pybedtools
    pip install driverpower


2: Make response tables
-----------------------

The response (y; dependent variable) table records the observed number of mutations, number of mutated samples and the length per genomic element.
This table is required for both ``model`` (training) and ``infer`` (test) sub-commands.

You can make them easily using our helper script ``script/prepare.py``. Inputs will be mutations and elements:

.. code-block:: bash

    # Training responses
    python script/prepare.py random_mutations.tsv train_elements.tsv whitelist.bed train_y.tsv
    # Test responses
    python script/prepare.py random_mutations.tsv test_elements.tsv whitelist.bed test_y.tsv


2: Build the background mutation rate model
-------------------------------------------
The background mutation rate (BMR) model is used to estimate the expected number of somatic mutations for each genomic element,
given its features. DriverPower sub-command ``model`` is used to train BMR models. To build the BMR model, training features
(X) and responses (y) are required. DriverPower supports two algorithms for the BMR model, generalized linear models (GLM)
and gradient boosting machines (GBM).

Here we show how to build a GBM with our example data:

.. code-block:: bash

    driverpower model \
        --feature train_features.hdf5 \
        --response train_y.tsv \
        --method GBM

Step 2: Build the background mutation rate model
------------------------------------------------

See :ref:`the model sub-command <model>` for all parameters and notes. Example code snippets are as follows:

.. code-block:: bash
    :caption: *Train a generalized linear model with feature selection by randomized lasso*
    :name: train-glm

    $ driverpower model \
        --feature PATH_TO_X \
        --response PATH_TO_y \
        --method GLM

.. code-block:: bash
    :caption: *Train a gradient boosting machine with default parameters*
    :name: train-gbm

    $ driverpower model \
        --feature PATH_TO_X \
        --response PATH_TO_y \
        --method GBM




Step 3: Call the driver candidates
----------------------------------

.. code-block:: bash
    :caption: *Call driver candidates with CADD scores*

    $ driverpower infer \
        --feature PATH_TO_X \
        --response PATH_TO_y \
        --modelInfo PATH_TO_model_info \
        --funcScore PATH_TO_func_score \
        --funcScoreCut 'CADD:0.01' \
        --name 'DriverPower' \
        --outDir ./output/