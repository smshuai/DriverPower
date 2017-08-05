Quick start
===========

Step 0: Define the task
-----------------------
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

Step 1: Prepare the data
------------------------
DriverPower runs on processed training and test data as described in :doc:`data format </data_format>`.

Step 2: Build the background mutation rate model
------------------------------------------------
The background mutation rate (BMR) model is used to estimate the expected number of somatic mutations for each genomic element,
given its features. DriverPower sub-command ``model`` is used to train BMR models. To build the BMR model, training features
(X) and responses (y) are required. DriverPower supports two algorithms for the BMR model, generalized linear models (GLM)
and gradient boosting machines (GBM).

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