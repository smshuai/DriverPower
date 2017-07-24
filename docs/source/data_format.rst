Data format
===========

Response table (``--response``)
-------------------------------
The response (y; dependent variable) table records the observed number of mutations, number of mutated samples and the length per genomic element. This table is required for both ``model`` (training) and ``infer`` (test) sub-commands.

1. **Format**: TSV (tab-delimited text) with header. Compressed tables are also accepted.
2. **Fields**:

    * ``binID``: identifier of genomic element and used as key.
    * ``length``: effective length of the element.
    * ``nMut``: number of observed mutations.
    * ``nSample``: number of observed samples with mutations.
    * ``N``: total number of samples.
3. **Example**:

========    ======  ====    ======= ===
binID       length  nMut    nSample N
========    ======  ====    ======= ===
TP53.CDS
KRAS.CDS
...         ...     ...     ...     ...
========    ======  ====    ======= ===

Feature table (``--feature``)
-----------------------------
The feature (X; independent variable) table records genomic features per element.
This table can be large (~10GB) for thousands of features and ~1M genomic elements.
This table is required for both ``model`` (training) and ``infer`` (test) sub-commands.

1. **Format**:

    * TSV (or compressed TSV) with header. Loading may be slow for large datasets.
    * `HDF5 <https://pandas.pydata.org/pandas-docs/stable/io.html#io-hdf5>`_ (\*.h5 or \*.hdf5). The HDF5 must contain key ``X``, which is the feature table in `pandas.DataFrame <https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html>`_. Used for fast loading.

2. **Fields**:

    * ``binID``: identifier of genomic element and used as key.
    * Column 2+: one feature per column. Unique feature names are required.

3. **Example**:

========    ======= ============    ===
binID       GERP    E128-H3K27ac    ...
========    ======= ============    ===
TP53.CDS    4.80287 1.19475         ...
KRAS.CDS    3.56563 2.53435         ...
========    ======= ============    ===

Feature importance table (``--featImp``)
----------------------------------------
The feature importance table records the feature importance score derived from randomized lasso, gradient boosting machine or other algorithms.
This table can be used together with the ``--featCut`` argument to specify features kept in the background model.
This table is only used and optional for the ``model`` sub-command with method ``GLM``.
If not provided, DriverPower will generate and output it automatically.

1. **Format**: TSV (or compressed TSV) with header.
2. **Fields**:

    * ``name``: name of features, should match the column name in feature table.
    * ``importance``: feature importance score.

3. **Example**:

============    ==========
name            importance
============    ==========
GERP            0.3
E128-H3K27ac    0.5
...             ...
============    ==========

Functional score table (``--funcScore``)
----------------------------------------
The functional score table records the functional impact score per genomic element.
DriverPower can work with any functional impact score scheme.
This table is only used and optional for the ``infer`` sub-command.

1. **Format**: TSV (compressed TSV) with header.
2. **Fields**:

    * ``binID``: identifier of genomic element and used as key.
    * Column 2+: one type of functional impact score per column.

3. **Example**:

========    ======= =======  ===
binID       CADD    EIGEN    ...
========    ======= =======  ===
TP53.CDS    4.80287 1.19475  ...
KRAS.CDS    3.56563 2.53435  ...
========    ======= =======  ===

Model and model information files (output of ``model``)
-------------------------------------------------------

Result table (output of ``infer``)
----------------------------------
