# DriverPower
![Github version](https://img.shields.io/badge/version-1.0.0-yellow.svg)
[![GitHub license](https://img.shields.io/badge/license-AGPL-blue.svg)](./LICENSE)
[![PyPI version](https://badge.fury.io/py/DriverPower.svg)](https://badge.fury.io/py/DriverPower)
[![Documentation Status](https://readthedocs.org/projects/driverpower/badge/?version=latest)](http://driverpower.readthedocs.io/en/latest/?badge=latest)


## Introduction

DriverPower is a tool used to discover potential coding and non-coding cancer driver elements from tumour whole-genome or whole-exome somatic mutations.
## Installation

DriverPower requires Python >= 3.5 and some other packages. If you don't have Python 3.5 or higher, we recommend to install Python with [Anaconda](https://www.continuum.io/downloads).

To install DriverPower, you can either download and install it from [Python Package Index (PyPI)](https://pypi.python.org/pypi/DriverPower/) by typing the following command in your terminal
```bash
$ pip install driverpower
```
or download the latest source to install:
```bash
$ git clone https://github.com/smshuai/DriverPower.git
$ cd DriverPower && pip install .
```

To view available sub-commands and help of DriverPower, type
```
$ driverpower -h
usage: driverpower [-h] [-v] {preprocess,select,model,detect} ...

DriverPower v0.5.0dev: Combined burden and functional impact test for coding
and noncoding cancer elements

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print the version of DriverPower

The DriverPower sub-commands include:
  {preprocess,select,model,detect}
    preprocess          Load and preprocess data
    select              Run feature selection on preprocessed data
    model               Find driver bins with preprocessed training and test
                        data (deprecated)
    detect              Detect cancer driver elements
```

To view help for a particular sub-commands like `select`, type
```
$ driverpower select -h
usage: driverpower select [-h] -d PATH_DATA [--scaler {robust,standard,none}]
                          [--sampling SAMPLING] [-o OUT]

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -d PATH_DATA, --data PATH_DATA
                        Path to the preprocessed training set (HDF5)

optional parameters:
  --scaler {robust,standard,none}
                        robust or standard (default: robust). Scaler used to
                        scale the data
  --sampling SAMPLING   Number > 0 (default: 1). Sampling the data based on
                        the provided value. Value in (0,1] is used as a
                        fraction. Value > 1 is used as the number of data
                        points.
  -o OUT, --output OUT  Path to the output file (default:
                        ./feature_select.tsv)
```

## A quick example

Data, script and expected output for this example can be found under the [example](https://github.com/smshuai/DriverPower/tree/master/example) folder. Input data requirements can be found after this example.

> **Task**: Identify potential protein-coding driver genes and non-coding driver pomoters from cancer somatic mutations.

In the [./example/data/](https://github.com/smshuai/DriverPower/tree/master/example/data) folder we have:

1. **example.mut.gz**: somatic mutations. Here we use CNS-Medullo data as an example.
2. **example.cds.features.gz** / **example.promCore.features.gz**: reference features for CDS and promCore elements. Only 38 important features are kept for file size issue.
3. **example.train.h5**: pre-processed training data generated with `driverpower preprocess`. Only 38 important features are kept.

In the [./example/annot/](https://github.com/smshuai/DriverPower/tree/master/example/annot) folder we have:

1. **callable.bed.gz**: genome whitelist regions or so called "callable" regions.
2. **cds.bed.gz** / **promCore.bed.gz**: genome coordinates of CDS and promCore elements.

In order to detect drivers for multiple element set simultaneously, list of data and annotation files should be written into a single file like the [./example/test_list.tsv](https://github.com/smshuai/DriverPower/blob/master/example/test_list.tsv):

| name     | element                 | feature                             | func_names           |
|----------|-------------------------|-------------------------------------|----------------------|
| CDS      | ./annot/cds.bed.gz      | ./data/example.cds.features.gz      | cadd,eigen_coding    |
| promoter | ./annot/promCore.bed.gz | ./data/example.promCore.features.gz | cadd,eigen_noncoding |

`func_names` here specifies types of functional impact scores to use.

Finally, we can run DriverPower with the bash script [./example/run_example.sh](https://github.com/smshuai/DriverPower/blob/master/example/run_example.sh):
```bash
#!/usr/bin/env bash
driverpower detect \
    --variant    "./data/example.mut.gz" \
    --trainH5    "./data/example.train.h5" \
    --testFile   "./test_list.tsv" \
    --callable   "./annot/callable.bed.gz" \
    --cohortName "example.CNS-Medullo" \
    --outDir     "./output/"
```
Running logs look like:
```
$ cd ./example/ && ./run_example.sh
03/13/2017 14:28:32 | INFO: DriverPower 0.5.0dev
03/13/2017 14:28:32 | INFO: Sub-command - Detect
03/13/2017 14:28:34 | INFO: Successfully load 223792 mutations from 141 donors
03/13/2017 14:28:43 | INFO: 219822 (98.23%) mutations are in callable regions
03/13/2017 14:28:43 | INFO: Load training data...
03/13/2017 14:28:43 | INFO: Successfully load data for 141 samples
03/13/2017 14:28:43 | INFO: Successfully load X with shape: (90190, 38)
03/13/2017 14:28:43 | INFO: Successfully load y with shape: (90190, 2)
03/13/2017 14:28:43 | INFO: Use robust scaler
===== Test Set [1/2]: CDS (n=20185) =====
03/13/2017 14:28:48 | INFO: Number of mutations in set: 2818
03/13/2017 14:28:48 | INFO: Successfully load 38 features for 20185 bins
03/13/2017 14:29:01 | WARNING: 915 (32.47%) NA values found in EIGEN_CODING scores. All NAs are ignored
03/13/2017 14:29:01 | WARNING: 5 (0.18%) NA values found in CADD scores. All NAs are ignored
03/13/2017 14:29:01 | INFO: Find 13 elements with q-value <=  0.1
===== Test Set [2/2]: promoter (n=20164) =====
03/13/2017 14:29:04 | INFO: Number of mutations in set: 924
03/13/2017 14:29:05 | INFO: Successfully load 38 features for 20164 bins
03/13/2017 14:29:17 | WARNING: 4 (0.43%) NA values found in CADD scores. All NAs are ignored
03/13/2017 14:29:17 | WARNING: 148 (16.02%) NA values found in EIGEN_NONCODING scores. All NAs are ignored
03/13/2017 14:29:18 | INFO: Find 1 elements with q-value <=  0.1
```
Two output files could be found in ./example/output/, which should match files under ./example/expected_output/.


## Input data requirements
#### Response table (`--response`)
1. **Format**: TSV (tab-delimited text) with header. Compressed tables are also accepted.
2. **Fields**:
- `binID`: identifier of genomic element and used as key.
- `length`: effective length of the element.
- `nMut`: number of observed mutations.
- `nSample`: number of observed samples with mutations.
- `N`: total number of samples.
3. **Example**:
 
| binID | length | nMut | nSample | N |
|-------|--------|------|---------|---|
| TP53.CDS |  |  | 100 | |
| KRAS.CDS |  |  | 100 | |

-----
#### Feature table (`--feature`)

1. **Format**:
- TSV (compressed TSV) with header. Loading may be slow for large datasets.
- [HDF5](https://pandas.pydata.org/pandas-docs/stable/io.html#io-hdf5) (*.h5 or *.hdf5). The HDF5 must contain key `X`, which is the feature table in [pandas.DataFrame](https://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html). Used for fast loading.
2. **Fields**:
- `binID`: identifier of genomic element and used as key.
- Column 2+: one feature per column. Unique feature names are required.

3. **Example**:

| binID | GERP    | E128-H3K27ac | ... |
|-------|---------|--------------|-----|
| TP53.CDS  | 4.80287 | 1.19475      | ... |
| KRAS.CDS  | 3.56563 | 2.53435      | ... |

-----
#### Feature importance table (`--featImp`)
1. **Format**: TSV (compressed TSV with header.
2. **Fields**:
- `name`: name of features, should match the column name in feature table.
- `importance`: feature importance score.
3. **Example**:

| name | importance |
|------|-------------|
| GERP         | 0.3 |
| E128-H3K27ac | 0.5 |

-----
#### Functional score table (`--funcScore`)

1. **Format**: TSV (compressed TSV) with header.
2. **Fields**:
- `binID`: identifier of genomic element and used as key.
- Column 2+: one type of functional impact score per column.
3. **Example**:

| binID | CADD    | EIGEN | ... |
|-------|---------|--------------|-----|
| TP53.CDS  | 4.80287 | 1.19475      | ... |
| KRAS.CDS  | 3.56563 | 2.53435      | ... |

## Parameters
#### Sub-command `model`
- `--method`: [required] method used for the background model. Must be GLM or GBM.
- `--featImpCut`: [optional] cutoff for feature importance score. Features with score >= cutoff will be used in the model. Only used for GLM. Default is 0.5.
- `--gbmParam`: [optional] path to the parameter pickle for GBM. The pickle file must contain a valid python dictionary for [XGBoost parameters](https://github.com/dmlc/xgboost/blob/master/doc/parameter.md). The default pickle file is generated as follow:
```python
import pickle
# Default parameters for XGBoost used in DriverPower
param = {'max_depth':8,
         'eta':0.05,
         'subsample':0.6,
         'nthread':15,
         'objective':'count:poisson',
         'verbose':0,
         'max_delta_step':1.2,
         'eval_metric':'poisson-nloglik'}
# Dump to pickle file ./xgb.param.pkl
with open('./xgb.param.pkl', 'wb') as f:
    pickle.dump(param, f)
```
- `--name`: [optional] Prefix for output files. Default is 'DriverPower'.
- `--outDir`: [optional] Directory for output files. Default is './output/'.
-----
#### Sub-command `infer`
- `--method`: [optional] probability distribution used to generate p-values (binom, negbinom, auto). Default is auto. Decision will be made automatically based on the dispersion test.
- `--featImpCut`: [optional] same as in sub-command `model`.
- `--scale`: [optional]
- `--funcScoreCut`: [optional]
- `--geoMean`: [optional] Use geometric mean of nMut and nSample in test. Default is `True`.
- `--name`: [optional] Prefix for output files. Default is 'DriverPower'.
- `--outDir`: [optional] Directory for output files. Default is './output/'.

-----
## LICENSE
DriverPower is distributed under the terms of the [GNU General Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.txt).

## Change Log
- 2017/03/02: Release version 0.4.0, which is used in the PCAWG driver analysis.
