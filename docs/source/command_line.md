# Command-line interface

DriverPower uses a command-line interface.
To view available sub-commands and usage of DriverPower, type
```bash
$ driverpower -h
usage: driverpower [-h] [-v] {model,infer} ...

DriverPower v1.0.0: Combined burden and functional impact tests for coding
and non-coding cancer driver discovery

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  Print the version of DriverPower

The DriverPower sub-commands include:
  {model,infer}
    model        Build the background mutation model
    infer        Infer driver elements
```

## The `model` sub-command

The `model` sub-command is used to build the background mutation model with training data.
Command-line options for this command can be viewed as follows:
```bash
$ driverpower model -h
usage: driverpower model [-h] --feature str --response str [--featImp str]
                         --method {GLM,GBM} [--featImpCut float] [--name str]
                         [--modelDir str]

optional arguments:
  -h, --help          show this help message and exit

Input data:
  --feature str       Path to the training feature table (default: None)
  --response str      Path to the training response table (default: None)
  --featImp str       Path to the feature importance table [optional]
                      (default: None)

Parameters:
  --method {GLM,GBM}  Algorithms to use (default: None)
  --featImpCut float  Cutoff of feature importance score [optional] (default:
                      0.5)
  --name str          Identifier for output files [optional] (default:
                      DriverPower)
  --modelDir str      Directory of output model files [optional] (default:
                      ./output/)
```

* **Input**
    * Required data are `--feature` and `--response`.
    * Optional data are `--featImp`.
* **Output**
    * GLM model is saved at `"[modelDir]/[name].GLM.model.pkl"`.
    The corresponding model information is saved at `"[modelDir]/[name].GLM.model_info.pkl"`.
    * GBM models are saved as`"[modelDir]/[name].GBM.model.fold[k]"`.
    The corresponding model information is saved at `"[modelDir]/[name].GBM.model_info.pkl`
    * Feature importance table: returned for both GLM and GBM when no input `--featImp`.
    For GLM, feature importance is the number of times a feature is used by randomized lasso.
    For GBM, feature importance is the average gain of the feature across all gradient boosting trees.
```eval_rst
.. important:: Please DO NOT rename the model files because their names are recorded in model information. Model files can be moved to another directory as long as ``--modelDir`` is specified in ``infer``.
```
* **Parameters**
    * `--method`: [*required*] method used for the background model. Must be GLM or GBM.
    * `--featImpCut`: [*optional*] cutoff for feature importance score. Features with score >= cutoff will be used in the model. Only used for GLM. Default is 0.5.
    * `--gbmParam`: [*optional*] path to the parameter pickle for GBM. The pickle file must contain a valid python dictionary for [XGBoost parameters](https://github.com/dmlc/xgboost/blob/master/doc/parameter.md).
    * `--gbmFold`: [*optional*] Number of model fold to train for GBM. Fold must be an integer >= 2. Default value is 3.
    * `--name`:  [*optional*] Prefix for output files. Default is 'DriverPower'.
    * `--modelDir`: [*optional*] Directory for output model and model information files. Default is './output/'.
* **Notes**
    * Input data requirements can be found at [data format](https://driverpower.readthedocs.org/en/latest/data_format.html).
    * `--gbmFold` controls the split of training data by k-fold cross-validation.
    For example, default 3-fold model means the training data are divided into 3 equal folds.
    Each time, two folds of data are used to train the model and
    the rest fold is used to validate the model.
    Hence, three model files will be generated at the end. The prediction will then be the average of three models.  
    * Training phase can take hours for large training set and consume a large amount of memories.
    For example, using our default training set (~1M elements and 1373 features), training randomized lasso plus GLM takes about xx hours and xx RAMs.
    Training 3-fold GBM with xx cores takes about xx hours and xx RAMs. 
    * The default pickle file for `--gbmParam` is generated as follow:

```python
import pickle
# Default parameters for XGBoost used in DriverPower
param = {'max_depth': 8,
         'eta': 0.05,
         'subsample': 0.6,
         'nthread': 15,
         'objective': 'count:poisson',
         'verbose': 0,
         'max_delta_step': 1.2,
         'eval_metric': 'poisson-nloglik'}
# Dump to pickle file ./xgb.param.pkl
with open('./xgb.param.pkl', 'wb') as f:
    pickle.dump(param, f)
    
``` 

## The `infer` sub-command

The `infer` sub-command is used to call drivers from test data with pre-trained models.
Command-line options for this command can be viewed as follows:
```bash
$ driverpower infer -h
usage: driverpower infer [-h] --feature str --response str --modelInfo str
                         [--modelDir str] [--funcScore str]
                         [--method {auto,binomial,negative_binomial}]
                         [--scale float] [--funcScoreCut str] [--geoMean bool]
                         [--name str] [--outDir str]

optional arguments:
  -h, --help            show this help message and exit

Required input data:
  --feature str         Path to the test feature table (default: None)
  --response str        Path to the test response table (default: None)
  --modelInfo str       Path to the model information (default: None)

Optional input data:
  --modelDir str        Directory of the trained model(s) (default: None)
  --funcScore str       Path to the functional score table (default: None)

Parameters:
  --method {auto,binomial,negative_binomial}
                        Test method to use [optional] (default: auto)
  --scale float         Scaling factor for theta in negative binomial
                        distribution [optional] (default: 1)
  --funcScoreCut str    Score name:cutoff pairs for all scores e.g.,
                        "CADD:0.01;DANN:0.03;EIGEN:0.3" [optional] (default:
                        None)
  --geoMean bool        Use geometric mean in test [optional] (default: True)
  --name str            Identifier for output files [optional] (default:
                        DriverPower)
  --outDir str          Directory of output files [optional] (default:
                        ./output/)
```

* **Input**
* **Output**
* **Parameters**
    * `--method`: [optional] probability distribution used to generate p-values (binom, negbinom, auto). Default is auto. Decision will be made automatically based on the dispersion test.
    * `--featImpCut`: [optional] same as in sub-command `model`.
    * `--scale`: [optional]
    * `--funcScoreCut`: [optional]
    * `--geoMean`: [optional] Use geometric mean of nMut and nSample in test. Default is `True`.
    * `--name`: [optional] Prefix for output files. Default is 'DriverPower'.
    * `--outDir`: [optional] Directory for output files. Default is './output/'.
* **Notes**