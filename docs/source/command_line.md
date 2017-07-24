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

The `model` sub-command is used to build the background mutation model.
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

## The `infer` sub-command

The `infer` sub-command is used to call drivers from test data.
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
