# Installation guide

DriverPower requires Python >= 3.5 and some other computing packages.
If you don't have Python 3.5 or higher, we recommend to install Python with [Anaconda](https://www.continuum.io/downloads).

## Install from pypi
You can install DriverPower from the
[Python Package Index (PyPI)](https://pypi.python.org/pypi/DriverPower/)
by simply typing the following command in your terminal:
```bash
$ pip install driverpower
```

```eval_rst
.. important:: DriverPower requires `XGBoost <https://github.com/dmlc/xgboost>`_ for building gradient boosting machines. ``pip`` may not install XGBoost correctly especially for OS X owing to C++ compiler issue. See `XGBoost installation guides <https://pypi.python.org/pypi/xgboost/>`_ for more information.
```


## Install from source
You can also download the latest source from [GitHub](https://github.com/smshuai/DriverPower/releases) and install DriverPower like:
```bash
$ wget link_to_latest_release
$ cd DriverPower && pip install .
```
See the importance note above if you have issue with the XGBoost installation.