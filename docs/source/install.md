# Installation guide

DriverPower requires Python >= 3.5 and some other computing packages.
If you don't have Python 3.5 or higher, we recommend to install Python with [Anaconda](https://www.continuum.io/downloads).

## Install from conda
The best way to install driverpower is to use `conda`:
```console
# Setup bioconda channels (see https://bioconda.github.io/user/install.html#set-up-channels)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Create a new env and install DriverPower
conda create -n driverpower
conda activate driverpower
conda install -c smshuai driverpower
```

## Install from pypi or source
### Install from pypi
You can install DriverPower from the
[Python Package Index (PyPI)](https://pypi.python.org/pypi/DriverPower/)
by simply typing the following command in your terminal:
```console
$ pip install driverpower
```

### Install from source
You can also download the [latest source](https://github.com/smshuai/DriverPower/releases/latest/) from GitHub to install.
For example, to install version 1.0.x (change x to the right version number):
```console
$ tar -xzf DriverPower-1.0.x.tar.gz
$ cd DriverPower-1.0.x/ && pip install .
```
