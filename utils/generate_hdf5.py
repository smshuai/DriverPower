"""
This script is used to convert a tsv feature table to HDF format for faster loading.
"""
import pandas as pd
from driverpower.dataIO import read_feature


X = read_feature(path='./train.feature.tsv.gz') # path to the feature table
store = pd.HDFStore('./train.feature.hdf5') # create output HDF5
store['X'] = X
store.close()
