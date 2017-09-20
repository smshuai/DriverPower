""" Convert a tsv feature table to hdf5 for faster loading
$ tsv2hdf5.py FEATURE OUT_PATH [USE_BINS]
"""

import pandas as pd
import numpy as np
import sys
from driverpower.dataIO import read_feature


if __name__ == '__main__':
    X_path = sys.argv[1]
    out_path = sys.argv[2]
    if len(sys.argv) == 4:
        use_bins = np.loadtxt(sys.argv[3], dtype=np.str_)
    else:
        use_bins = None
    X = read_feature(X_path)
    if use_bins is not None:
        X = X.loc[use_bins]
    store = pd.HDFStore(out_path)
    store['X'] = X
    store.close()
