import argparse
import os
import logging
import sys
import pandas as pd
import numpy as np
from scipy.special import logit
from driverpower.load import load_memsave
from driverpower.preprocess import get_response, scaling, sampling
from driverpower.feature_select import run_lasso, run_rndlasso, run_spearmanr, run_fregression
from driverpower import __version__
# from driverpower.model import model
# from driverpower.func_adj import func_adj


# logging config
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
    format='%(asctime)s | %(name)s | %(levelname)s: %(message)s',
    datefmt='%m/%d/%Y %H:%M:%S')
# create logger
logger = logging.getLogger('DP')


def get_args():
    parser = argparse.ArgumentParser(prog='driverpower')
    # global argument
    parser.add_argument('-v', '--version', dest='version', action="store_true",
        help='Print the version of DriverPower')
    subparsers = parser.add_subparsers(title='The DriverPower sub-commands include',
        dest='subparser_name')
    #
    # Load and preprocess data
    #
    parser_preprocess = subparsers.add_parser('preprocess',
        help='Load and preprocess data')
    # required parameters
    re_parser = parser_preprocess.add_argument_group(title="required arguments")
    # required 3 tables for test data
    re_parser.add_argument('-c', '--count', dest='path_ct', required=True, type=str,
        help='Path to the count table')
    re_parser.add_argument('-f', '--feature', dest='path_cv', required=True, type=str,
        help='Path to the feature table')
    re_parser.add_argument('-l', '--length', dest='path_cg', required=True, type=str,
        help='Path to the effective length table')
    # Optional parameters
    op_parser = parser_preprocess.add_argument_group(title="optional parameters")
    op_parser.add_argument('--len_threshold', dest='len_threshold', type=int, default=500,
        help='Integer (default: 500). Bins with length < len_threshold will be discarded')
    op_parser.add_argument('--recur_threshold', dest='recur_threshold', type=int, default=2,
        help='Integer (default: 2). Bins having mutations in < recur_threshold samples will be discarded')
    op_parser.add_argument('-o', '--output', dest='out', type=str, default='data.h5',
        help='Path to the output file (default: ./data.h5)')
    #
    # Feature selection
    #
    parser_select = subparsers.add_parser('select',
        help='Run feature selection on preprocessed data')
    # required parameters
    re_select = parser_select.add_argument_group(title="required arguments")
    re_select.add_argument('-d', '--data', dest='path_data', required=True, type=str,
        help='Path to the preprocessed training set (HDF5)')
    # optinal parameters
    op_select = parser_select.add_argument_group(title="optional parameters")
    op_select.add_argument('--scaler', choices=['robust', 'standard', 'none'],
        type=str, dest='scaler', default='robust',
        help='robust or standard (default: robust). Scaler used to scale the data')
    op_select.add_argument('--sampling',
        type=float, dest='sampling', default=1.0,
        help='Number > 0 (default: 1). Sampling the data based on the provided value. Value in (0,1] is used as a fraction. Value > 1 is used as the number of data points.')
    op_select.add_argument('-o', '--output', dest='out', type=str, default='feature_select.tsv',
        help='Path to the output file (default: ./feature_select.tsv)')
    #
    #
    #
    parser_model = subparsers.add_parser('model', help='Modelling module')

    args = parser.parse_args()
    # no argument, print main help instead
    if len(sys.argv)==1:
        parser.print_help()
        parser.exit(1)
    if args.version:
        print("DriverPower", __version__)
    return args

def run_preprocess(args):
    # initial output HDF5
    store = pd.HDFStore(args.out, mode='w')
    ct, cg, cv = load_memsave(args.path_ct,
        args.path_cg, args.path_cv,
        args.len_threshold, args.recur_threshold)
    # get response
    ybinom = get_response(ct, cg)
    # y to pd.DF
    ybinom = pd.DataFrame(ybinom, columns=['ct','len_ct'], index=cg.index)
    # write to store
    store.append('X', cv, chunksize=50000)
    store['y'] = ybinom
    store.close()
    logger.info('Data pre-process done!')

def run_select(args):
    # load data from HDF5
    X = pd.read_hdf(args.path_data, 'X')
    y = pd.read_hdf(args.path_data, 'y')
    # check index (binID)
    assert np.array_equal(X.index, y.index), 'X and y have different row indexes'
    logger.info('Successfully load X with shape: {}'.format(X.shape))
    logger.info('Successfully load y with shape: {}'.format(y.shape))
    # Sampling data
    X, y = sampling(X, y, args.sampling)
    # y to np.array
    y = y.as_matrix()
    # feature names
    fnames = X.columns.values 
    # scale for Xtrain only
    X = scaling(Xtrain=X.as_matrix(), scaler_type=args.scaler)
    # spearmanr
    rho = run_spearmanr(X, y)
    # f regression
    freg = run_fregression(X, y)
    # run LassoCV
    lasso = run_lasso(X, y)
    # run rndlasso
    rndlasso = run_rndlasso(X, y, lasso.alpha_)
    # results
    res = pd.DataFrame(np.array([rho, freg, np.abs(lasso.coef_), rndlasso.scores_]).T,
        index=fnames,
        columns=['rho', 'freg','lasso', 'rndlasso'])
    res.index.name = 'fname'
    res.to_csv(args.out, sep='\t')
    logger.info('Feature selection done!')

def main():
    args = get_args()
    if args.subparser_name == 'preprocess':
        run_preprocess(args)
    elif args.subparser_name == 'select':
        run_select(args)


if __name__ == '__main__':
    main()
