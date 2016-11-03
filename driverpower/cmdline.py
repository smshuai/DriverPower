import argparse
import os
import logging
import sys
import pandas as pd
from driverpower.load import load_memsave
from driverpower.preprocess import get_response, scaling
# from driverpower.feature_select import fselect
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
        help='String (default: data.h5). Output file name (HDF5 format)')
    #
    #
    #
    parser_select = subparsers.add_parser('select', help='Feature selection module')
    #
    #
    #
    parser_model = subparsers.add_parser('model', help='Modelling module')

    args = parser.parse_args()
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

def main():
    args = get_args()
    if args.subparser_name == 'preprocess':
        run_preprocess(args)

if __name__ == '__main__':
    main()
