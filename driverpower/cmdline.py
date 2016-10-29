"""
"""
import argparse
import os
import logging
import sys
import pandas as pd
from driverpower.load import load_all
from driverpower.preprocess import preprocess
from driverpower.feature_select import fselect
from driverpower.model import model
from driverpower.func_adj import func_adj


# logging config
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
    format='%(asctime)s | %(name)s | %(levelname)s: %(message)s',
    datefmt='%m/%d/%Y %H:%M:%S')
# create logger
logger = logging.getLogger('DP')


def get_params():
    parser = argparse.ArgumentParser()

    re_parser = parser.add_argument_group(title="required arguments")
    # required 4 tables for test data
    re_parser.add_argument('-c', '--count', dest='path_ct_test', required=True, type=str,
        help='Path to the count table for test set')
    re_parser.add_argument('-f', '--feature', dest='path_cv_test', required=True, type=str,
        help='Path to the feature table for test set')
    re_parser.add_argument('-l', '--length', dest='path_cg_test', required=True, type=str,
        help='Path to the effective length table for test set')
    # required 3 tables for train data
    re_parser.add_argument('-C', '--Count', dest='path_ct_train', required=True, type=str,
        help='Path to the count table for training set')
    re_parser.add_argument('-F', '--Feature', dest='path_cv_train', required=True, type=str,
        help='Path to the feature table for training set')
    re_parser.add_argument('-L', '--Length', dest='path_cg_train', required=True, type=str,
        help='Path to the effective length table for training set')
    # Optional parameters
    op_parser = parser.add_argument_group(title="optional parameters")
    op_parser.add_argument('-o', '--outfolder', dest='out', type=str, default='output',
        help='String (default: output). Output folder name')
    op_parser.add_argument('-p', '--prefix', dest='pref', type=str, default='driverpower',
        help='String (default: driverpower). Output file prefix')
    op_parser.add_argument('--len_threshold', dest='len_threshold', type=int, default=500,
        help='Integer (default: 500). Bins with length < len_threshold will be discarded')
    op_parser.add_argument('--recur_threshold', dest='recur_threshold', type=int, default=2,
        help='Integer (default: 2). Bins having mutations in < recur_threshold samples will be discarded')
    op_parser.add_argument('--scaler', choices=['robust', 'standard', 'none'],
        type=str, dest='scaler', default='robust',
        help='robust or standard (default: robust). Scaler used to scale the data')
    op_parser.add_argument('--fselect', choices=['rndlasso', 'lasso', 'none'],
        type=str, dest='fselect', default='rndlasso',
        help='rndlasso | lasso | none (default: rndlasso). Feature selection method')
    op_parser.add_argument('--coding', dest='is_coding', action="store_true",
        help='Allow test for coding bins')
    op_parser.add_argument('--funcadj', dest='funcadj', type=int, default=85,
        help='Integer between 1 and 99 (default: 85). Strength of functional adjustment. Integer outside of (0, 100) will disable functional adjustment')
    op_parser.add_argument('-m', '--mutation', dest='path_mut', type=str, default=None,
        help='Path to the mutation table for test set')
    #
    args = parser.parse_args()

    # check all paths
    check_file(args.path_ct_test)
    check_file(args.path_cv_test)
    check_file(args.path_cg_test)
    check_file(args.path_ct_train)
    check_file(args.path_cv_train)
    check_file(args.path_cg_train)
    if 0 < args.funcadj < 100:
        # func adj ON
        logger.info('Functional adjustment is ON')
        check_file(args.path_mut)
    return args


def check_file(path):
    if not os.path.isfile(os.path.expanduser(path)):
        # no such file
        logger.error('File not found in {}'.format(path))
        sys.exit(1)


def main():
    args = get_params()
    if 0 < args.funcadj < 100:
        is_funcadj = True # turn on func adj
    else:
        is_funcadj = False
    if args.fselect == 'none':
        is_fselect = False
    else:
        is_fselect = True
    #
    # load all the data. fnames is feature names
    logger.info('Start data loading')
    (cg_test, ct_test, X_test, mut,
        cg_train, ct_train, X_train,
        fnames) = load_all(
        args.path_cg_test, args.path_ct_test, args.path_cv_test, args.path_mut,
        args.path_cg_train, args.path_ct_train, args.path_cv_train)
    #
    # get response, filter and scale.
    logger.info('Start data prepocessing')
    ## gnames is filtered binIDs in test data
    ## grecur is filtered recur in test data (index by gnames)
    (X_train, ybinom_train, X_test,
        ybinom_test, gnames, grecur) = preprocess(cg_test, ct_test, X_test, cg_train,
    ct_train, X_train, args.len_threshold, args.recur_threshold, args.scaler)
    #
    # feature selection
    if is_fselect:
        logger.info('Start feature selection')
        X_train, X_test, fscores = fselect(X_train, X_test, ybinom_train,
            fnames, method=args.fselect)
        fs_res = pd.DataFrame({'fname': fnames, 'fscore': fscores})
        fs_res.sort_values('fscore', ascending=False)
    else:
        logger.info('Skip feature selection')
        fs_res = None
    #
    # model
    logger.info('Start modelling')
    res = model(X_train, ybinom_train, X_test, ybinom_test, gnames, grecur, method='glm')
    #
    # func adj
    if is_funcadj:
        logger.info('Start functional adjustment')
        res = func_adj(res, mut, method='eigen',
            path_eigen='/u/sshuai/sshuai/func_score/eigen/v1.1',
            is_coding=args.is_coding, cutoff=85)
        res.sort_values('Pval', inplace=True)
    else:
        logger.info('Skip functional adjustment')
        res.sort_values('PvalM', inplace=True)
    #
    # Output
    out_dir = args.out
    if not os.path.exists(out_dir):
        # create output folder
        os.mkdir(out_dir)
    # main result
    res.to_csv(os.path.join(out_dir, args.pref+'.res.tsv'), sep='\t')
    # feature selection result
    if fs_res is not None:
        fs_res.to_csv(os.path.join(out_dir, args.pref+'.fscore.tsv'), sep='\t')

if __name__ == '__main__':
    main()
