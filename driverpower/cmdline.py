import argparse
import os
import logging
import sys
import pandas as pd
import numpy as np
from driverpower import __version__
from driverpower.load import load_mut
from driverpower.preprocess import preprocess_v1
from driverpower.feature_select import fselect_v1
from driverpower.model import model, get_gmean
from driverpower.func_adj import func_adj
from driverpower.detect import detect

# logging config
# logging.basicConfig(stream=sys.stdout, level=logging.INFO,
#    format='%(asctime)s | %(name)s | %(levelname)s: %(message)s',
#    datefmt='%m/%d/%Y %H:%M:%S')
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                    format='%(asctime)s | %(levelname)s: %(message)s',
                    datefmt='%m/%d/%Y %H:%M:%S')
# create logger
logger = logging.getLogger('DP')


def get_args():
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                          argparse.MetavarTypeHelpFormatter):
        pass
    parser = argparse.ArgumentParser(prog='driverpower', description='DriverPower v{}: Combined burden and functional impact test for coding and noncoding cancer elements'.format(__version__))
    # global argument
    parser.add_argument('-v', '--version', dest='version', action="store_true",
        help='Print the version of DriverPower')
    subparsers = parser.add_subparsers(title='The DriverPower sub-commands include',
        dest='subcommand')
    #
    # Load and preprocess data
    #
    parser_preprocess = subparsers.add_parser('preprocess',
        help='Load and preprocess data', formatter_class=CustomFormatter)
    # required parameters
    re_parser = parser_preprocess.add_argument_group(title="required arguments")
    # required 3 tables for test data
    re_parser.add_argument('--variant', dest='mut_path', required=True, type=str,
        help='Path to the variant table')
    re_parser.add_argument('--feature', dest='feature_path', required=True, type=str,
        help='Path to the feature table')
    re_parser.add_argument('--element', dest='bin_path', required=True, type=str,
        help='Path to the element set (BED)')
    # Optional parameters
    op_parser = parser_preprocess.add_argument_group(title="optional parameters")
    op_parser.add_argument('--callable', dest='callable_path', required=False, type=str,
        help='Path to the whitelist regions', default=None)
    # op_parser.add_argument('--len_threshold', dest='len_threshold', type=int, default=500,
    #     help='Integer (default: 500). Bins with length < len_threshold will be discarded')
    # op_parser.add_argument('--recur_threshold', dest='recur_threshold', type=int, default=2,
    #     help='Integer (default: 2). Bins having mutations in < recur_threshold samples will be discarded')
    # op_parser.add_argument('--sampling',
    #     type=float, dest='sampling', default=1.0,
    #     help='Number > 0 (default: 1). Sampling the data based on the provided value. Value in (0,1] is used as a fraction. Value > 1 is used as the number of data points.')
    op_parser.add_argument('--output', dest='out_path', type=str, default='train.h5',
        help='Path to the output file')
    #
    # Feature selection
    #
    parser_select = subparsers.add_parser('select',
        help='Run feature selection on preprocessed data', formatter_class=CustomFormatter)
    # required parameters
    re_select = parser_select.add_argument_group(title="required arguments")
    re_select.add_argument('--trainH5', dest='h5_path', required=True, type=str,
        help='Path to the preprocessed training set (HDF5)')
    # optional parameters
    op_select = parser_select.add_argument_group(title="optional parameters")
    op_select.add_argument('--scaler', choices=['robust', 'standard', None],
        type=str, dest='scaler_type', default='robust',
        help='Scaler used to scale the features')
    # op_select.add_argument('--sampling',
    #     type=float, dest='sampling', default=1.0,
    #     help='Number > 0 (default: 1). Sampling the data based on the provided value. Value in (0,1] is used as a fraction. Value > 1 is used as the number of data points.')
    op_select.add_argument('--noGmean', dest='no_gmean', required=False, action="store_true",
        help='Do not geometric mean of nMut and nSample as response')
    op_select.add_argument('--output', dest='out_path', type=str, default='feature_select.tsv',
        help='Path to the output file')
    #
    # Model
    #
    parser_model = subparsers.add_parser('model', help='Find driver bins with preprocessed training and test data (deprecated)')
    # required parameters
    re_model = parser_model.add_argument_group(title="required arguments")
    re_model.add_argument('--train', dest='path_train', required=True, type=str,
        help='Path to the preprocessed training set (HDF5)')
    re_model.add_argument('--test', dest='path_test', required=True, type=str,
        help='Path to the preprocessed test set (HDF5)')
    # optional parameters
    ## other options
    op_model = parser_model.add_argument_group(title="optional parameters")
    op_model.add_argument('--coding', dest='is_coding', action="store_true",
        help='Test for coding bins')
    op_model.add_argument('--fold',
        type=int, dest='fold', default=1,
        help='Int (default: 1). Split the data into k-folds by effective length')
    op_model.add_argument('--gmean', dest='is_gmean', action="store_true",
        help='Use geometric mean of nMut and nSample as response')
    op_model.add_argument('--scaler', choices=['robust', 'standard', 'none'],
        type=str, dest='scaler', default='robust',
        help='robust or standard (default: robust). Scaler used to scale the data')
    op_model.add_argument('-o', '--output', dest='out', type=str, default='driverpower_result.tsv',
        help='Path to the output file (default: ./driverpower_result.tsv)')
    ## feature select
    op_select_model = parser_model.add_argument_group(title="optional parameters for feature selection")
    op_select_model.add_argument('--select', dest='path_select', type=str,
        help='Path to the feature selection table')
    op_select_model.add_argument('--select_criteria', dest='criteria', type=str, default='rndlasso',
        help='(default: rndlasso). Feature selection criteria')
    op_select_model.add_argument('--select_cutoff', dest='cutoff', type=float, default=0.5,
        help='(default: 0.5). Feature selection cutoff')
    ## functional adjustment
    op_func_model = parser_model.add_argument_group(title="optional parameters for functional adjustment")
    op_func_model.add_argument('--func', choices=['cadd', 'eigen', 'none'],
        type=str, dest='func', default='none',
        help='(default: none). Type of functional scores to use')
    op_func_model.add_argument('--func_dir', dest='dir_func', type=str, default='~/dp_func/',
        help='Directory of functional scores (default: ~/dp_func/)')
    op_func_model.add_argument('--func_cutoff', dest='funcadj', type=int, default=85,
        help='Integer between 1 and 99 (default: 85). Strength of functional adjustment.')
    op_func_model.add_argument('--mut', dest='path_mut', type=str,
        help='Path to the mutation table')
    #
    # Detect
    #
    parser_detect = subparsers.add_parser(
        'detect', help='Detect cancer driver elements',
        description='Detect cancer driver elements', formatter_class=CustomFormatter)
    # required parameters
    re_detect = parser_detect.add_argument_group(title="required parameters")
    re_detect.add_argument('--variant', dest='path_mut', required=True, type=str,
        help='Path to the variant table', default=argparse.SUPPRESS)
    re_detect.add_argument('--testFile', dest='path_test', required=True, type=str,
        help='Path to the test file list', default=argparse.SUPPRESS)
    re_detect.add_argument('--trainH5', dest='path_train', required=True, type=str,
        help='Path to the training HDF5', default=argparse.SUPPRESS)
    # optional parameters
    op_detect = parser_detect.add_argument_group(title="optional parameters")
    op_detect.add_argument('--callable', dest='path_call', required=False, type=str,
        help='Path to the whitelist regions', default=None)
    op_detect.add_argument('--funcPool', dest='func_pool', required=False, type=str,
        choices=['mean', 'meanpool', 'maxpool'], default='mean',
        help='Pooling method used to calculate functional impact score per element')
    op_detect.add_argument('--funcConf', dest='path_conf', required=False, type=str,
        help='Path to the functional score configuration file', default=None)
    op_detect.add_argument('--selectPath', dest='select_path', required=False, type=str,
        help='Path to the feature selection result', default=None)
    op_detect.add_argument('--selectName', dest='select_name', required=False, type=str,
        help='Feature selection method name', default=None)
    op_detect.add_argument('--selectCut', dest='select_cut', required=False, type=float,
        help='Feature importance cutoff', default=0)
    op_detect.add_argument('--scaler', choices=['robust', 'standard', None],
        type=str, dest='scaler', default='robust',
        help='Scaler used to scale the features')
    op_detect.add_argument('--noGmean', dest='no_gmean', required=False, action="store_true",
        help='Do not geometric mean of nMut and nSample as response')
    op_detect.add_argument('--cohortName', dest='tumor_name', required=False,
        help='Cohort name to use in output files', default=None, type=str)
    op_detect.add_argument('--outDir', dest='out_dir', required=False,
        type=str, default='./',
        help='Directory of output files')
    args = parser.parse_args()
    #
    # Parameters check
    #
    # no argument, print main help instead
    if len(sys.argv)==1:
        parser.print_help()
        parser.exit(1)
    if args.version:
        print("DriverPower", __version__)
        sys.exit(0)
    # check for preprocess
    if args.subcommand == 'preprocess':
        # check input files
        check_file(args.mut_path)
        check_file(args.feature_path)
        check_file(args.bin_path)
        # check output file
        args.out_path = check_out(args.out_path)
        # check sampling value
        # if args.sampling < 0:
        #     logger.error('Sampling value must be greater than 0. You enter {}'.format(args.sampling))
        #     sys.exit(1)
    # check for select
    elif args.subcommand == 'select':
        # check input file
        check_file(args.h5_path)
        # check output file
        args.out_path = check_out(args.out_path)
        # check sampling
        # if args.sampling < 0:
        #     logger.error('Sampling value must be greater than 0. You enter {}'.format(args.sampling))
        #     sys.exit(1)
    # check for model
    elif args.subcommand == 'model':
        # check input file
        check_file(args.path_train)
        check_file(args.path_test)
        if args.path_select is not None:
            # use selection
            check_file(args.path_select)
        if args.func != 'none':
            # use func adj
            if args.path_mut is None:
                logger.error('Please specify mutation table (--mut) for functional adjustment')
                sys.exit(1)
            else:
                check_file(args.path_mut)
            if args.funcadj < 1 or args.funcadj > 99:
                logger.error('--func_cutoff must be an int between 1 and 99. You enter {}'.format(args.func_cutoff))
                sys.exit(1)
        # check output file
        args.out = check_out(args.out)
        # check fold
        if args.fold < 0:
            logger.error('Fold value must be greater than 0. You enter {}'.format(args.fold))
            sys.exit(1)
    return args


def check_file(path):
    if not os.path.isfile(os.path.expanduser(path)):
        # no such file
        logger.error('File not found in {}'.format(path))
        sys.exit(1)
    else:
        return True


def check_out(path):
    ''' Check the output file path.
    Return a path that is not used and creatable
    '''
    path = os.path.expanduser(path)
    if os.path.isfile(path):
        # exist already?
        logger.warning('The output file {} already exists'.format(path))
        i = 1
        while os.path.isfile(path):
            path += '.' + str(i)
            i += 1
        logger.warning('Use {} instead'.format(path))
    else:
        # creatable?
        try:
            open(path, 'w').close()
        except OSError:
            logger.error('The output path {} is not valid'.format(path))
            sys.exit(1)
    return path


def run_preprocess(args):
    logger.info('Sub-command Preprocess'.format(__version__))
    # version 1 - preprocess
    preprocess_v1(mut_path=args.mut_path, callable_path=args.callable_path, bin_path=args.bin_path,
                  feature_path=args.feature_path, out_path=args.out_path)
    logger.info('Pre-process done!')


def run_select(args):
    logger.info('Sub-command Select'.format(__version__))
    use_gmean = False if args.no_gmean else True
    # use train h5
    fselect_v1(h5_path=args.h5_path, scaler_type=args.scaler_type, use_gmean=use_gmean, out_path=args.out_path)
    logger.info('Feature selection done!')


def run_model(args):
    logger.info('Sub-command - Model'.format(__version__))
    # load training data
    ytrain = pd.read_hdf(args.path_train, 'y')
    if ytrain.shape[0] == 0:
        logger.error('No training data in hdf5 file')
        sys.exit(1)
    Xtrain = pd.read_hdf(args.path_train, 'X')
    brecur = pd.read_hdf(args.path_train, 'recur')
    bsid = pd.read_hdf(args.path_train, 'sid')
    Ntrain = bsid.unique().shape[0]
    ytrain.len_ct = (ytrain.len_ct + ytrain.ct) * Ntrain - ytrain.ct
    assert np.array_equal(Xtrain.index, ytrain.index), 'Training X and y have different row indexes'
    assert np.array_equal(brecur.index, ytrain.index), 'Training recur and y have different row indexes'
    logger.info('Successfully find training data for {} samples'.format(Ntrain))
    logger.info('Successfully load X train with shape: {}'.format(Xtrain.shape))
    logger.info('Successfully load y train with shape: {}'.format(ytrain.shape))
    # load test data
    ytest = pd.read_hdf(args.path_test, 'y')
    if ytest.shape[0] == 0:
        logger.error('No test data in hdf5 file')
        sys.exit(1)
    Xtest = pd.read_hdf(args.path_test, 'X')
    grecur = pd.read_hdf(args.path_test, 'recur')
    glength = ytest.len_ct + ytest.ct
    gsid = pd.read_hdf(args.path_test, 'sid')
    Ntest = gsid.unique().shape[0]
    ytest.len_ct = (ytest.len_ct + ytest.ct) * Ntest - ytest.ct
    assert np.array_equal(Xtest.index, ytest.index), 'Test X and y have different row indexes'
    assert np.array_equal(grecur.index, ytest.index), 'Test recur and y have different row indexes'
    logger.info('Successfully find test data for {} samples'.format(Ntest))
    logger.info('Successfully load X test with shape: {}'.format(Xtest.shape))
    logger.info('Successfully load y test with shape: {}'.format(ytest.shape))
    # make sure fnames match
    Xtrain.sort_index(1, inplace=True)
    Xtest.sort_index(1, inplace=True)
    assert np.array_equal(Xtest.columns, Xtrain.columns), 'Training and test X have different feature names'
    # obtain feature selection
    if args.path_select is not None:
        is_select = True
        select_tb = pd.read_table(args.path_select, index_col='fname')
        if args.criteria in select_tb.columns.values:
            logger.info('Use {} as criteria in feature selection'.format(args.criteria))
            fset, fidx = feature_score(select_tb[args.criteria].abs(), select_tb.index.values, args.cutoff)
            logger.info('At cutoff={}, {} features are selected'.format(args.cutoff, len(fset)))
        else:
            logger.error('Feature selection criteria {} is not in selection table'.format(args.criteria))
            sys.exit(1)
    else:
        is_select = False
        logger.info('Use all features')
    # get gnames
    gnames  = ytest.index.values
    if is_select: # feature selection ON
        # select features based on names
        Xtrain = Xtrain.loc[:, fset]
        Xtest = Xtest.loc[:, fset]
    # scaling
    logger.info('Training set shapes: X {} and y {}'.format(Xtrain.shape, ytrain.shape))
    logger.info('Test set shapes: X {} and y {}'.format(Xtest.shape, ytest.shape))
    Xtrain, Xtest = scaling(Xtrain=Xtrain.as_matrix(),
        Xtest=Xtest.as_matrix(), scaler_type=args.scaler)
    # glm
    res = model(Xtrain, ytrain.as_matrix(), Xtest,
        ytest.as_matrix(), gnames, grecur, brecur,
        args.is_gmean, method='glm', fold=args.fold)
    # functional adjustment
    if args.func != 'none':
        # read mutation table
        mut = load_mut(args.path_mut)
        logger.info('Start functional adjustment')
        res = func_adj(res, mut, method=args.func,
            dir_func=os.path.expanduser(args.dir_func),
            is_coding=args.is_coding, cutoff=args.funcadj)
    res['Length'] = glength
    res.to_csv(args.out, sep='\t')
    logger.info('Model done!')


def run_detect(args):
    logger.info('Sub-command - Detect')
    use_gmean = False if args.no_gmean else True
    detect(args.path_mut, args.path_call, args.path_test,
          args.path_train, args.select_path, args.select_name,
          args.select_cut, args.path_conf, args.func_pool,
          use_gmean, args.scaler, args.tumor_name, args.out_dir)

def main():
    args = get_args()
    logger.info('DriverPower {}'.format(__version__))
    if args.subcommand == 'preprocess':
        run_preprocess(args)
    elif args.subcommand == 'select':
        run_select(args)
    elif args.subcommand == 'model':
        run_model(args)
    elif args.subcommand == 'detect':
        run_detect(args)


if __name__ == '__main__':
    main()
