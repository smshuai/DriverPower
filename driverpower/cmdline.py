import argparse
import os
import logging
import sys
import pandas as pd
import numpy as np
from driverpower.load import load_memsave, load_mut
from driverpower.preprocess import get_response, scaling, sampling
from driverpower.feature_select import run_lasso, run_rndlasso, run_spearmanr, run_fregression
from driverpower.feature_select import feature_score
from driverpower import __version__
from driverpower.model import model, get_gmean
from driverpower.func_adj import func_adj


# logging config
logging.basicConfig(stream=sys.stdout, level=logging.INFO,
    format='%(asctime)s | %(name)s | %(levelname)s: %(message)s',
    datefmt='%m/%d/%Y %H:%M:%S')
# create logger
logger = logging.getLogger('DP')


def get_args():
    logger.info('DriverPower {}'.format(__version__))
    parser = argparse.ArgumentParser(prog='driverpower')
    # global argument
    parser.add_argument('-v', '--version', dest='version', action="store_true",
        help='Print the version of DriverPower')
    subparsers = parser.add_subparsers(title='The DriverPower sub-commands include',
        dest='subcommand')
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
    op_parser.add_argument('--sampling',
        type=float, dest='sampling', default=1.0,
        help='Number > 0 (default: 1). Sampling the data based on the provided value. Value in (0,1] is used as a fraction. Value > 1 is used as the number of data points.')
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
    op_select.add_argument('--gmean', dest='is_gmean', action="store_true",
        help='Use geometric mean of nMut and nSample as response')
    op_select.add_argument('-o', '--output', dest='out', type=str, default='feature_select.tsv',
        help='Path to the output file (default: ./feature_select.tsv)')
    #
    # Model
    #
    parser_model = subparsers.add_parser('model', help='Find driver bins with preprocessed training and test data')
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

    # optinal parameters
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
    # check for preprocess
    if args.subcommand == 'preprocess':
        # check input files
        check_file(args.path_ct)
        check_file(args.path_cg)
        check_file(args.path_cv)
        # check output file
        args.out = check_out(args.out)
        # check sampling value
        if args.sampling < 0:
            logger.error('Sampling value must be greater than 0. You enter {}'.format(args.sampling))
            sys.exit(1)
    elif args.subcommand == 'select':
        # check input file
        check_file(args.path_data)
        # check output file
        args.out = check_out(args.out)
        # check sampling
        if args.sampling < 0:
            logger.error('Sampling value must be greater than 0. You enter {}'.format(args.sampling))
            sys.exit(1)
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
            if args.func_cutoff < 1 or args.func_cutoff > 99:
                logger.error('--func_cutoff must be an int between 1 and 99. You enter {}'.format(args.func_cutoff))
                sys.exit(1)
        # check output file
        check_out(args.out)
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
    # initial output HDF5
    store = pd.HDFStore(args.out, mode='w')
    ct, cg, cv, grecur = load_memsave(args.path_ct,
        args.path_cg, args.path_cv,
        args.len_threshold, args.recur_threshold)
    # sample IDs
    sid = pd.Series(ct.sid.unique())
    sid.name = 'sid'
    # get response
    ybinom = get_response(ct, cg)
    # y to pd.DF
    ybinom = pd.DataFrame(ybinom, columns=['ct','len_ct'], index=cg.index)
    # sampling
    cv, ybinom = sampling(cv, ybinom, args.sampling)
    # write to store
    store.append('X', cv, chunksize=50000)
    store['y'] = ybinom
    store['grecur'] = grecur
    store['sid'] = sid
    store.close()
    logger.info('Pre-process done!')


def run_select(args):
    logger.info('Sub-command Select'.format(__version__))
    # load data from HDF5
    X = pd.read_hdf(args.path_data, 'X')
    y = pd.read_hdf(args.path_data, 'y')
    recur = pd.read_hdf(args.path_data, 'grecur')
    # check index (binID)
    assert np.array_equal(X.index, y.index), 'X and y have different row indexes'
    logger.info('Successfully load X with shape: {}'.format(X.shape))
    logger.info('Successfully load y with shape: {}'.format(y.shape))
    # Sampling data
    X, y = sampling(X, y, args.sampling)
    # silent delete logCG if exist
    if 'logCG' in X.columns.values:
        del X['logCG']
    # y to np.array
    y = y.as_matrix()
    if args.is_gmean:
        y = get_gmean(y, recur)
    # feature names
    fnames = X.columns.values
    # scale for Xtrain only
    X = X.as_matrix()
    logger.info('Training set X shape: {}'.format(X.shape))
    logger.info('Training set y shape: {}'.format(y.shape))
    X = scaling(Xtrain=X, scaler_type=args.scaler)
    # spearmanr
    rho = run_spearmanr(X, y)
    # f regression
    freg = run_fregression(X, y)
    # run LassoCV
    lasso = run_lasso(X, y)
    # run rndlasso
    rndlasso = run_rndlasso(X, y, lasso.alpha_)
    # results
    res = pd.DataFrame(np.array([rho, freg, lasso.coef_, rndlasso.scores_]).T,
        index=fnames,
        columns=['rho', 'freg','lasso', 'rndlasso'])
    res.index.name = 'fname'
    res.to_csv(args.out, sep='\t')
    logger.info('Feature selection done!')


def run_model(args):
    logger.info('Sub-command - Model'.format(__version__))
    # load training data
    Xtrain = pd.read_hdf(args.path_train, 'X')
    ytrain = pd.read_hdf(args.path_train, 'y')
    brecur = pd.read_hdf(args.path_train, 'grecur')
    assert np.array_equal(Xtrain.index, ytrain.index), 'Training X and y have different row indexes'
    logger.info('Successfully load X train with shape: {}'.format(Xtrain.shape))
    logger.info('Successfully load y train with shape: {}'.format(ytrain.shape))
    # load test data
    Xtest = pd.read_hdf(args.path_test, 'X')
    ytest = pd.read_hdf(args.path_test, 'y')
    grecur = pd.read_hdf(args.path_test, 'grecur')
    assert np.array_equal(Xtest.index, ytest.index), 'Test X and y have different row indexes'
    assert np.array_equal(grecur.index, ytest.index), 'Test recur and y have different row indexes'
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
            logger.info('At cutoff={}, {} selected features are: {}'.format(args.cutoff, len(fset), ", ".join(fset)))
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
    res.to_csv(args.out, sep='\t')
    logger.info('Model done!')


def main():
    args = get_args()
    if args.subcommand == 'preprocess':
        run_preprocess(args)
    elif args.subcommand == 'select':
        run_select(args)
    elif args.subcommand == 'model':
        run_model(args)


if __name__ == '__main__':
    main()
