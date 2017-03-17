''' Detect module for DriverPower
'''

import os
import sys
import logging
import pandas as pd
import numpy as np
from pybedtools import BedTool
from sklearn.preprocessing import StandardScaler, RobustScaler
from statsmodels.sandbox.stats.multicomp import multipletests
from driverpower.load import load_testFile, load_fselect, load_hdf5, load_covar, load_mut_bed, load_func_scores
from driverpower.model import get_model, estimate_bgmr, do_binom_test
from driverpower.func_adj import func_adj_new

logger = logging.getLogger('DETECT')


def getMutCtCg(mut, bed, mut_cnames, bed_cnames):
    ''' Get mutations in bed and mutation counts per element.
    Args:
        mut_cnames - list, column names of mut
        bed_cnames - list, column names of bed
    Return:
        mut_in, ct, cg, recur
    '''
    # mut table
    colnames = mut_cnames + bed_cnames
    mut_in = mut.intersect(bed, wa=True, wb=True)\
        .to_dataframe(names=colnames, dtype={'chrom': str, 'chrom_bin': str})
    ct = mut_in.pivot_table(index='binID', values='chrom', aggfunc=len)
    bed_df = bed.to_dataframe(names=bed_cnames, dtype={'chrom_bin': str})
    bed_df['length'] = bed_df.end_bin - bed_df.start_bin
    cg = bed_df.pivot_table(index='binID', values='length', aggfunc=sum)
    recur = mut_in.loc[:, ('binID', 'sid')]\
        .drop_duplicates().\
        pivot_table(index='binID', values='sid', aggfunc=len)
    return mut_in, ct, cg, recur


def calc_bin_fscore(mut, func_cols, agg_method):
    ''' Calculate functional score per bin
    Args:
        mut - pd.DF
        func_cols - list, colnames of func scores
        agg_method - str, aggregate method from ['mean', 'maxpool', 'meanpool']
    Return:
        fscores - pd.DF, indexed by binID, each column is one type of score
    '''
    fscores = pd.DataFrame(index=mut.binID.unique(), columns=func_cols)
    for ix in func_cols:
        num_na = mut[ix].isnull().sum()
        if num_na > 0:
            logger.warning('{} ({:.2f}%) NA values found in {} scores. All NAs are ignored'\
                           .format(num_na, num_na/mut.shape[0]*100, ix))
            mut_nafree = mut[mut[ix].notnull()] # remove NA variants instead of filling 0
        else:
            mut_nafree = mut
        if agg_method == 'mean':
            # mean of scores per element
            fscores[ix] = mut_nafree.pivot_table(index='binID', values=ix, aggfunc='mean')
        elif agg_method == 'maxpool':
            # max score per sample per element
            tmp = mut_nafree.pivot_table(index=['binID', 'sid'], values=ix, aggfunc='max').reset_index()
            # mean of max scores per element
            fscores[ix] = tmp.pivot_table(index='binID', values=ix, aggfunc='mean')
        elif agg_method == 'meanpool':
            # mean score per sample per element
            tmp = mut_nafree.pivot_table(index=['binID', 'sid'], values=ix, aggfunc='mean').reset_index()
            # mean of mean scores per element
            fscores[ix] = tmp.pivot_table(index='binID', values=ix, aggfunc='mean')
        else:
            logger.error('Functional score aggregate method not recognized (You enter: {})'.format(agg_method))
            sys.exit(1)
    return fscores


def format_res(res, func_cols):
    ''' Format output result
    Args:
        res - pd.DF, unformatted
        func_cols - set, set of functional names
    Return:
        res - pd.DF, formatted
    '''
    # data type
    res.length = res.length.astype(int)
    res.nMut = res.nMut.astype(int)
    res.nSample = res.nSample.astype(int)
    func_cols = list(func_cols)
    func_cols.sort()
    # get p.min and q.min
    if func_cols:
        pval_cols = ['p.' + name for name in func_cols]
        res['p.min'] = res[pval_cols].min(axis=1)
        not_na_pval = res['p.min'].notnull()
        res['q.min'] = np.nan
        res.loc[not_na_pval, 'q.min'] = multipletests(res['p.min'][not_na_pval], method='fdr_bh')[1]
        # sort
        res.sort_values('q.min', inplace=True)
        logger.info('Find {} elements with q-value <=  0.1'.format(sum(res['q.min']<=0.1)))
        func_res_names = []
        for name in func_cols:
            func_res_names += [name, 'p.'+name, 'q.'+name]
        res = res[['length', 'nMut', 'nSample', 'BGMR', 'p.raw', 'q.raw']+func_res_names+['p.min','q.min']]
    else:
        logger.info('Find {} elements with q-value <=  0.1'.format(sum(res['q.raw']<=0.1)))
        res = res[['length', 'nMut', 'nSample', 'BGMR', 'p.raw', 'q.raw']]
    return res


def detect(mut_path, callable_path, testFile_path, trainH5_path,
           fselect_path, fselect_name, fselect_cutoff,
           func_conf_path, agg_method, use_gmean, scaler_type,
           tumor_name, out_dir):
    ''' Main wrapper for detection
    Args:
        mut_path  - str, path to the mutation table
        callable_path - str, path to the whitelist regions
        testFile_path - str, path to the test file list
        trainH5_path   - str, path to the training data
        fselect_path   - str, path to the feature selection file
        fselect_name   - str, name of the feature selection method
        fselect_cufoff - float, cutoff for the feature selection
        func_conf_path - str, path to the functional scores location file
        agg_method - str, method used in calculating element score
        use_gmean - bool, use gmean of nMut and nSample as response if True
        scaler_type - str, feature scaling method
        tumor_name - str, cohort name used in output files
        out_dir - str, output file directory
    '''
    result_list = []
    mut_df, ndonor = load_mut_bed(mut_path)
    mut_bed = BedTool.from_dataframe(mut_df, na_rep='NA')
    if callable_path:
        callable_bed = BedTool(callable_path)
        mut_bed = mut_bed.intersect(callable_bed, wa=True)
        ncall = mut_bed.count()
        logger.info('{} ({:.2f}%) mutations are in callable regions'\
                    .format(ncall, ncall/mut_df.shape[0]*100))
    # Process feature selection result
    usefeatures = load_fselect(fselect_path, fselect_name, fselect_cutoff) if fselect_path else None
    # load training data and train the model
    logger.info('Load training data and train the model...')
    Xtrain_df, ytrain_df, Ntrain = load_hdf5(trainH5_path, usefeatures)
    train_index = Xtrain_df.index.values
    train_columns = Xtrain_df.columns.values
    Xtrain_mat = Xtrain_df.as_matrix()
    # scaling if need
    if scaler_type == 'robust':
        logger.info('Use robust scaler')
        scaler = RobustScaler()
        scaler.fit(Xtrain_mat)
        Xtrain_mat = scaler.transform(Xtrain_mat)
    elif scaler_type == 'standard':
        logger.info('Use standard scaler')
        scaler = StandardScaler()
        scaler.fit(Xtrain_mat)
        Xtrain_mat = scaler.transform(Xtrain_mat)
    else:
        scaler = None
    # train the model
    model = get_model(Xtrain_mat, ytrain_df, Ntrain, use_gmean, 'glm')
    # load test file list
    test_files = load_testFile(testFile_path, mut_df.columns.values, func_conf_path)
    bed_cnames = ['chrom_bin', 'start_bin', 'end_bin', 'binID']
    nset = len(test_files)
    curr_set = 1
    for test_set in test_files:
        (name, bin_path, feature_path, func_tuples) = test_set
        binIDs = pd.read_table(bin_path, sep='\t', header=None,\
                               names=bed_cnames, usecols=['binID'])
        binIDs = binIDs.binID.unique()
        print('===== Test Set [{}/{}]: {} (n={}) ====='.format(curr_set, nset, name, binIDs.shape[0]))
        curr_set += 1
        # result table
        result = pd.DataFrame(index=binIDs, columns=['length', 'nMut', 'nSample'])
        result.index.names = ['binID']
        # find mut in bin
        bin_bed = BedTool(bin_path)
        if callable_path:
            bin_bed = bin_bed.intersect(callable_bed)
        mut_cnames = list(mut_df.columns.values)
        mut, ct, cg, recur = getMutCtCg(mut_bed, bin_bed, mut_cnames, bed_cnames)
        logger.info('Number of mutations in set: {}'.format(mut.shape[0]))
        # add to result
        result['length'] = cg
        result['nSample'] = recur
        result['nMut'] = ct
        # fill na with 0
        result['length'].fillna(0, inplace=True)
        result['nSample'].fillna(0, inplace=True)
        result['nMut'].fillna(0, inplace=True)
        # load features
        Xtest_df = load_covar(feature_path, usefeatures)
        # estimate BGMR
        Xtest_df.sort_index(inplace=True)
        result.sort_index(inplace=True)
        # check names
        assert np.array_equal(train_index, ytrain_df.index), 'X and y indexes of training data do not match'
        assert np.array_equal(Xtest_df.index, result.index), 'X and y indexes of test data do not match'
        assert np.array_equal(train_columns, Xtest_df.columns), 'Train and test feature names do not match'
        if scaler_type in ['robust', 'standard']:
            Xtest_mat = scaler.transform(Xtest_df.as_matrix())
        else:
            Xtest_mat = Xtest_df.as_matrix()
        result['BGMR'] = estimate_bgmr(model, Xtest_mat, 'glm')
        p_raw, q_raw = do_binom_test(result, ndonor, result['BGMR'], use_gmean)
        result['p.raw'] = p_raw
        result['q.raw'] = q_raw
        result.sort_values('p.raw', inplace=True)
        # func adj
        if func_tuples:
            # retrive scores if necessary
            func_cols = set([i[0] for i in func_tuples])
            to_retrived = func_cols.difference(mut.columns)
            if to_retrived:
                mut = load_func_scores(func_conf_path, to_retrived, mut)
            # bin level functional scores
            fscores = calc_bin_fscore(mut, func_cols, agg_method)
            for ix, score in fscores.iteritems():
                result[ix] = score
            for fcol in func_tuples:
                pvals, qvals = func_adj_new(result, fcol, ndonor, use_gmean)
                result['p.'+fcol[0]] = pvals
                result['q.'+fcol[0]] = qvals
        else:
            func_cols = None
        # write result
        if tumor_name:
            out_path = os.path.join(out_dir, tumor_name + '.' + name + '.DriverPower.res.tsv')
        else:
            out_path = os.path.join(out_dir, name + '.DriverPower.res.tsv')
        res = format_res(result, func_cols)
        res.to_csv(out_path, na_rep='NA', sep='\t')
