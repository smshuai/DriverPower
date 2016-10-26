"""
Load data from DriverPower.
"""

import pandas as pd
import numpy as np
import logging

def load_coverage(path_coverage):
    ''' Load covergae table. Return a pd.DF with binID as the key
    '''
    # logging.info('Loading coverage table')
    cg = pd.read_csv(path_coverage, sep='\t', header=0, index_col='binID')
    for col in ['chrom', 'start', 'end']:
        if col in cg.columns.values:
            del cg[col]
    cg.sort_index(inplace=True)
    # check unique binID
    assert len(cg.index.values) == len(cg.index.unique()), "binID in coverage table is not unique."
    # check column number
    assert cg.shape[1] == 32, "Number of triple nucleotide contexts for coverage table is not 32." # ncol = 32
    logging.info('Successfully load coverages for {} bins'.format(cg.shape[0]))
    return cg

def load_covar(path_covar):
    ''' Load pre-processed covar. table. Return a pd.DF with binID as the key
    '''
    # logging.info('Loading feature table')
    cv = pd.read_csv(path_covar, sep='\t', header=0, index_col='binID')
    # check unique binID
    assert len(cv.index.values) == len(cv.index.unique()), "binID in feature table is not unique."
    cv.sort_index(inplace=True)
    na_count = cv.isnull().sum()
    if na_count.sum() > 0:
        na_fnames = na_count.index.values[na_count>0]
        logging.warning('NA values found in features: {}'.format(', '.join(na_fnames)))
        logging.warning('Fill NA with 0')
        cv.fillna(0, inplace=True)
    logging.info('Successfully load features for {} bins'.format(cv.shape[0]))
    return cv

def load_count(path_count):
    ''' Load mutation count data.
    Input 4 cols TSV file ['binID', 'sid', 'categ', 'ct']
    Output pd.DF without index
    '''
    # logging.info('Loading count table')
    ct = pd.read_csv(path_count, sep='\t', header=0)
    # check column names
    cname_check=np.array_equal(np.sort(ct.columns), ['binID', 'categ', 'ct', 'sid'])
    assert cname_check, 'Header of count table is {}. Need binID, categ, ct, sid'.format(", ".join(str(i) for i in ct.columns))
    ct.sort_values('binID', inplace=True)
    nmut = ct.ct.sum()
    nsample = len(ct.sid.unique())
    logging.info('Successfully load counts for {} mutations across {} samples'.format(nmut, nsample))
    return ct

def load_all(path_cg_test, path_ct_test, path_cv_test,
    path_cg_train, path_ct_train, path_cv_train):
    ''' Load train and test data
    '''
    logging.info('Loading test data')
    cg_test = load_coverage(path_cg_test)
    ct_test = load_count(path_ct_test)
    cv_test = load_covar(path_cv_test)
    assert np.array_equal(cv_test.index, cg_test.index), 'binIDs in test feature and coverage tables do not match'
    logging.info('Loading train data')
    cg_train = load_coverage(path_cg_train)
    ct_train = load_count(path_ct_train)
    cv_train = load_covar(path_cv_train)
    assert np.array_equal(cv_train.index, cg_train.index), 'binIDs in train feature and coverage tables do not match'
    return (cg_test, ct_test, cv_test, cg_train, ct_train, cv_train)