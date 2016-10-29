"""
Load data from DriverPower.
"""

import pandas as pd
import numpy as np
import logging


logger = logging.getLogger('LOAD')
# logger.setLevel(logging.INFO)
# # create console handler with a higher log level
# ch = logging.StreamHandler()
# ch.setLevel(logging.INFO)
# # create formatter and add it to the handlers
# formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s',
#     datefmt='%m/%d/%Y %H:%M:%S')
# ch.setFormatter(formatter)
# # add the handlers to the logger
# logger.addHandler(ch)

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
    logger.info('Successfully load coverages for {} bins'.format(cg.shape[0]))
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
        na_fnames = na_count.index.values[np.where(na_count>0)]
        logger.warning('NA values found in features [{}]'.format(', '.join(na_fnames)))
        logger.warning('Fill NA with 0')
        cv.fillna(0, inplace=True)
    logger.info('Successfully load features for {} bins'.format(cv.shape[0]))
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
    logger.info('Successfully load counts for {} mutations across {} samples'.format(nmut, nsample))
    return ct


def load_mut(path_mut):
    ''' Load mutation table.
    Eight columns are required ['chrom', 'start', 'end', 'type', 'ref', 'alt', 'sid', 'binID']
    Output pd.DF without index
    '''
    if path_mut is None:
        return None # read nothing
    mut = pd.read_table(path_mut, sep='\t', header=0, index_col=None)
    # check column names
    need_cols = pd.Series(['chrom', 'start', 'end', 'type', 'ref', 'alt', 'sid', 'binID'])
    assert np.sum(~need_cols.isin(mut.columns))==0, 'Mutation table needs the following columns: {}'.format(", ".join(need_cols))
    mut = mut.loc[:,need_cols] # only keep need cols
    logger.info('Successfully load {} mutations'.format(mut.shape[0]))
    return mut


def load_all(path_cg_test, path_ct_test, path_cv_test, path_mut,
    path_cg_train, path_ct_train, path_cv_train):
    ''' Load all train and test data
    '''
    logger.info('Loading test data')
    cg_test = load_coverage(path_cg_test)
    ct_test = load_count(path_ct_test)
    cv_test = load_covar(path_cv_test)
    mut = load_mut(path_mut) # only load mut for test data
    assert np.array_equal(cv_test.index, cg_test.index), 'binIDs in test feature and coverage tables do not match'
    logger.info('Loading train data')
    cg_train = load_coverage(path_cg_train)
    ct_train = load_count(path_ct_train)
    cv_train = load_covar(path_cv_train)
    assert np.array_equal(cv_train.index, cg_train.index), 'binIDs in train feature and coverage tables do not match'
    assert np.array_equal(cv_train.columns, cv_test.columns), 'Feature names in train and test sets do not match'
    fnames = cv_train.columns.values
    return (cg_test, ct_test, cv_test.as_matrix(), mut,
        cg_train, ct_train, cv_train.as_matrix(), fnames)
