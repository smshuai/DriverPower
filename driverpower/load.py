"""
Load data for DriverPower.
"""

import pandas as pd
import numpy as np
import tabix
import logging
import sys
from driverpower.feature_select import feature_score
from tabix import TabixError
from driverpower.helperDP import get_filter

logger = logging.getLogger('LOAD')

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


def load_covar(path_covar, usefeatures=None):
    ''' Load pre-processed covar. table. Return a pd.DF with binID as the key
    '''
    # logging.info('Loading feature table')
    if usefeatures:
        cv = pd.read_csv(path_covar, sep='\t', header=0, index_col='binID',
                         usecols=['binID'] + usefeatures)
    else:
        cv = pd.read_csv(path_covar, sep='\t', header=0, index_col='binID')
    # check unique binID
    assert len(cv.index.values) == len(cv.index.unique()), "binID in feature table is not unique."
    cv.sort_index(inplace=True) # sort row index
    cv.sort_index(1, inplace=True) # sort column index
    na_count = cv.isnull().sum()
    if na_count.sum() > 0:
        na_fnames = na_count.index.values[np.where(na_count>0)]
        logger.warning('NA values found in features [{}]'.format(', '.join(na_fnames)))
        logger.warning('Fill NA with 0')
        cv.fillna(0, inplace=True)
    logger.info('Successfully load {} features for {} bins'.format(cv.shape[1], cv.shape[0]))
    return cv


def load_count(path_count):
    ''' Load mutation count data.
    Input 3 cols TSV file ['binID', 'sid', 'ct'], other columns will be ignored
    Output pd.DF without index
    '''
    # logging.info('Loading count table')
    ct = pd.read_csv(path_count, sep='\t', header=0)
    # check column names
    need_cols = pd.Series(['binID', 'sid', 'ct'])
    assert np.sum(~need_cols.isin(ct.columns))==0, 'Count table needs the following columns: {}'.format(", ".join(need_cols))
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
    need_cols = ['chrom', 'start', 'end', 'type', 'ref', 'alt', 'sid', 'binID']
    try:
        mut = pd.read_table(path_mut, sep='\t', header=0, index_col=None,
                            usecols=need_cols,
                            dtype={'chrom': str})
    except ValueError:
        logger.error('Mutation table needs the following columns **with** header: {}'.format(", ".join(need_cols)))
        sys.exit(1)
    # check column names
    mut = mut.loc[:,need_cols] # reorder columns
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


def load_memsave(path_ct, path_cg, path_cv, len_threshold=500, recur_threshold=2):
    ''' Load CT, CG and CV in a memsave way. Return filtered CT, CG, CV and grecur
    '''
    cg = load_coverage(path_cg)
    ct = load_count(path_ct)
    # pre-filter CV
    keep, tab = get_filter(ct, cg, len_threshold=len_threshold,
                           recur_threshold=recur_threshold, return_tab=True)
    keep_bin = tab.index.values[keep]
    Nbin  = tab.shape[0]
    Nchunk = int(Nbin / 50000)
    chunk_idx = 1
    # read in chunk and filter
    logger.info('Start to load and filter features')
    cv_reader = pd.read_table(path_cv, index_col='binID', chunksize=50000)
    cv = cv_reader.get_chunk()
    cv = cv[cv.index.isin(keep_bin)] # cv container
    for chunk in cv_reader:
        logger.info('Load features chunk {}/{}'.format(chunk_idx, Nchunk))
        chunk_idx += 1
        chunk = chunk[chunk.index.isin(keep_bin)]
        cv = cv.append(chunk)
    # check unique binID
    assert len(cv.index.values) == len(cv.index.unique()), "binID in feature table is not unique."
    cv.sort_index(inplace=True) # sort row index
    cv.sort_index(1, inplace=True) # sort column index
    na_count = cv.isnull().sum()
    if na_count.sum() > 0:
        na_fnames = na_count.index.values[np.where(na_count>0)]
        logger.warning('NA values found in features [{}]'.format(', '.join(na_fnames)))
        logger.warning('Fill NA with 0')
        cv.fillna(0, inplace=True)
    logger.info('Successfully load features for {} bins'.format(cv.shape[0]))
    cg = cg[cg.index.isin(cv.index)] # filter cg
    assert np.array_equal(cv.index.sort_values(), cg.index.sort_values()), 'binIDs in feature and coverage tables do not match'
    assert np.array_equal(cv.index, cg.index), 'binIDs in feature and coverage tables are not sorted'
    ct = ct[ct.binID.isin(cv.index)] # filter ct as well
    # get recur
    grecur = tab.recur.loc[keep_bin]
    grecur.sort_index(inplace=True)
    assert np.array_equal(cv.index, grecur.index), 'binIDs in feature and recur tables are not equal'
    return ct, cg, cv, grecur


def load_all_memsave(path_cg_test, path_ct_test, path_cv_test, path_mut,
    path_cg_train, path_ct_train, path_cv_train, len_threshold, recur_threshold):
    ''' Load all data in a memsave way
    '''
    logger.info('Loading test data') # no change to test data
    ct_test, cg_test, cv_test = load_memsave(path_ct_test, path_cg_test, path_cv_test,
        len_threshold, recur_threshold)
    mut = load_mut(path_mut) # only load mut for test data
    logger.info('Loading train data')
    ct_train, cg_train, cv_train = load_memsave(path_ct_train, path_cg_train, path_cv_train,
        len_threshold, recur_threshold)
    assert np.array_equal(cv_train.columns, cv_test.columns), 'Feature names in train and test sets do not match'
    fnames = cv_train.columns.values
    return (cg_test, ct_test, cv_test.as_matrix(), mut,
        cg_train, ct_train, cv_train.as_matrix(), fnames)

###
# v0.5.0 Most functions added for DETECT
###

# chrom 1-22,X,Y
VALID_CHROMS = set(['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                    '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X', 'Y'])
## For the detect module
def load_mut_bed(mut_path):
    ''' Load mutation table with at least 6 columns (header=True).
    ['chrom', 'start', 'end', 'ref', 'alt', 'sid']
    The first three columns must be chrom, start, end
    Args:
        mut_path - str, path to the mutation table
    '''
    mut = pd.read_table(mut_path, sep='\t', header=0,
                        dtype={'chrom': str, 'start': np.int64, 'end': np.int64,
                               'ref': str, 'alt': str, 'sid': str})
    # check chrom names
    chroms = set(mut.chrom.unique())
    assert chroms.issubset(VALID_CHROMS), 'Mutation table contains invalid chromosome names. Valid names are 1-22,X,Y'
    # check ref != alt
    nerror = sum(mut.ref == mut.alt)
    assert nerror == 0, 'Mutation table contains {} row(s) where REF == ALT'.format(nerror)
    # check start end
    nerror = sum(mut.start >= mut.end)
    assert nerror == 0, 'Mutation table contains {} row(s) where start >= end'.format(nerror)
    # check ref and alt 'ATCG-'
    alphabet_ref = set(''.join(mut.ref.unique()))
    alphabet_alt = set(''.join(mut.alt.unique()))
    valid_letters = set('ATCG-')
    assert alphabet_ref.issubset(valid_letters), 'Mutation table contains invalid letters ({}) in REF. Valid letters are "ATCG-"'.format(alphabet_ref)
    assert alphabet_alt.issubset(valid_letters), 'Mutation table contains invalid letters ({}) in ALT. Valid letters are "ATCG-"'.format(alphabet_alt)
    # Assign DEL, INS and SNP
    ndonor = mut.sid.unique().shape[0]
    if 'type' not in mut.columns:
        mut = assign_variant_type(mut)
    logger.info('Successfully load {} mutations from {} donors'.format(mut.shape[0], ndonor))
    return mut, ndonor

def load_testFile(testFile_path, mut_cols, func_conf_path):
    ''' Load and check the test file lists
    Args:
        testFile_path - str, path to the test file list
        mut_cols - list, column names in mutation table
        func_conf_path - str, path to the functional scores location file
    Return:
        file_list - list, each element is one row in testFile. (name, element_path, feature_path, func_tuples or None). func_tuples are list of (func_name, func_thresh, func_cut)
    '''
    file_list = []
    testFile = pd.read_table(testFile_path, sep='\t', header=0)
    # at least three columns 
    assert set(testFile.columns).issuperset(['element', 'feature', 'name']), "The test files list need at least three columns ['element', 'feature', 'name']"
    # name must be unique
    assert testFile['name'].unique().shape[0] == testFile['name'].shape[0], 'Name in the test files list must be unique'
    # check if func_names provide
    mut_cols = [name.upper() for name in mut_cols]
    if 'func_names' in testFile.columns:
        func_names_uniq = set([name.upper() for i in testFile['func_names'] for name in i.split(",")])
        if func_names_uniq.difference(mut_cols):
            # some func_names not found in mut_cols, search func_conf.
            conf_names = pd.read_csv(func_conf_path, header=0, usecols=['name'])
            conf_names_uniq = set([name.upper() for name in conf_names['name'].unique()])
            left = func_names_uniq.difference(mut_cols).difference(conf_names_uniq)
            if left:
                logger.error("Cannot find {} functional scores in mutation table nor configuration file".format(left))
                sys.exit(1)
        # when func_names provide, check func_thresh and func_cut
        if 'func_thresh' in testFile.columns:
            for ix, row in testFile.iterrows():
                assert row['func_names'].count(',') == row['func_thresh'].count(','), 'func_thresh and func_names must have the same length for test name {}'.format(row['name'])
                func_tuples = [ (name.upper(), float(thresh), None) for name, thresh in zip(row['func_names'].split(","), row['func_thresh'].split(","))]
                file_list.append((row['name'], row['element'], row['feature'], func_tuples))
        elif 'func_cut' in testFile.columns:
            for ix, row in testFile.iterrows():
                assert row['func_names'].count(',') == row['func_cut'].count(','), 'func_cut and func_names must have the same length for test name {}'.format(row['name'])
                func_tuples = [ (name.upper(), None, float(cutoff)) for name, cutoff in zip(row['func_names'].split(","), row['func_cut'].split(","))]
                file_list.append((row['name'], row['element'], row['feature'], func_tuples))
        else:
            # use default 85
            for ix, row in testFile.iterrows():
                func_tuples = [ (name.upper(), 85, None) for name in row['func_names'].split(",")]
                file_list.append((row['name'], row['element'], row['feature'], func_tuples))
    else:
        # no func_names
        for ix, row in testFile.iterrows():
            file_list.append((row['name'], row['element'], row['feature'], None))
    return file_list


def assign_variant_type(mut):
    ''' Assign INS, DEL, SNP, other to mutations
    '''
    logger.info('Infer variant type (SNP, INS, DEL, other) from mutation table')
    length_ref = np.array([len(i) for i in mut.ref])
    length_alt = np.array([len(i) for i in mut.alt])
    length_var = np.array(mut.end - mut.start)
    is_ref_minus = np.array(mut.ref == '-')
    is_alt_minus = np.array(mut.alt == '-')
    # Find SNPs:
    # 1. length_ref == length_alt == 1
    # 2. ref != alt != '-'
    # 3. end - start == 1
    test1 = np.logical_and(length_ref==1, length_alt==1)
    test2 = np.logical_and(np.logical_not(is_ref_minus), np.logical_not(is_alt_minus))
    test3 = length_var == 1
    is_snp = [ x&y&z for x,y,z in zip(test1, test2, test3)]
    # Find INSs:
    # 1. ref == '-'
    # 2. end - start == 2
    is_ins = [ x&y for x,y in zip(is_ref_minus, length_var==2)]
    # Find DELs:
    # 1. alt == '-'
    # 2. end - start == length_ref
    is_del = [ x&y for x,y in zip(is_alt_minus, length_var==length_ref)]
    mut['type'] = 'other'
    mut.loc[is_snp, 'type'] = 'SNP'
    mut.loc[is_del, 'type'] = 'DEL'
    mut.loc[is_ins, 'type'] = 'INS'
    return mut


def load_hdf5(h5_path, usefeatures=None):
    ''' Load HDF5 data from preprocess. Keys are ['y', 'X', 'recur', 'sid']
    Args:
        h5_path - str, path to the hdf5 path
        usefeatures - list, a list of feature names to use. Default: None
    Return:
        X - pd.DF, indexed by binID, each column is a feature
        y - pd.DF, indexed by binID, columns are nMut, length, nSample
        N - int, number of samples
    '''
    y = pd.read_hdf(h5_path, 'y')
    X = pd.read_hdf(h5_path, 'X')
    if usefeatures:
        X = X.loc[:, usefeatures]
    recur = pd.read_hdf(h5_path, 'recur')
    sids = pd.read_hdf(h5_path, 'sid')
    N = sids.unique().shape[0]
    logger.info('Successfully load data for {} samples'.format(N))
    # check index (binID)
    assert np.array_equal(X.index, y.index), 'X and y have different row indexes'
    assert np.array_equal(y.index, recur.index), 'recur and y have different row indexes'
    logger.info('Successfully load X with shape: {}'.format(X.shape))
    logger.info('Successfully load y with shape: {}'.format(y.shape))
    y['length']  = y.len_ct + y.ct
    y['nSample'] = recur.astype(np.int_)
    del y['len_ct']
    y.columns = ['nMut', 'length', 'nSample']
    return X, y, N


def load_fselect(fs_path, fs_name, fs_cutoff):
    ''' Load the feature selection result
    Args:
        fs_path - str, path to the feature selection result
        fs_name - str, colname of the feature selection method to use
        fs_cutoff - float, cutoff for the feature selection
    Return:
        fset - list, a list of feature names that will be used in modelling
    '''
    fs = pd.read_table(fs_path, sep='\t', header=0, usecols=['fname', fs_name])
    fset, fidx = feature_score(fs[fs_name].abs(), fs.fname, fs_cutoff)
    logger.info('Use {} at cutoff={}, {} features are selected'.format(fs_name, fs_cutoff, len(fset)))
    return list(fset)


def load_func_scores(func_conf_path, methods, mut):
    ''' Load the configure file for functional scores
    Args:
        func_conf_path - str, path to the conf file
        methods - list, name of methods will be used
        mut - pd.DF, mutation table
    '''
    conf = pd.read_csv(func_conf_path, header=0,
                       dtype={'ref_ix':np.int_, 'alt_ix':np.int_, 'score_ix':np.int_})
    support_methods = set([i.upper() for i in conf.name.unique()]) # all in upper case
    methods = set([i.upper() for i in methods])
    assert support_methods.issuperset(methods), \
            'Cannot find configuration for the following method(s): {}'.format(', '.join(methods-support_methods))
    conf_in = conf[conf.name.str.upper().isin(methods)]
    # new columns to be added in the mut if not exist.
    uniq_names = [cname.upper() for cname in conf_in.name.unique()]
    for cname in uniq_names:
        if cname not in mut.columns:
            logger.info('Load functional score: {}'.format(cname))
            conf_cname = conf[conf.name == cname]
            mut[cname] = retrive_score(mut, conf_cname)
        else:
            logger.info('Skip loading {}. Already exist in mutation table.'.format(cname))
    return mut

def retrive_score(mut, conf):
    ''' Obtain functional scores based on mut and conf
    '''
    conf = conf.sort_values('order')
    score = np.empty(shape=mut.shape[0])
    score[:] = np.NAN
    for ix_conf, conf_row in conf.iterrows():
        # logger.info('Retriving {} - {} - chrom {}'.format(conf_row['name'], conf_row['type'], conf_row['chroms']))
        tb = tabix.open(conf_row['path'])
        for ix, var in mut.iterrows():
            if (var['type'] == conf_row['type'] or \
                (var['type'] in ['INS', 'DEL'] and conf_row['type'] == 'INDEL')):
                try:
                    query_res = tb.query(var.chrom, var.start, var.end)
                except TabixError:
                    query_res = []
                    # Known error, eigen coding has no chrom X, Y data.
                    if not (conf_row['name'] in ['EIGEN_CODING', 'EIGEN_NONCODING'] and var.chrom in ['X', 'Y']):
                        logger.warning('Retriving {} - {} score error for {}:{}-{}'\
                                       .format(conf_row['name'], conf_row['type'],
                                               var.chrom, var.start, var.end))
                for res in query_res:
                    if var.ref == res[conf_row.ref_ix] and var.alt == res[conf_row.alt_ix]:
                        score[ix] = float(res[conf_row.score_ix])
    return score
