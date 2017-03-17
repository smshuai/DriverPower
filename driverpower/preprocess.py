"""
Pre-process data for DriverPower.
"""
import pandas as pd
import numpy as np
import sys
import logging
from sklearn.preprocessing import StandardScaler, RobustScaler
from pybedtools import BedTool
from driverpower.load import load_mut_bed
from driverpower.detect import getMutCtCg

# create logger
logger = logging.getLogger('PREPROCESS')


def sampling(X, y, recur, N):
    ''' Sample X and y, which are both pd.DF and have the same index.
    If N in (0,1] sample a fraction of data.
    If N is interger > 1, sample N data points.
    '''
    # check index
    assert np.array_equal(X.index, y.index), 'X and y have different indexes'
    assert np.array_equal(recur.index, y.index), 'recur and y have different row indexes'
    if N <= 0:
        logger.error('Sampling value must > 0. You enter {}'.format(N))
        sys.exit(1)
    elif 0 < N < 1:
        logger.info('Sample {}% of data'.format(N*100))
        y = y.sample(frac=N, replace=False)
        y.sort_index(inplace=True)
        recur = recur[recur.index.isin(y.index)]
        recur.sort_index(inplace=True)
        X = X[X.index.isin(y.index)]
        X.sort_index(inplace=True)
    elif N==1:
        # no sampling needed
        pass
    elif N>1:
        N = int(N)
        if N >= y.shape[0]:
            logger.warning('Sampling value >= number of data points {}'.format(y.shape[0]))
            logger.warning('Use all data')
        else:
            logger.info('Sample {} data points'.format(N))
            y = y.sample(n=N, replace=False)
            y.sort_index(inplace=True)
            recur = recur[recur.index.isin(y.index)]
            recur.sort_index(inplace=True)
            X = X[X.index.isin(y.index)]
            X.sort_index(inplace=True)
    assert np.array_equal(X.index, y.index), 'X and y have different row indexes after sampling'
    assert np.array_equal(y.index, recur.index), 'recur and y have different row indexes after sampling'
    return X, y, recur


def get_response(ct, cg):
    ''' Obtain response from count and coverage table
    '''
    # count per bin
    ct_byb = ct.pivot_table(values='ct', index='binID', aggfunc=sum)
    ct_byb.sort_index(inplace=True)
    # coverage per bin
    cg_sum = cg.sum(axis=1).copy()
    cg_sum.name = 'cg'
    # ct.index should be a subset of cg.index
    assert np.sum(~ct.binID.isin(cg.index)) == 0, 'Count table binID is not a subset of Coverage table binID'
    # concat ct per bin and cg per bin
    ct_byb = pd.concat([ct_byb, cg_sum], axis=1)
    # Since ct binID is only a subset of cg binID, fill NA with 0.
    ct_byb.fillna(0, inplace=True)
    assert np.array_equal(ct_byb.index, cg.index), 'ct index does not match cg index'
    # binom GLM responses (# success, # failure)
    ybinom = np.array([ct_byb.ct, ct_byb.cg - ct_byb.ct]).T
    return ybinom


def get_filter(ct, cg, len_threshold=500, recur_threshold=2, return_recur=False, return_tab=False):
    ''' Obtain filter based on length and recurrence.
    Args:
        ct - pd.DF. Count table, 4 cols, no index
        cg - pd.DF. Coverage table, indexed by binID
        len_threshold   - int. Bins with effective length < len_threshold will be discarded
        recur_threshold - int. Bins mutated in < recur_threshold samples will be discarded
        return_recur    - bool. Whether or not return the recur
        return_tab      - bool. Whether or not return the filter table.
    Return:
        keep  - np.array. Index of bins passed all filters
        recur - pd.Series. Indexed by binID. Recurrence of all bins.
        filter_tab - pd.DF. Indexed by binID. Two columns, cg and recur
    '''
    # ct pivot by binID and sid
    ct_byb_bys = ct.pivot_table(values='ct', index='binID', columns='sid', aggfunc=sum).fillna(0)
    # Sum of cg per binID
    cg_sum = cg.sum(axis=1).copy()
    nbin = cg_sum.shape[0]
    logger.info('{} bins before filter'.format(nbin))
    # Recurrence table
    recur = (ct_byb_bys > 0).sum(axis=1)
    filter_tab = pd.concat([cg_sum, recur], axis=1)
    filter_tab.columns = ['cg', 'recur']
    filter_tab.fillna(0, inplace=True) # recur can be 0
    assert np.array_equal(filter_tab.index, cg.index), "Filter index does not match CG index"
    # 1. at least mutated in recur_threshold samples
    keep1 = np.where(filter_tab.recur >= recur_threshold)[0] # index of bins pass
    logger.info('{} ({:.2f}%) bins have mutations in at least {} samples'.format(keep1.shape[0], keep1.shape[0]/nbin*100, recur_threshold))
    # 2. at least len_threshold bp
    keep2 = np.where(filter_tab.cg >= len_threshold)[0] # index of bins pass
    logger.info('{} ({:.2f}%) bins have effective length >= {} bp'.format(keep2.shape[0], keep2.shape[0]/nbin*100, len_threshold))
    # intersection of keep1 and keep2
    keep = np.intersect1d(keep1, keep2)
    logger.info('{} ({:.2f}%) bins passed all filters'.format(keep.shape[0], keep.shape[0]/nbin*100))
    if return_recur:
        return keep, filter_tab.recur
    if return_tab:
        return keep, filter_tab
    return keep # index


def apply_filter(keep, xy_list):
    ''' Apply filter from get_filter() to list of tables
    Args:
        keep - return value from get_filter. Row indexes that pass filters
        xy_list - list of numpy ndarray. Could be 1D or 2D. Only rows will be selected.
    Return:
        xy_filtered - list of numpy ndarray.
    '''
    xy_filtered = []
    for myarr in xy_list:
        xy_filtered.append(myarr[keep])
    return xy_filtered


def filter(ct_train, ct_test, cg_train, cg_test,
    X_train, X_test, ybinom_train, ybinom_test, len_threshold=500, recur_threshold=2):
    # obtain filter
    logger.info('Filtering test data')
    keep_test, grecur  = get_filter(ct_test, cg_test, len_threshold, recur_threshold, True)
    logger.info('Filtering training data')
    keep_train = get_filter(ct_train, cg_train, len_threshold, recur_threshold)
    # apply filter
    X_test,  ybinom_test  = apply_filter(keep_test,  [X_test,  ybinom_test])
    X_train, ybinom_train = apply_filter(keep_train, [X_train, ybinom_train])
    gnames, grecur = apply_filter(keep_test, [cg_test.index.values, np.array(grecur)])
    return (X_test, ybinom_test, X_train, ybinom_train, gnames, grecur)


def scaling(Xtrain, Xtest=None, scaler_type='robust'):
    ''' Scale data with sklearn robust or standard scaler.
    Args:
        Xtrain, Xtest - numpy ndarray, each column represents one feature.
        scaler_type - str, either 'robust' or 'standard'
    Return:
        scaled np.ndarray of Xtrain and Xtest
    '''
    assert scaler_type in ['standard', 'robust', 'none'], 'Scaler must be chosen from [standard, robust, none]'
    if scaler_type == 'none':
        logger.info('Skip data scaling')
        return None
    logger.info('Scaling data with {} scaler'.format(scaler_type))
    if scaler_type == 'robust':
        scaler = RobustScaler()
    elif scaler_type == 'standard':
        scaler = StandardScaler()
    # same scaler fitted with Xtrain
    scaler.fit(Xtrain)
    if Xtest is None: # no Xtest
        return scaler.transform(Xtrain)
    return scaler.transform(Xtrain), scaler.transform(Xtest)


def get_gmean(y, recur):
    ''' Use binomial response y and recur to produce a new gmean response.
    '''
    logger.info('Use geometric mean as response')
    ynew = np.zeros(y.shape, dtype=int)
    ynew[:,0] = np.sqrt(recur * y[:,0])
    ynew[:,1] = y.sum(1) - ynew[:,0]
    return ynew

def preprocess(cg_test, ct_test, X_test, cg_train,
    ct_train, X_train, len_threshold, recur_threshold, scaler_type):
    ''' Main wrapper function for preprocess.
    '''
    # Generate binom responses
    ybinom_train = get_response(ct_train, cg_train)
    ybinom_test  = get_response(ct_test, cg_test)
    # filter
    X_test, ybinom_test, X_train, ybinom_train, gnames, grecur = filter(ct_train, ct_test, cg_train, cg_test,
    X_train, X_test, ybinom_train, ybinom_test, len_threshold, recur_threshold)
    # Scale
    if scaler_type != 'none':
        X_train, X_test = scaling(X_train, X_test, scaler_type)
    return (X_train, ybinom_train, X_test, ybinom_test, gnames, grecur)

##
# v0.5.0
##


def preprocess_v0(args):
    ''' Deprecated. Used in PCAWG Freeze.
    In the hdf5:
        X
        y
        recur
        sid
    '''
    # initial output HDF5
    store = pd.HDFStore(args.out, mode='w')
    ct, cg, cv, recur = load_memsave(args.path_ct,
        args.path_cg, args.path_cv,
        args.len_threshold, args.recur_threshold)
    if cv.shape[0] == 0:
        logger.warning('No bin left')
    # sample IDs
    sid = pd.Series(ct.sid.unique())
    sid.name = 'sid'
    # get response
    ybinom = get_response(ct, cg)
    # y to pd.DF
    ybinom = pd.DataFrame(ybinom, columns=['ct','len_ct'], index=cg.index)
    # sampling
    cv, ybinom, recur = sampling(cv, ybinom, recur, args.sampling)
    # write to store
    store.append('X', cv, chunksize=50000)
    store['y'] = ybinom
    store['recur'] = recur
    store['sid'] = sid
    store.close()


def preprocess_v1(mut_path, callable_path, bin_path, feature_path, out_path):
    ''' Preprocess version 1.1
    Args:
        mut_path - path to variants
        callable_path - path to whitelist regions
        bin_path - path to element set
        feature_path - path to features
        out_path - output h5 path
    '''
    mut_df, ndonor = load_mut_bed(mut_path)
    mut_bed = BedTool.from_dataframe(mut_df, na_rep='NA')
    bin_bed = BedTool(bin_path)
    if callable_path:
        callable_bed = BedTool(callable_path)
        mut_bed = mut_bed.intersect(callable_bed, wa=True)
        ncall = mut_bed.count()
        bin_bed = bin_bed.intersect(callable_bed)
        logger.info('{} ({:.2f}%) mutations are in callable regions'\
                    .format(ncall, ncall/mut_df.shape[0]*100))
    # get nMut, nSample and Length
    mut_cnames = list(mut_df.columns.values)
    bed_cnames = ['chrom_bin', 'start_bin', 'end_bin', 'binID']
    mut_df, ct, cg, recur = getMutCtCg(mut_bed, bin_bed, mut_cnames, bed_cnames)
    # y
    binIDs = pd.read_table(bin_path, sep='\t', header=None,\
                           names=bed_cnames, usecols=['binID'])
    binIDs = binIDs.binID.unique()
    y = pd.DataFrame(index=binIDs, columns=['length', 'nMut', 'nSample'])
    y.index.names = ['binID']
    # add to y
    y['length'] = cg
    y['nSample'] = recur
    y['nMut'] = ct
    # fill na with 0
    y['length'] = y['length'].fillna(0).astype(int)
    y['nSample'] = y['nSample'].fillna(0).astype(int)
    y['nMut'] = y['nMut'].fillna(0).astype(int)
    # default keep bins with nMut > 0 and Length > 100bp
    keep = np.logical_and(y['nMut']>0, y['length']>100)
    y = y[keep]
    y.sort_index(inplace=True)
    keep_bin = y.index.values
    # load X by chunks
    logger.info('Start to load and filter features')
    cv_reader = pd.read_table(feature_path, index_col='binID', chunksize=50000)
    cv = cv_reader.get_chunk()
    cv = cv[cv.index.isin(keep_bin)] # cv container
    Nchunk = int(binIDs.shape[0] / 50000)
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
    assert np.array_equal(cv.index, y.index), 'binIDs in X and y tables are not matched'
    # output
    store = pd.HDFStore(out_path, mode='w')
    store['meta'] = pd.Series({'version': 'v1', 'N':1024})
    store['y'] = y
    store.append('X', cv, chunksize=50000)
    store.close()