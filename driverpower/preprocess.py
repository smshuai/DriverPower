"""
Pre-process data for DriverPower.

Pre-process steps are:
1. Calculate response
2. Filter
3. Scale 
"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, RobustScaler
import logging


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

    # concat ct per bin and cg per bin
    ct_byb = pd.concat([ct_byb, cg_sum], axis=1)
    # Since ct binID is only a subset of cg binID, fill NA with 0.
    ct_byb.fillna(0, inplace=True)
    assert np.array_equal(ct_byb.index, cg.index), 'ct index does not match cg index'
    # binom GLM responses (# success, # failure)
    ybinom = np.array([ct_byb.ct, ct_byb.cg - ct_byb.ct]).T
    print('- Binom Response shape: ', ybinom.shape)
    return ybinom


def get_filter(ct, cg, len_threshold=500, recur_threshold=2, return_recur=False):
    ''' Obtain filter based on length and recurrence.
    Args:
        ct - pd.DF. Count table, 4 cols, no index
        cg - pd.DF. Coverage table, indexed by binID
        len_threshold   - int. Bins with effective length < len_threshold will be discarded
        recur_threshold - int. Bins mutated in < recur_threshold samples will be discarded
        return_recur    - bool. Whether or not return the recur
    Return:
        keep  - np.array. Index of bins passed all filters
        recur - np.Series. Indexed by binID. Recurrence of all bins.
    '''
    # ct pivot by binID and sid
    ct_byb_bys = ct.pivot_table(values='ct', index='binID', columns='sid', aggfunc=sum).fillna(0)
    # Sum of cg per binID
    cg_sum = cg.sum(axis=1).copy()
    nbin = cg_sum.shape[0]
    logging.info('{} bins before filter'.format(nbin))
    # Recurrence table
    recur = (ct_byb_bys > 0).sum(axis=1)
    filter_tab = pd.concat([cg_sum, recur], axis=1)
    filter_tab.columns = ['cg', 'recur']
    filter_tab.fillna(0, inplace=True) # recur can be 0
    assert np.array_equal(filter_tab.index, cg.index), "Filter index does not match CG index"
    # 1. at least mutated in recur_threshold samples
    keep1 = np.where(filter_tab.recur >= recur_threshold)[0] # index of bins pass
    logging.info('{} ({:.2f}%) bins have mutations in at least {} samples'.format(keep1.shape[0], keep1.shape[0]/nbin*100, recur_threshold))
    # 2. at least len_threshold bp
    keep2 = np.where(filter_tab.cg >= len_threshold)[0] # index of bins pass
    logging.info('{} ({:.2f}%) bins have effective length >= {} bp'.format(keep2.shape[0], keep2.shape[0]/nbin*100, len_threshold))
    # intersection of keep1 and keep2
    keep = np.intersect1d(keep1, keep2)
    logging.info('{} ({:.2f}%) bins passed all filters'.format(keep.shape[0], keep.shape[0]/nbin*100))
    if return_recur:
        return keep, filter_tab.recur
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
    logging.info('Filtering test data')
    keep_test, grecur  = get_filter(ct_test, cg_test, len_threshold, recur_threshold, True)
    logging.info('Filtering training data')
    keep_train = get_filter(ct_train, cg_train, len_threshold, recur_threshold)
    # apply filter
    X_test,  ybinom_test  = apply_filter(keep_test,  [X_test,  ybinom_test])
    X_train, ybinom_train = apply_filter(keep_train, [X_train, ybinom_train])
    gnames, grecur = apply_filter(keep_test, [cg_test.index.values, np.array(grecur)])
    return (X_test, ybinom_test, X_train, ybinom_train, gnames, grecur)


def scaling(Xtrain, Xtest, scaler_type):
    ''' Scale data with sklearn robust or standard scaler.
    Args:
        Xtrain, Xtest - numpy ndarray, each column represents one feature.
        scaler_type - str, either 'robust' or 'standard'
    Return:
        scaled np.ndarray of Xtrain and Xtest
    '''
    assert scaler_type in ['standard', 'robust'], 'Scaler must be chosen from "standard" or "robust"'
    if scaler_type == 'robust':
        scaler = RobustScaler()
    elif scaler_type == 'standard':
        scaler = StandardScaler()
    # same scaler fitted with Xtrain
    scaler.fit(Xtrain)
    return scaler.transform(Xtrain), scaler.transform(Xtest)


def preprocess():
    pass