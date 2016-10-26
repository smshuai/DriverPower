import pandas as pd
import numpy as np
import scipy as sp

from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.linear_model import LassoCV, RandomizedLasso

from statsmodels.sandbox.stats.multicomp import multipletests
import statsmodels.api as sm
import warnings

import sys
import os



# helper functions
def load_coverage(path_coverage):
    ''' Load covergae table. Return a pd.DF with binID as the key
    '''
    cg = pd.read_csv(path_coverage, sep='\t', header=0)
    for col in ['chrom', 'start', 'end']:
        if col in cg.columns.values:
            del cg[col]
    cg = cg.set_index('binID')
    cg.sort_index(inplace=True)
    assert cg.shape[1] == 32, "Number of columns (3-mer) for coverage table is not 32." # ncol = 32
    return cg

def load_covar(path_covar):
    ''' Load pre-processed covar. table. Return a pd.DF with binID as the key
    '''
    cv = pd.read_csv(path_covar, sep='\t', header=0)
    cv = cv.set_index('binID')
    cv.sort_index(inplace=True)
    na_count = cv.isnull().sum()
    print(na_count[na_count>0])
    cv.fillna(0, inplace=True)
    return cv

def load_count(path_count):
    ''' Load mutatin count data. 4 cols TSV file ['binID', 'sid', 'categ', 'ct']
    '''
    ct = pd.read_csv(path_count, sep='\t', header=0)
    ct.sort_values('binID', inplace=True)
    return ct

def load_all(path_cg, path_ct, path_cv):
    cg = load_coverage(path_cg)
    print('- CG: ', cg.shape)
    ct = load_count(path_ct)
    print('- CT: ', ct.shape)
    cv = load_covar(path_cv)
    print('- CV: ', cv.shape)
    assert np.array_equal(cv.index, cg.index), 'CV and CG binIDs do not match'
    return (cg, ct, cv)

def get_response(ct, cg):
    # Response - Test data
    ct_byb = ct.pivot_table(values='ct', index='binID', aggfunc=sum) # count by binID
    ct_byb.sort_index(inplace=True)
    cg_sum = cg.sum(axis=1).copy()
    cg_sum.name = 'cg'
    ct_byb = pd.concat([ct_byb, cg_sum], axis=1)
    ct_byb.fillna(0, inplace=True)
    assert np.array_equal(ct_byb.index, cg.index), 'ct index does not match cg index'
    # binom GLM responses (# success, # failure)
    ybinom = np.array([ct_byb.ct, ct_byb.cg - ct_byb.ct]).T
    print('- Binom Response shape: ', ybinom.shape)
    return ybinom

def get_filter(ct, cg, len_threshold=500, recur_threshold=2, return_recur=False):
    # ct by binID and by sid
    ct_byb_bys = ct.pivot_table(values='ct', index=['binID'], columns=['sid'], aggfunc=sum)
    ct_byb_bys.fillna(0, inplace=True)
    # Sum of cg per binID
    cg_sum = cg.sum(axis=1).copy()
    print('{} bins in input'.format(cg_sum.shape[0]))
    # Recurrence table
    recur = (ct_byb_bys > 0).sum(axis=1)
    filter_tab = pd.concat([cg_sum, recur], axis=1)
    filter_tab.columns = ['cg', 'recur']
    filter_tab.fillna(0, inplace=True) # recur can be 0
    assert np.array_equal(filter_tab.index, cg.index), "Filter index does not match CG index"
    # 1. at least len_threshold bp
    keep1 = np.where(filter_tab.cg >= len_threshold)[0] # index of bins pass
    print('- {} bins >= {} bp'.format(keep1.shape[0], len_threshold))
    # 2. at least mutated in recur_threshold samples
    keep2 = np.where(filter_tab.recur >= recur_threshold)[0] # index of bins pass
    print('- {} bins have mutations in at least {} samples'.format(keep2.shape[0], recur_threshold))
    # intersection of keep1 and keep2
    keep = np.intersect1d(keep1, keep2)
    print('- {} ({:.2f}%) bins passed all filters'.format(keep.shape[0], keep.shape[0]/cg_sum.shape[0] * 100))
    if return_recur:
        return keep, filter_tab.recur
    else:
        return keep

def apply_filter(keep, xy_list):
    ''' Apply filter from get_filter() to X and Y list
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


def split_by_cg(cg_train, cg_test, fold=4):
    '''
    '''
    # Output format. Each row corresponds to one fold.
    train_spliter = np.zeros((fold, cg_train.shape[0]), dtype=bool)
    test_spliter  = np.zeros((fold, cg_test.shape[0]), dtype=bool)
    q   = np.linspace(0,100,fold+1)
    cut = np.percentile(cg_train, q)
    cut[0] = cut[0] - 1
    cut[fold] = cut[fold] + 1
    for i in np.arange(fold):
        train_spliter[i, :] = np.logical_and(cg_train >= cut[i], cg_train < cut[i+1])
        test_spliter[i, :]  = np.logical_and(cg_test >= cut[i], cg_test < cut[i+1])
    assert np.sum(train_spliter.sum(1)) == cg_train.shape[0]
    assert np.sum(test_spliter.sum(1)) == cg_test.shape[0]
    return train_spliter, test_spliter
