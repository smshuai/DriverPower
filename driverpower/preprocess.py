"""
Pre-process data for DriverPower.
Pre-process steps are:
1. Filter
2. Scale
3. Calculate response
"""
import pandas as pd
import numpy as np

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

def filter(ct_train, ct_test, cg_train, cg_test,
    Xtrain, Xtest, ybinom_train, ybinom_test):
    # obtain filter
    print('Test data')
    keep_test, grecur  = get_filter(ct=ct_test, cg=cg_test, return_recur=True)
    print('Train data')
    keep_train = get_filter(ct=ct_train, cg=cg_train)
    # apply filter
    X_test,  ybinom_test  = apply_filter(keep_test,  [Xtest,  ybinom_test])
    X_train, ybinom_train = apply_filter(keep_train, [Xtrain, ybinom_train])
    gnames, grecur = apply_filter(keep_test, [cg_test.index.values, np.array(grecur)])
    return (X_test, ybinom_test, X_train, ybinom_train, gnames, grecur)
