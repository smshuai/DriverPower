''' Model module for DriverPower
'''

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import binom_test
from driverpower.preprocess import get_gmean
import logging


# create logger
logger = logging.getLogger('MODEL')


def split_by_cg(cg_train, cg_test=None, fold=4):
    ''' Generate spliters that split the data into k-fold equal-size sets by coverage.
    Args:
        cg_train - Vector of coverage in training set
        cg_test  - Vector of coverage in test set
    '''
    # Output format. Each row corresponds to one fold.
    train_spliter = np.zeros((fold, cg_train.shape[0]), dtype=bool)
    # percentiles
    q   = np.linspace(0,100,fold+1)
    cut = np.percentile(cg_train, q)
    # Expand end points. Make sure all data are used.
    cut[0] = -1
    cut[fold] = cut[fold] + 1e11
    # split training set for each fold
    for i in np.arange(fold):
        train_spliter[i, :] = np.logical_and(cg_train >= cut[i], cg_train < cut[i+1])
    # make sure all data points are used once.
    assert np.sum(train_spliter.sum(1)) == cg_train.shape[0]
    # test data
    if cg_test is not None:
        # split test set as well
        test_spliter  = np.zeros((fold, cg_test.shape[0]), dtype=bool)
        for i in np.arange(fold):
            test_spliter[i, :]  = np.logical_and(cg_test >= cut[i], cg_test < cut[i+1])
        assert np.sum(test_spliter.sum(1)) == cg_test.shape[0], \
        '{} != {}'.format(np.sum(test_spliter.sum(1)), cg_test.shape[0])
        return train_spliter, test_spliter
    return train_spliter


def run_glm(X_train, ybinom_train, X_test):
    ''' Run binomial glm in statsmodels
    '''
    # Add const mannully. sm.add_constant cannot add 1 for shape (1, n)
    X_train = np.c_[X_train, np.ones(X_train.shape[0])]
    X_test  = np.c_[X_test, np.ones(X_test.shape[0])]
    glm = sm.GLM(ybinom_train, X_train, family=sm.families.Binomial())
    glm_res = glm.fit()
    mu_pred = glm_res.predict(X_test)
    return mu_pred


def run_glm_fold(X_train, ybinom_train, X_test, cg_train=None, cg_test=None, fold=1):
    if fold == 1:
        mu_pred = run_glm(X_train, ybinom_train, X_test)
    else:
        # output vector
        mu_pred = np.zeros(X_test.shape[0])
        train_spliter, test_spliter = split_by_cg(cg_train, cg_test, fold)
        for i in np.arange(fold): # pred for each fold
            idx_train = np.where(train_spliter[i, :])[0]
            idx_test  = np.where(test_spliter[i, :])[0]
            mu_pred[idx_test] = \
                run_glm(X_train[idx_train,:], ybinom_train[idx_train,:], X_test[idx_test,:])
    return mu_pred


def raw_test(mu_pred, ybinom_test, gnames, grecur=None):
    ''' Perform binomial test with functional adjustment
    '''
    res = pd.DataFrame({'binID': gnames, 'Length': ybinom_test.sum(1).astype(np.int),
        'nMut':ybinom_test[:, 0].astype(np.int), 'Mu': mu_pred})
    logger.info('Calculate Pval')
    # Raw pval. binom(nMut, length, mu_pred)
    res['PvalRaw'] = [binom_test(x, n, p, 'greater') for x, n, p in zip(res.nMut, res.Length, res.Mu)]
    res['QvalRaw'] = multipletests(res['PvalRaw'], method='fdr_bh')[1]
    res.sort_values('QvalRaw', inplace=True)
    # indexed by binID
    res.set_index('binID', inplace=True)
    if grecur is not None:
        res['nSample'] = grecur.astype(np.int)
        # Mean = (nSample + nMut) / 2
        # res['nMutSample'] = (res.nSample + res.nMut)/2
        # Geometric mean
        res['nMutSample'] = np.sqrt(res.nSample * res.nMut)
        res['nMutSample'] = res.nMutSample.astype(np.int)
        # PvalM. binom(nMutSample, length, mu_pred)
        res['PvalM'] = [binom_test(x, n, p, 'greater') for x, n, p in zip(res.nMutSample, res.Length, res.Mu)]
        res['QvalM'] = multipletests(res.PvalM, method='fdr_bh')[1]
        res.sort_values('QvalM', inplace=True)
    return res


def model(X_train, ybinom_train,
    X_test, ybinom_test,
    gnames, grecur=None, brecur=None, use_gmean=False, method='glm', fold=1):
    support_method = ['glm']
    assert method in support_method, 'Invalid model type. Must be chosen from {}'.format(support_method)
    if use_gmean:
        y_train = get_gmean(ybinom_train, brecur)
    else:
        y_train = ybinom_train
    logger.info('Build the model with {}-fold'.format(fold))
    if method == 'glm':
        # obtain cg_train and cg_test
        cg_train = ybinom_train.sum(1)
        cg_test = ybinom_test.sum(1)
        mu_pred = run_glm_fold(X_train, y_train, X_test, cg_train, cg_test, fold)
    res = raw_test(mu_pred, ybinom_test, gnames, grecur)
    return res

###
# For v0.5.0 DETECT
###
def get_model(X_train, y_train, N_train, use_gmean=True, method='glm'):
    '''
    Args:
        X_train - np.array, scaled
        y_train - pd.DF, indexed by binID, three columns ['length', 'nMut', 'nSample']
        N_train - int, number of donors in training set
        use_gmean - bool, use gmean of nMut and nSample as response if True
        method - str, method to use in prediction
    Return:
        model - fitted model that has predict method
    '''
    support_method = ['glm']
    assert method in support_method, 'Invalid model type. Must be chosen from {}'.format(support_method)
    # make two columns response (# success, # failure)
    y_train2d = make2dy(y_train, N_train, use_gmean)
    if method == 'glm':
        # Add const mannully. sm.add_constant cannot add 1 for shape (1, n)
        X_train = np.c_[X_train, np.ones(X_train.shape[0])]
        glm = sm.GLM(y_train2d, X_train, family=sm.families.Binomial())
        model = glm.fit()
    return model


def estimate_bgmr(model, X_test, method='glm'):
    ''' Estimate the background mutation rate
    Args:
        model - fitted model that has predict() method
        X_test - np.array, scaled
        method - str, method to use in prediction
    Return:
        mu_pred - np.array, predicted background mutation rate
    '''
    support_method = ['glm']
    assert method in support_method, 'Invalid model type. Must be chosen from {}'.format(support_method)
    # make two columns response (# success, # failure)
    if method == 'glm':
        # Add const
        X_test  = np.c_[X_test, np.ones(X_test.shape[0])]
        mu_pred = model.predict(X_test)
    return mu_pred


def do_binom_test(y, N, mu, use_gmean, nsample_thresh=1, nmut_thresh=1, len_thresh=100):
    ''' Perform binomial test and BH correction
    Args:
        y  - pd.DF, indexed by binID, three columns ['length', 'nMut', 'nSample']
        N  - int, number of donors
        mu - pd.Series, vector of mutation rate
        use_gmean - bool, use gmean of nMut and nSample as response if True
        nsample_thresh - int, keep pvals for nsample >= nsample_thresh. Default: 1
        nmut_thresh - int, keep pvals for nmut >= nmut_thresh. Default: 1
        len_thresh - int, keep pvals for length >= len_thresh. Default: 1
    Return:
        pvals - np.array, a list of p-values, one for each element, with NA masked
        qvals - np.array, BH q-values, with NA masked
    '''
    y2d = make2dy(y, N, use_gmean)
    pvals = np.array([binom_test(x, n, p, 'greater') if p<1 else 1 for x, n, p in zip(y2d[:,0], y2d.sum(1), mu)])
    # filter
    keep = np.logical_and(y['length']>=len_thresh, y['nSample']>=nsample_thresh)
    keep = np.array(np.logical_and(keep, y['nMut']>=nmut_thresh))
    remove = np.logical_not(keep)
    pvals[remove] = np.nan
    qvals = pvals.copy()
    qvals[keep] = multipletests(pvals[keep], method='fdr_bh')[1]
    return pvals, qvals

def make2dy(y, N, use_gmean=True):
    ''' Make two columns response (# success, # failure)
    Args:
        y - pd.DF, indexed by binID, three columns ['length', 'nMut', 'nSample']
        N - int, number of donors
        use_geman - bool, use gmean of nMut and nSample as response if True
    Return:
        y2d - np.array, two columns
    '''
    y2d = np.zeros((y.shape[0], 2), dtype=np.int_)
    y2d[:,0] = np.sqrt(y.nMut * y.nSample) if use_gmean else y.nMut
    y2d[:,1] = y['length'] * N - y2d[:,0]
    return y2d
