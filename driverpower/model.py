''' Model module for DriverPower
'''

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import binom_test
import logging


# create logger
logger = logging.getLogger('MODEL')
# logger.setLevel(logging.INFO)
# # create console handler
# ch = logging.StreamHandler()
# ch.setLevel(logging.INFO)
# # create formatter and add it to the handlers
# formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s',
#     datefmt='%m/%d/%Y %H:%M:%S')
# ch.setFormatter(formatter)
# # add the handlers to the logger
# logger.addHandler(ch)


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
    X_train = sm.add_constant(X_train, prepend=False)
    X_test  = sm.add_constant(X_test, prepend=False)
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


def get_gmean(y, recur):
    ''' Use binomial response y and recur to produce a new gmean response.
    '''
    ynew = np.zeros(y.shape, dtype=int)
    ynew[:,0] = np.sqrt(recur * y[:,0])
    ynew[:,1] = y.sum(1) - ynew[:,0]
    return ynew


def model(X_train, ybinom_train, X_test, ybinom_test, gnames,
    grecur=None, brecur=None, use_gmean=False, method='glm', fold=1):
    support_method = ['glm']
    assert method in support_method, 'Invalid model type. Must be chosen from {}'.format(support_method)
    if use_gmean:
        logger.info('Use geometric mean as response')
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
