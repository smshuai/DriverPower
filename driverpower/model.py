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


def run_glm(X_train, ybinom_train, X_test):
    ''' Run binomial glm in statsmodels
    '''
    X_train = sm.add_constant(X_train, prepend=False)
    X_test  = sm.add_constant(X_test, prepend=False)
    glm = sm.GLM(ybinom_train, X_train, family=sm.families.Binomial())
    glm_res = glm.fit()
    mu_pred = glm_res.predict(X_test)
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


def model(X_train, ybinom_train, X_test, ybinom_test, gnames, grecur=None, method='glm'):
    support_method = ['glm']
    assert method in support_method, 'Invalid model type. Must be chosen from {}'.format(support_method)
    logger.info('Build the model')
    if method == 'glm':
        mu_pred = run_glm(X_train, ybinom_train, X_test)
    res = raw_test(mu_pred, ybinom_test, gnames, grecur)
    return res
