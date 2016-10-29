''' Feature selection module for DriverPower
'''
import logging
import numpy as np
from sklearn.linear_model import LassoCV, RandomizedLasso
from scipy.special import logit


# create logger
logger = logging.getLogger('FEATURE SELECT')
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


def run_lasso(X_train, ybinom_train, max_iter=3000, cv=5):
    ''' Implement LassoCV provided by sklearn
    Note:
        To generate alpha for rnalasso, use max_iter=3000 and cv=5
        For feature selection purpose, use max_iter=3000 and cv=10
    '''
    logger.info('Implementing LassoCV with {} iter. and {}-fold CV'.format(max_iter, cv))
    # generate logit response
    ylogit_train = logit(ybinom_train[:,0]/ybinom_train.sum(1))
    clf = LassoCV(max_iter=max_iter, cv=cv)
    lassocv = clf.fit(X_train, ylogit_train)
    logger.info('LassoCV alpha = {}'.format(lassocv.alpha_))
    return lassocv


def run_rndlasso(X_train, ybinom_train, alpha,
    n_resampling=200, sample_fraction=0.4):
    ''' Implement RandomizedLasso provided by sklearn
    '''
    logger.info('Implementing Randomized Lasso with alpha={}, n_resampling={} and sample_fraction={}'.format(alpha, n_resampling, sample_fraction))
    # generate logit response
    ylogit_train = logit(ybinom_train[:,0]/ybinom_train.sum(1))
    clf = RandomizedLasso(alpha=alpha, n_resampling=n_resampling,
        sample_fraction=sample_fraction,
        selection_threshold=1e-3, normalize=False)
    rndlasso = clf.fit(X_train, ylogit_train)
    return rndlasso

def feature_score(fscores, fnames, cutoff):
    ''' Select features based on scores and cutoff.
    Args:
        fscores - np.array. Feature importance scores
        fnames  - np.array. Feature names (same length as fscores)
        cutoff  - cutoff for feature selection
    Return:
        fset - np.array. Names of selected features
        fidx - np.array. Numerical index of selected features
    '''
    assert len(fnames) == len(fscores), 'Feature importance scores and names have diffent length'
    fidx = np.where(fscores > cutoff)[0]
    fset = fnames[fidx]
    return fset, fidx


def fselect(X_train, X_test, ybinom_train, fnames, method='rndlasso', cutoff_rndlasso=0.5, cutoff_lasso=0.001):
    ''' Main wrapper function for feature selection
    '''
    support_method = ['rndlasso', 'lassocv']
    assert method in support_method, 'Invalid feature selection method. Must be chosen from {}'.format(support_method)
    if method == 'rndlasso':
        # find alpha
        lassocv  = run_lasso(X_train, ybinom_train, max_iter=3000, cv=5)
        # run rndlasso
        rndlasso = run_rndlasso(X_train, ybinom_train, alpha=lassocv.alpha_)
        # feature importance
        fscores = rndlasso.scores_
        cutoff = cutoff_rndlasso
    if method == 'lassocv':
        lassocv = run_lasso(X_train, ybinom_train, max_iter=3000, cv=10)
        fscores = np.abs(lassocv.coef_)
        cutoff = cutoff_lasso
    fset, fidx = feature_score(fscores, fnames, cutoff)
    logger.info('At cutoff={}, {} selected features are: {}'.format(cutoff, len(fset), ", ".join(fset)))
    return X_train[:, fidx], X_test[:, fidx], fscores
