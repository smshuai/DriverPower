''' Feature selection module for DriverPower
'''
import logging
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LassoCV, RandomizedLasso
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.feature_selection import f_regression
from scipy.special import logit
from scipy.stats import spearmanr
from driverpower.load import load_hdf5
from driverpower.model import make2dy


# create logger
logger = logging.getLogger('FEATURE SELECT')


def run_spearmanr(X, ybinom):
    ''' Implement univariate spearmanr
    '''
    ylogit = logit(ybinom[:,0]/ybinom.sum(1))
    rho = np.apply_along_axis(lambda x: spearmanr(x, ylogit)[0], 0, X)
    return rho

def run_fregression(X, ybinom):
    ''' Implement uniariate F regression
    '''
    ylogit = logit(ybinom[:,0]/ybinom.sum(1))
    freg = f_regression(X, ylogit, center=False)
    return freg[0]

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
        selection_threshold=1e-3, max_iter=3000, normalize=False)
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
    if len(fset) == 0:
        logger.error('No feature is selected. Please try another cutoff')
        sys.exit(1)
    return fset, fidx


def fselect(X_train, X_test, ybinom_train, fnames, method='rndlasso', cutoff_rndlasso=0.5, cutoff_lasso=0.001):
    ''' Main wrapper function for feature selection
    '''
    support_method = ['rndlasso', 'lasso']
    assert method in support_method, 'Invalid feature selection method. Must be chosen from {}'.format(support_method)
    if method == 'rndlasso':
        # find alpha
        lassocv  = run_lasso(X_train, ybinom_train, max_iter=3000, cv=5)
        # run rndlasso
        rndlasso = run_rndlasso(X_train, ybinom_train, alpha=lassocv.alpha_)
        # feature importance
        fscores = rndlasso.scores_
        cutoff = cutoff_rndlasso
    if method == 'lasso':
        lassocv = run_lasso(X_train, ybinom_train, max_iter=3000, cv=10)
        fscores = np.abs(lassocv.coef_)
        cutoff = cutoff_lasso
    fset, fidx = feature_score(fscores, fnames, cutoff)
    logger.info('At cutoff={}, {} selected features are: {}'.format(cutoff, len(fset), ", ".join(fset)))
    return X_train[:, fidx], X_test[:, fidx], fscores

###
# v0.5.0
###
def fselect_v1(h5_path, scaler_type, use_gmean, out_path):
    ''' Run feature selection for preprocess HDF5 (v1)
    '''
    logger.info('Loading training data ...')
    Xtrain_df, ytrain_df, Ntrain = load_hdf5(trainH5_path, usefeatures)
    train_columns = Xtrain_df.columns.values
    Xtrain_mat = Xtrain_df.as_matrix()
    # scaling if need
    if scaler_type == 'robust':
        logger.info('Use robust scaler')
        scaler = RobustScaler()
        scaler.fit(Xtrain_mat)
        Xtrain_mat = scaler.transform(Xtrain_mat)
    elif scaler_type == 'standard':
        logger.info('Use standard scaler')
        scaler = StandardScaler()
        scaler.fit(Xtrain_mat)
        Xtrain_mat = scaler.transform(Xtrain_mat)
    ytrain2d = make2dy(ytrain_df, Ntrain, use_gmean)
    lassocv = run_lasso(Xtrain_mat, ytrain2d)
    rndlasso = run_rndlasso(Xtrain_mat, ytrain2d, alpha=lassocv.alpha_)
    fscores = rndlasso.scores_
    res = pd.DataFrame(fscores, index=train_columns, columns=['rndlasso'])
    res.index.name = 'fname'
    res.to_csv(out_path, sep='\t')
