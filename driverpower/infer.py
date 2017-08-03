""" Making inferences about driver.

This module supports:
1. Predict nMut
2. Burden test
3. Functional adjusted test

"""

import logging
import sys
import os
import numpy as np
from scipy.stats import binom_test, nbinom
from driverpower.dataIO import read_model_info, read_feature, read_response, read_gbm, read_glm, read_scaler, read_fs
from driverpower.dataIO import save_result
from driverpower.model import scale_data
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import xgboost as xgb
    from statsmodels.sandbox.stats.multicomp import multipletests


logger = logging.getLogger('INFER')


def make_inference(model_dir, model_info_path,
                   X_path, y_path,
                   fs_path=None, fs_cut=None,
                   test_method='auto', scale=1, use_gmean=True,
                   project_name='DriverPower', out_dir='./output'):
    """ Main wrapper function for inference

    Args:
        model_path (str): path to the model
        model_info_path (str): path to the model information
        X_path (str): path to the X
        y_path (str): path to the y
        fs_path (str): path to the functional score file
        fs_cut (str): "CADD:0.01;DANN:0.03;EIGEN:0.3"
        test_method (str): 'binomial', 'negative_binomial' or 'auto'.
        scale (float): scaling factor used in negative binomial distribution.
        use_gmean (bool): use geometric mean in test.

    Returns:

    """
    model_info = read_model_info(model_info_path)
    logger.info('Model type: {}'.format(model_info['model_name']))
    model_dir = model_dir if model_dir is not None else model_info['model_dir']
    model_name = model_info['model_name']
    # Load data
    X = read_feature(X_path)
    # order X by feature names of training data
    X = X.loc[:, model_info['feature_names']]
    y = read_response(y_path)
    # use bins with both X and y
    use_bins = np.intersect1d(X.index.values, y.index.values)
    X = X.loc[use_bins, :].values  # X is np.array now
    y = y.loc[use_bins, :]
    # scale X for GLM
    if model_name == 'GLM':
        scaler_path = os.path.join(model_info['model_dir'], model_info['project_name'] + '.feature_scaler.pkl')
        scaler = read_scaler(scaler_path)
        X = scale_data(X, scaler)
        X = X[:, np.isin(model_info['feature_names'], model_info['use_features'])]
    # make prediction
    if model_name == 'GLM':
        y['nPred'] = predict_with_glm(X, y, model_dir, model_info)
    elif model_name == 'GBM':
        y['nPred'] = predict_with_gbm(X, y, model_dir, model_info)
    else:
        logger.error('Unknown background model: {}. Please use GLM or GBM'.format(model_name))
        sys.exit(1)
    # burden test
    count = np.sqrt(y.nMut * y.nSample) if use_gmean else y.nMut
    offset = y.length * y.N + 1
    y['raw_p'] = burden_test(count, y.nPred, offset,
                             test_method, model_info, scale)
    y['raw_q'] = bh_fdr(y.raw_p)
    # functional adjustment
    y = functional_adjustment(y, fs_path, fs_cut, test_method,
                              model_info, scale, use_gmean)
    # save to disk
    save_result(y, project_name, out_dir)
    logger.info('Job done!')


def predict_with_glm(X, y, model_dir, model_info):
    """ Predict number of mutation with GLM.

    Args:
        X (np.array): feature matrix.
        y (pd.df): response.
        model_dir (str): directory to the model pkl.
        model_info (dict): model meta-data.

    Returns:
        np.array: array of predictions.

    """
    assert model_info['model_name'] == 'GLM',\
        'Wrong model name in model info: {}. Need GLM.'.format(model_info['model_name'])
    model_path = os.path.join(model_dir, model_info['project_name']+'.GLM.pkl')
    model = read_glm(model_path)
    # Add const. to X
    X = np.c_[X, np.ones(X.shape[0])]
    pred = np.array(model.predict(X) * y.length * y.N)
    return pred


def predict_with_gbm(X, y, model_dir, model_info):
    """

    Args:
        X:
        y:
        model_dir:
        model_info:

    Returns:

    """
    assert model_info['model_name'] == 'GBM',\
        'Wrong model name in model info: {}. Need GBM.'.format(model_info['model_name'])
    testData = xgb.DMatrix(data=X, label=y.nMut.values, feature_names=model_info['feature_names'])
    testData.set_base_margin(np.array(np.log(y.length+1/y.N) + np.log(y.N)))
    kfold = model_info['kfold']
    pred = np.zeros(y.shape[0])
    for k in range(1, kfold+1):
        model = read_gbm(k, model_info['project_name'], model_dir, model_info['params'])
        pred += model.predict(testData)
    pred = pred / kfold
    return pred


def burden_test(count, pred, offset, test_method, model_info, s):
    """ Perform burden test.

    Args:
        count:
        pred:
        offset:
        test_method:
        model_info:
        s:
        use_gmean:

    Returns:

    """
    if test_method == 'auto':
        test_method = 'binomial' if model_info['pval_dispersion'] > 0.05 else 'negative_binomial'
    if test_method == 'negative_binomial':
        logger.info('Using negative binomial test with s={}, theta={}'.format(s, model_info['theta']))
        theta = s * model_info['theta']
        pvals = np.array([negbinom_test(x, mu, theta)
                          for x, mu in zip(count, pred)])
    elif test_method == 'binomial':
        logger.info('Using binomial test')
        pvals = np.array([binom_test(x, n, p, 'greater')
                          for x, n, p in zip(count, offset,
                                             pred/offset)])
    else:
        logger.error('Unknown test method: {}. Please use binomial, negative_binomial or auto'.format(test_method))
        sys.exit(1)
    return pvals


def negbinom_test(x, mu, theta):
    """ Test with negative binomial distribution

    Convert mu and theta to scipy parameters n and p:

    p = 1 / (theta * mu + 1)
    n = mu * p / (1 - p)

    Args:
        x (float): observed number of mutations (or gmean).
        mu (float): predicted number of mutations (mean of negative binomial distribution).
        theta (float): dispersion parameter of negative binomial distribution.

    Returns:
        float: p-value from NB CDF. pval = 1 - F(n<x)

    """
    p = 1 / (theta * mu + 1)
    n = mu * p / (1 - p)
    pval = 1 - nbinom.cdf(x, n, p, loc=1)
    return pval


def functional_adjustment(y, fs_path, fs_cut, test_method,
                          model_info, scale, use_gmean=True):
    """

    Args:
        y:
        fs_path:
        fs_cut:
        test_method:
        model_info:
        scale:
        use_gmean:

    Returns:

    """
    if fs_path is None:
        return y
    # Convert fs_cut to a dict. "CADD:0.01;DANN:0.03;EIGEN:0.3"
    fs_cut_dict = dict([i.split(':') for i in fs_cut.strip().split(';')])
    fs = read_fs(fs_path, fs_cut_dict)
    # merge with y
    y = y.join(fs)
    ct = 0  # number of score used
    avg_weight = np.zeros(y.shape[0])
    for score, cutoff in fs_cut_dict.items():
        logger.info('Using {} scores'.format(score))
        if float(cutoff) < 0  or float(cutoff) > 1:
            logger.warning('Unused score {} due to invalid cutoff {}, must between 0 and 1'.format(score, cutoff))
            continue
        if float(cutoff) == 0:
            cutoff = 0.001  # add a small number for 0 cutoff
        threshold = -10*np.log10(float(cutoff))  # convert to phred-scale
        # Calculate weight for near-significant elements
        weight = score+'_weight'
        y[weight] = y[score] / threshold
        y.loc[y.raw_q>.25, weight] = 1
        # Calculate nMut * weight for near-significant elements
        n_score = score+'_nMut'
        y[n_score] = y[weight] * y.nMut
        # Calculate p-values (q-values) for near-significant elements
        count = np.sqrt(y.loc[y.raw_q<=.25, n_score] * y.loc[y.raw_q<=.25, 'nSample'])\
            if use_gmean else y.loc[y.raw_q<=.25, n_score]
        offset = y.loc[y.raw_q<=.25, 'length'] * y.loc[y.raw_q<=.25, 'N'] + 1
        p_adj = score+'_p'
        y[p_adj] = y.raw_p
        y.loc[y.raw_q<=.25, p_adj] = burden_test(count, y.loc[y.raw_q<=.25, 'nPred'],
                                                 offset, test_method, model_info, scale)
        y[score+'_q'] = bh_fdr(y[p_adj])
        avg_weight += y[weight]
        ct += 1
    # Use combined weights if more than 2 scores are used
    if ct >= 2:
        y['avg_weight'] = avg_weight / ct
        y['avg_p'] = y.raw_p
        y['avg_nMut'] = y.avg_weight * y.nMut
        count = np.sqrt(y.loc[y.raw_q<=.25, 'avg_nMut'] * y.loc[y.raw_q<=.25, 'nSample'])\
            if use_gmean else y.loc[y.raw_q<=.25, 'avg_nMut']
        offset = y.loc[y.raw_q<=.25, 'length'] * y.loc[y.raw_q<=.25, 'N'] + 1
        y.loc[y.raw_q<=.25, 'avg_p'] = burden_test(count, y.loc[y.raw_q<=.25, 'nPred'],
                                                   offset, test_method, model_info, scale)
        y['avg_q'] = bh_fdr(y.avg_p)
    return y


def bh_fdr(pvals):
    """ BH FDR correction

    Args:
        pvals (np.array): array of p-values

    Returns:
        np.array: array of q-values
    """
    return multipletests(pvals, method='fdr_bh')[1]
