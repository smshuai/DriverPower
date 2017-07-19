""" Background mutation rate model

Two types of BMR model are supported:
1. Randomized lasso + GLM
2. Gradient boosting machine

"""

import logging
import sys
import numpy as np
import xgboost as xgb
import statsmodels.api as sm
from sklearn.preprocessing import RobustScaler
from sklearn.linear_model import LassoCV, RandomizedLasso
from sklearn.model_selection import KFold
from scipy.special import logit
from driverpower.dataIO import read_feature, read_response, read_fi, read_param
from driverpower.dataIO import save_scaler, save_fi, save_glm, save_gbm, save_model_info

logger = logging.getLogger('BMR')


def run_bmr(model_name, X_path, y_path,
            fi_cut=0.5, fi_path=None,
            kfold=3, param_path=None,
            project_name='DriverPower', out_dir='./DriverPower.output/'):
    """ Wrapper function for BMR model.

    Args:
        model_name (str): 'GLM' or 'GBM'
        X_path (str): path to training X
        y_path (str): path to training y
        fi_cut (float): cutoff for feature importance
        fi_path (str): path to the feature importance
        kfold (int): K fold CV for GBM
        project_name (str): name of the project
        out_dir (str): directory for saving output files

    Returns:

    """
    use_features = read_fi(fi_path, fi_cut)
    run_feature_select = False if use_features else True
    X = read_feature(X_path, use_features)
    y = read_response(y_path)
    # use bins with both X and y
    use_bins = np.intersect1d(X.index.values, y.index.values)
    X = X.loc[use_bins, :]
    y = y.loc[use_bins, :]
    feature_names = X.columns.values
    # convert X to np.array
    X = X.values
    if model_name == 'GLM':
        # Scale data is necessary for GLM
        X, scaler = scale_data(X)
        save_scaler(scaler, project_name, out_dir)
        if run_feature_select:
            # Run lasso to get alpha
            alpha = run_lasso(X, y)
            # Run rnd lasso to get feature importance
            fi_scores = run_rndlasso(X, y, alpha)
            fi = save_fi(fi_scores, feature_names, project_name, out_dir)
            # Remove unimportant features
            use_features = (fi.importance >= fi_cut).values
            X = X[:, feature_names.isin(use_features)]
        # Run GLM to get trained model
        model = run_glm(X, y)
        save_glm(model, project_name, out_dir)
        # Run dispersion test
        pval, theta = dispersion_test(y.nMut, model.fittedvalues)
        # Save model info.
        model_info = {'model_name': model_name,
                      'pval_dispersion': pval,
                      'theta': theta,
                      'feature_names': feature_names,
                      'use_features': use_features,
                      'project_name': project_name,
                      'model_dir': out_dir}
    elif model_name == 'GBM':
        # make xgb data
        X = xgb.DMatrix(data=X, label=y.nMut, feature_names=feature_names)
        X.set_base_margin(np.log(y.length+1/y.N) + np.log(y.N))
        # k-fold CV
        ks = KFold(n_splits=kfold)
        param = read_param(param_path)
        y['nPred'] = np.nan
        k = 1  # idx of model fold
        for train, valid in ks.split(range(X.num_row())):
            # train the model
            model = run_gbm(X.slice(train), X.slice(valid), param)
            # predict on valid
            y.iloc[valid, 'nPred'] = model.predict(X.slice(valid))
            save_gbm(model, k, project_name, out_dir)
            k += 1
        # Run dispersion test
        pval, theta = dispersion_test(y.nMut, y.nPred)
        model_info = {'model_name': model_name,
                      'pval_dispersion': pval,
                      'theta': theta,
                      'kfold': kfold,
                      'feature_names': feature_names}
    else:
        logger.error('Unknown background model: {}. Please use GLM or GBM'.format(model_name))
        sys.exit(1)
    save_model_info(model_info, project_name, out_dir)


def scale_data(X, scaler=None):
    """ Scale X with robust scaling.
    
    Args:
        X (pd.df): feature matrix indexed by binID.
        scaler (RobustScaler): pre-trained scaler. Default is None
        
    Returns:
        np.array: normalized feature matrix.
        RobustScaler: robust scaler fitted with training data,
            only returned when there is no pre-trained scaler.
    
    """
    if scaler:
        return scaler.transform(X)
    else:
        scaler = RobustScaler(copy=False)
        scaler.fit(X)
        return scaler.transform(X), scaler


def run_lasso(X, y, max_iter=3000, cv=5, n_threads=3):
    """ Implement LassoCV in sklearn
    
    Args:
        X (np.array): scaled X.
        y (pd.df): four columns response table. 
        max_iter (int): max iteration. 
        cv (int): CV fold.
        n_threads (int): Number of threads to use for parallel computing.

    Returns:
        float: trained alpha value.

    """
    logger.info('Implementing LassoCV with {} iter. and {}-fold CV'.format(max_iter, cv))
    # generate logit response
    y_logit = logit((y.nMut + 0.5) / (y.length * y.N))
    # sub-sampling X and y (300,000)
    use_ix = np.random.choice(y_logit.shape[0], 300000, replace=False)
    Xsub = X[use_ix, :]
    ysub = y_logit[use_ix]
    reg = LassoCV(max_iter=max_iter, cv=cv, copy_X=False, n_jobs=n_threads)
    lassocv = reg.fit(Xsub, ysub)
    logger.info('LassoCV alpha = {}'.format(lassocv.alpha_))
    return lassocv.alpha_


def run_rndlasso(X, y, alpha,
    n_resampling=500, sample_fraction=0.1, n_threads=1):
    """  Implement Randomized Lasso in sklearn

    Args:
        X (np.array): scaled X. 
        y (pd.df): four columns response table. 
        alpha (float): parameter trained from lassoCV 
        n_resampling (int): number of times for resampling 
        sample_fraction (float): fraction of data to use at each resampling

    Returns:
        np.array: feature importance scores

    """
    logger.info('Implementing Randomized Lasso with alpha={}, n_resampling={} and sample_fraction={}'.
                format(alpha, n_resampling, sample_fraction))
    # generate logit response
    y_logit = logit((y.nMut + 0.5) / (y.length * y.N))
    reg = RandomizedLasso(alpha=alpha,
                          n_resampling=n_resampling,
                          sample_fraction=sample_fraction,
                          selection_threshold=1e-3,
                          max_iter=3000,
                          normalize=False,
                          n_jobs=n_threads)
    rndlasso = reg.fit(X, y_logit)
    fi_scores = rndlasso.scores_
    return fi_scores


def run_glm(X, y):
    """ Train the binomial GLM
    
    Args:
        X (np.array): scaled X. 
        y (pd.df): four columns response table.

    Returns:
        sm.model: trained GLM models.
        
    """
    # make two columns response (# success, # failure)
    y_binom = np.zeros((y.shape[0], 2), dtype=np.int_)
    y_binom[:,0] = y.nMut
    y_binom[:,1] = y.length * y.N - y.nMut
    # Add const manually. sm.add_constant cannot add 1 for shape (1, n)
    X = np.c_[X, np.ones(X.shape[0])]
    glm = sm.GLM(y_binom, X, family=sm.families.Binomial())
    model = glm.fit()
    return model


def run_gbm(dtrain, dvalid, param):
    # specify validations set to watch performance
    watchlist = [(dvalid, 'eval'), (dtrain, 'train')]
    bst = xgb.train(params=param,
                    dtrain=dtrain,
                    num_boost_round=5000,
                    evals=watchlist,
                    early_stopping_rounds=5)
    return bst


def dispersion_test(pred, obs):
    pval = 0
    theta = 1
    return pval, theta
