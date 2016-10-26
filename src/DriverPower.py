#!/usr/bin/env python
import helperDP as dp
import os
import sys
import warnings
from datetime import datetime

import pandas as pd
import numpy as np
import scipy as sp
from scipy.stats import binom_test

import sklearn
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.linear_model import LassoCV, RandomizedLasso

import statsmodels
from statsmodels.sandbox.stats.multicomp import multipletests
import statsmodels.api as sm

''' DriverPower - Find the cancer driver potential for genomic elements

Steps:
1. Read cg, ct, cv for training and test data
2. Scaling
3. Filter
4. Feature selection
5. Modelling
'''

def get_param():
    path_cg_test  = './coverage/gc19.promDomain.bed.cg'
    path_ct_test  = './data/CRC/promDomain/CRC.promDomain.ct.tsv'
    path_cv_test  = './covariates/combined/cv.promDomain.tsv'
    path_cg_train = './coverage/strat.exon.bed.cg'
    path_ct_train = './data/CRC/strat/CRC.strat.ct.tsv'
    path_cv_train = './covariates/combined/cv.strat.tsv'
    scaler_type   = 'robust'
    skip_fs       = True # bool, skip feature selection
    output_prefix = 'CRC_promDomain_robust'
    output_dir    = './output/'
    save_data     = True # bool, save pre-processed data
    preprocessed  = False # whether or not the data is pre-processed (in .npy format)
    scale_first   = False # True: scale before filter; False: filter before scale
    assert scaler_type in ['standard', 'robust'], 'Scaler must be chosen from "standard" or "robust"'

    # print params
    print('Input parameters: ')
    print('''
        path_cg_test  = {}
        path_ct_test  = {}
        path_cv_test  = {}
        path_cg_train = {}
        path_ct_train = {}
        path_cv_train = {}
        scaler_type   = {}
        skip_fs       = {}
        output_prefix = {}
        output_dir    = {}
        save_data     = {}
        preprocessed  = {}
        scale_first   = {}
        '''.format(path_cg_test, path_ct_test, path_cv_test, path_cg_train, path_ct_train,
            path_cv_train, scaler_type, skip_fs, output_prefix, output_dir, save_data, preprocessed, scale_first))
    print('Versions:')
    print('- Python: ', sys.version)
    print('- Statsmodels: ', statsmodels.__version__)
    print('- Sklearn: ', sklearn.__version__)
    print('- Numpy: ', np.__version__)
    print('- Scipy: ', sp.__version__)
    print('- Pandas: ', pd.__version__)
    return (path_cg_test, path_ct_test, path_cv_test, path_cg_train, path_ct_train,
        path_cv_train, scaler_type, skip_fs, output_prefix, output_dir, save_data, preprocessed, scale_first)

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

def feature_selection(Xtrain, ytrain, fs_method, max_iter=2000, cv=5, alpha=0.00459451108436):
    '''
    '''
    if fs_method == 'lassocv':
        print('Implementing LassoCV with {} iter. and {}-fold CV'.format(max_iter, cv))
        clf = LassoCV(max_iter=max_iter, cv=10)
        lassocv = clf.fit(Xtrain, ytrain)
        print('LassoCV alpha = ', lassocv.alpha_)
        return lassocv, abs(lassocv.coef_)
    elif fs_method == 'rndlasso':
        # apply rndlasso
        print('Implementing Random Lasso with alpha={}'.format(alpha))
        clf = RandomizedLasso(alpha=alpha, sample_fraction=0.4, selection_threshold=1e-3, normalize=False)
        rndlasso = clf.fit(Xtrain, ytrain)
        return rndlasso, rndlasso.scores_

def feature_score(fscore, fnames, cutoff):
    ''' Return selected feature names and index vector
    '''
    fidx = np.where(fscore > cutoff)[0]
    fset = fnames[fidx]
    return fset, fidx

def read_raw_data(path_cg_test, path_ct_test, path_cv_test,
    path_cg_train, path_ct_train, path_cv_train):
    ''' Read cg ct and cv data from TSV or CSV files.
    '''
    print('Test Data')
    cg_test,  ct_test,  cv_test  = dp.load_all(path_cg=path_cg_test,  path_ct=path_ct_test,  path_cv=path_cv_test)
    print('Training Data')
    cg_train, ct_train, cv_train = dp.load_all(path_cg=path_cg_train, path_ct=path_ct_train, path_cv=path_cv_train)
    assert np.array_equal(cv_test.columns, cv_train.columns), 'Feature names do not match in Training and Testing data'
    fnames = cv_test.columns.values
    ##
    # Remove strat bins that overlap exons!
    ##
    blackbinID = np.genfromtxt('/u/sshuai/current/annotation/strat.exon.blacklist.txt', dtype=np.str)
    keep_cv = np.logical_not(cv_train.index.isin(blackbinID))
    keep_cg = np.logical_not(cg_train.index.isin(blackbinID))
    keep_ct = np.logical_not(ct_train.binID.isin(blackbinID))
    cv_train = cv_train[keep_cv]
    cg_train = cg_train[keep_cg]
    ct_train = ct_train[keep_ct]
    print('After blacklist train binID:')
    print('CV ', cv_train.shape)
    print('CG ', cg_train.shape)
    print('CT ', ct_train.shape)
    return (cg_test, ct_test, cv_test.as_matrix(),
        cg_train, ct_train, cv_train.as_matrix(), fnames)

def read_npy_data(path_X_train, path_X_test, path_y_train,
    path_y_test, path_gnames, path_fnames, path_grecur):
    ''' Read preprocessed X and y_binom .npy files.
    '''
    X_train = np.load(path_X_train)
    X_test  = np.load(path_X_test)
    ybinom_train = np.load(path_y_train)
    ybinom_test  = np.load(path_y_test)
    grecur = np.loadtxt(path_grecur)
    gnames = np.genfromtxt(path_gnames, dtype=np.str)
    fnames = np.genfromtxt(path_fnames, dtype=np.str)
    return X_train, X_test, ybinom_train, ybinom_test, gnames, fnames, grecur

def save_processed_data(X_train, X_test, ybinom_train, ybinom_test,
    gnames, fnames, grecur, output_dir, output_prefix):
    # Save X
    Xtrain_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'Xtrain.npy']))
    np.save(Xtrain_outpath, X_train)
    Xtest_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'Xtest.npy']))
    np.save(Xtest_outpath, X_test)
    # Save y
    ytrain_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'ytrain.npy']))
    np.save(ytrain_outpath, ybinom_train)
    ytest_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'ytest.npy']))
    np.save(ytest_outpath, ybinom_test)
    # Save gnames, fnames, grecur
    gnames_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'test_bin_names.txt']))
    pd.Series(gnames).to_csv(gnames_outpath, index=False)
    fnames_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'feature_names.txt']))
    pd.Series(fnames).to_csv(fnames_outpath, index=False)
    grecur_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'test_bin_recur.txt']))
    pd.Series(grecur).to_csv(grecur_outpath, index=False)
    print('\n'.join([Xtrain_outpath, Xtest_outpath, ytrain_outpath, ytest_outpath, gnames_outpath, fnames_outpath, grecur_outpath]))
    return 0

def filter(ct_train, ct_test, cg_train, cg_test,
    Xtrain, Xtest, ybinom_train, ybinom_test):
    # obtain filter
    print('Test data')
    keep_test, grecur  = dp.get_filter(ct=ct_test, cg=cg_test, return_recur=True)
    print('Train data')
    keep_train = dp.get_filter(ct=ct_train, cg=cg_train)
    # apply filter
    X_test,  ybinom_test  = dp.apply_filter(keep_test,  [Xtest,  ybinom_test])
    X_train, ybinom_train = dp.apply_filter(keep_train, [Xtrain, ybinom_train])
    gnames, grecur = dp.apply_filter(keep_test, [cg_test.index.values, np.array(grecur)])
    return (X_test, ybinom_test, X_train, ybinom_train, gnames, grecur)

def main():
    print('DriverPower V0.2')
    start = datetime.now()
    print('Start at {}'.format(start))
    ##
    ## Get parameters
    ##
    print('='*5, 'Obtaining Parameters', '='*5)
    (path_cg_test, path_ct_test, path_cv_test, path_cg_train,
        path_ct_train, path_cv_train, scaler_type, skip_fs,
        output_prefix, output_dir, save_data, preprocessed, scale_first) = get_param()
    ##
    ## Read data
    ##
    print('='*5, 'Reading data', '='*5)
    if preprocessed:
        (X_train, X_test, ybinom_train, ybinom_test,
            gnames, fnames, grecur) = read_npy_data(path_X_train, path_X_test, path_y_train,
            path_y_test, path_gnames, path_fnames, path_grecur)
    else:
        (cg_test, ct_test, X_test, cg_train,
            ct_train, X_train, fnames) = read_raw_data(path_cg_test, path_ct_test,
            path_cv_test, path_cg_train, path_ct_train, path_cv_train)

    if not preprocessed:
        ##
        ## Get response
        ##
        print('='*5, 'Getting response', '='*5)
        print('Test data')
        ybinom_test = dp.get_response(ct_test, cg_test)
        print('Training data')
        ybinom_train = dp.get_response(ct_train, cg_train)
        ##
        ## Scale and Filter
        ##
        if scale_first:
            # scale
            print('='*5, 'Applying {} scaler'.format(scaler_type), '='*5)
            X_train, X_test = scaling(X_train, X_test, scaler_type)
            # filter
            print('='*5, 'Applying filter', '='*5)
            (X_test, ybinom_test, X_train, ybinom_train, gnames, grecur) = filter(ct_train, ct_test, cg_train, cg_test,
                X_train, X_test, ybinom_train, ybinom_test)
        else:
            # filter
            print('='*5, 'Applying filter', '='*5)
            (X_test, ybinom_test, X_train, ybinom_train, gnames, grecur) = filter(ct_train, ct_test, cg_train, cg_test,
                X_train, X_test, ybinom_train, ybinom_test)
            # scale
            print('='*5, 'Applying {} scaler'.format(scaler_type), '='*5)
            X_train, X_test = scaling(X_train, X_test, scaler_type)
        ##
        ## Save data
        ##
        if save_data:
            print('='*5, 'Saving prepocessed data', '='*5)
            save_processed_data(X_train, X_test, ybinom_train, ybinom_test,
                gnames, fnames, grecur, output_dir, output_prefix)
    # More reponse
    ymu_train = ybinom_train[:,0] * 1.0 / ybinom_train.sum(1)
    ylogit_train = sp.special.logit(ymu_train)
    # ymu_test = ybinom_test[:,0] * 1.0 / ybinom_test.sum(1)
    # ylogit_test = sp.special.logit(ymu_test)
    ##
    ## Feature selection
    ##
    if skip_fs:
        print('='*5, 'Skip feature selection', '='*5)
        fs_res = pd.read_table('./output/CRC_exon_robust_feature_selection_res.tsv', header=0)
        assert np.array_equal(fs_res.fname, fnames), 'fnames error'
        cutoff = 0.5
        fset_rndlasso, fidx_rndlasso = feature_score(fs_res.rndlasso, fnames, cutoff=cutoff)
        print('RNDLasso - at cutoff={}, {} selected features are: {}'.format(cutoff, fidx_rndlasso.shape[0], ', '.join(fset_rndlasso)))
    else:
        print('='*5, 'Feature Selection', '='*5)
        lassocv,  fscore_lassocv   = feature_selection(X_train, ylogit_train, 'lassocv')
        rndlasso, fscore_rndlasso  = feature_selection(X_train, ylogit_train, 'rndlasso', alpha=lassocv.alpha_)
        # select feature by cutoff
        cutoff_rndlasso = 0.5
        fset_rndlasso, fidx_rndlasso = feature_score(fscore_rndlasso, fnames, cutoff=cutoff_rndlasso)
        print('RNDLasso - at cutoff={}, {} selected features are: {}'.format(cutoff_rndlasso, fidx_rndlasso.shape[0], ', '.join(fset_rndlasso)))
        cutoff_lassocv  = 0.001
        fset_lassocv, fidx_lassocv = feature_score(fscore_lassocv, fnames, cutoff=cutoff_lassocv)
        print('LassoCV - at cutoff={}, {} selected features are: {}'.format(cutoff_lassocv, fidx_lassocv.shape[0], ', '.join(fset_lassocv)))
        # write fs results to file (tsv)
        fs_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'feature_selection_res.tsv']))
        print('Saving feature selection results to {}'.format(fs_outpath))
        fs_res = pd.DataFrame({'fname': fnames, 'lassocv': fscore_lassocv, 'rndlasso': fscore_rndlasso})
        fs_res.to_csv(fs_outpath, sep='\t', index=False)
    ##
    ## Model
    ##
    print('='*5, 'Predicting mutation rate', '='*5)
    X_train = sm.add_constant(X_train[:, fidx_rndlasso], prepend=False)
    X_test  = sm.add_constant(X_test[:,  fidx_rndlasso], prepend=False)
    glm = sm.GLM(ybinom_train, X_train, family=sm.families.Binomial())
    glm_res = glm.fit()
    mu_pred = glm_res.predict(X_test)
    pval = [binom_test(x, n, p, 'greater') for x, n, p in zip(ybinom_test[:, 0], ybinom_test.sum(1), mu_pred)]
    qval = multipletests(pval, method='fdr_bh')[1]
    glm_res = pd.DataFrame({'binID': gnames, 'Score': -np.log10(pval), 'Pval': pval, 'Qval': qval,
                     'Length': ybinom_test.sum(1).astype(np.int), 'nMut':ybinom_test[:, 0].astype(np.int), 'Recur': grecur.astype(np.int), 'Mu': mu_pred})
    glm_res['Mean'] = (glm_res.Recur + glm_res.nMut)/2
    glm_res['Mean'] = glm_res.Mean.astype(np.int)
    glm_res['PvalM'] = [binom_test(x, n, p, 'greater') for x, n, p in zip(glm_res.Mean, glm_res.Length, glm_res.Mu)]
    glm_res['QvalM'] = multipletests(glm_res.PvalM, method='fdr_bh')[1]
    glm_res = glm_res.sort_values('PvalM')
    glm_res.set_index('binID', inplace=True)
    # write result to file
    mod_outpath = os.path.join(output_dir, '_'.join([output_prefix, 'glm_res.tsv']))
    print('Saving GLM results to {}'.format(mod_outpath))
    glm_res.to_csv(mod_outpath, sep='\t', index=True)
    ##
    ## END
    ##
    end = datetime.now()
    print('End at {}'.format(end))
    print('Use {}'.format(abs(end-start)))

if __name__ == '__main__':
    main()
