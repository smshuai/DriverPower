#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import binom_test


def main():
    # input
    path_cadd = sys.argv[1]
    path_eigen = sys.argv[2]
    path_output = sys.argv[3]
    res_cadd = pd.read_table(path_cadd, sep='\t', header=0, index_col='binID')
    res_eigen = pd.read_table(path_eigen, sep='\t', header=0, index_col='binID')
    assert np.array_equal(res_cadd.shape, res_eigen.shape)
    cutoff_cadd = float(sys.argv[4])
    cutoff_eigen = float(sys.argv[5])
    N = int(sys.argv[6])
    length = res_cadd.Length.astype(int)
    nMut = res_cadd.nMut.astype(int)
    nSample = res_cadd.nSample.astype(int)
    if cutoff_cadd > 0 and cutoff_eigen <= 0:
        # use cadd only
        print('CADD')
        res = tune_cutoff(res_cadd, cutoff_cadd, N)
    elif cutoff_cadd > 0 and cutoff_eigen > 0:
        # use both
        print('CADD')
        res_cadd = tune_cutoff(res_cadd, cutoff_cadd, N)
        print('EIGEN')
        res_eigen = tune_cutoff(res_eigen, cutoff_eigen, N)
        res_cadd.columns = ['Pval_CADD', 'Qval_CADD']
        res_eigen.columns = ['Pval_EIGEN', 'Qval_EIGEN']
        res = res_cadd.join(res_eigen)
        res['Pval_min'] = res.apply(lambda x: min(x.Pval_CADD, x.Pval_EIGEN), axis=1)
        not_na = np.logical_not(res.Pval_min.isnull())
        res.loc[not_na, 'Qval_min'] = multipletests(res.Pval_min[not_na], method='fdr_bh')[1]
        res = res.loc[:, ('Pval_min', 'Qval_min')]
        res.sort_values('Pval_min', inplace=True)
    else:
        # use eigen only
        print('EIGEN')
        res = tune_cutoff(res_eigen, cutoff_eigen, N)
    res.columns = ['p-value', 'q-value']
    res['length'] = length
    res['nMut'] = nMut
    res['nSample'] = nSample
    res.to_csv(path_output, sep='\t', na_rep='NA')


def tune_cutoff(res, cutoff, N):
    # use trimmed fscore (nMut>0)
    res.index.name = 'element_ID'
    threshold = res.fscore[res.nMut>0].fillna(0).quantile(cutoff/100)
    print('threshold:', threshold)
    res['MuAdj'] = res.Mu * threshold / res.fscore
    res['Pval'] = [binom_test(x, n, p, 'greater') if p<1 else 1 for x, n, p in zip(res.nMutSample, res.Length*N, res.MuAdj)]
    # set NA
    na_filter = np.logical_or(res.nSample < 1, res.Length <= 100)
    # Pval 1
    res.loc[res.nSample == 1, 'Pval'] = 1
    res.loc[na_filter, 'Pval'] = np.nan
    res.loc[np.logical_not(na_filter), 'Qval'] = multipletests(res.Pval[np.logical_not(na_filter)], method='fdr_bh')[1]
    res.sort_values('Pval', inplace=True)
    res.loc[na_filter, 'Qval'] = np.nan
    return res.loc[:, ('Pval', 'Qval')]

if __name__ == '__main__':
    main()
