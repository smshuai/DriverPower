#!/usr/bin/env python
import sys
import os
import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import binom_test, combine_pvalues


def main():
    # input
    tumor = sys.argv[1]
    element = sys.argv[2]
    indir = sys.argv[3]
    outdir = sys.argv[4]
    res_cadd = pd.read_table(os.path.join(indir, '.'.join([tumor, element, 'cadd', 'tsv'])),
        sep='\t', header=0, index_col='binID')
    res_eigen = pd.read_table(os.path.join(indir, '.'.join([tumor, element, 'eigen', 'tsv'])),
        sep='\t', header=0, index_col='binID')
    cutoff_cadd = int(sys.argv[5])
    cutoff_eigen = int(sys.argv[6])
    bed = pd.read_table(sys.argv[7], sep='\t', header=None, usecols=[0,1,2,3])
    bed.columns = ['chr', 'start', 'end', 'binID']
    bed.set_index('binID', inplace=True)
    N = int(sys.argv[8])
    dataset = sys.argv[9]
    res_cadd = tune_cutoff(res_cadd, cutoff_cadd, bed, N)
    res_eigen = tune_cutoff(res_eigen, cutoff_eigen, bed, N)
    res_cadd.columns = ['Pval_CADD', 'Qval_CADD']
    res_eigen.columns = ['Pval_EIGEN', 'Qval_EIGEN']
    res_comb = res_cadd.join(res_eigen)
    res_comb = bed.join(res_comb)
    # use the minimal of CADD and Eigen
    res_comb['Pval_min'] = res_comb.apply(lambda x: min(x.Pval_CADD, x.Pval_EIGEN), axis=1)
    res_comb['Qval_min'] = multipletests(res_comb.Pval_min, method='fdr_bh')[1]
    res_comb.sort_values('Pval_min', inplace=True)
    # fisher
    # res_comb['Pval_fisher'] = res_comb.apply(lambda x: combine_pvalues([x.Pval_CADD, x.Pval_EIGEN],
    #    method='fisher')[1], axis=1)
    # res_comb['Qval_fisher'] = multipletests(res_comb.Pval_fisher, method='fdr_bh')[1]
    # stouffer
    # res_comb['Pval_stouffer'] = res_comb.apply(lambda x: combine_pvalues([x.Pval_CADD, x.Pval_EIGEN],
    #    method='stouffer')[1], axis=1)
    # res_comb['Qval_stouffer'] = multipletests(res_comb.Pval_stouffer, method='fdr_bh')[1]
    res_comb = res_comb.reset_index().loc[:, ('chr', 'start', 'end', 'binID', 'Pval_min', 'Qval_min')]
    res_comb.columns = ('chr', 'start', 'end', 'element_ID', 'p-value', 'q-value')
    res_comb.to_csv(os.path.join(outdir, '.'.join([tumor, element, 'DriverPower', dataset, 'txt'])), sep='\t', index=False)

def tune_cutoff(res, cutoff, bed, N):
    threshold = np.percentile(res.fscore, cutoff)
    res['MuAdj'] = res.Mu * threshold / res.fscore
    res['Pval'] = [binom_test(x, n, p, 'greater') if p<1 else 1 for x, n, p in zip(res.nMutSample, res.Length*N, res.MuAdj)]
    res = bed.join(res)
    res.Pval.fillna(1, inplace=True)
    res['Qval'] = multipletests(res.Pval, method='fdr_bh')[1]
    res.sort_values('Pval', inplace=True)
    return res.loc[:, ('Pval', 'Qval')]

if __name__ == '__main__':
    main()
