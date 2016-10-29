''' Functional adjustment module for DriverPower
'''

import os
import tabix
import logging
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import binom_test

# create  logger
logger = logging.getLogger('FUNC ADJ')
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


def get_eigen(mut, path='/u/sshuai/sshuai/func_score/eigen/v1.1', coding=True):
    ''' Get eigen scores with tabix
    '''
    # make chrom string
    mut['chrom'] = mut.chrom.astype(str)
    # eigen only support chr 1-22
    valid_chrom = [str(i) for i in range(1,23)]
    mut = mut[mut.chrom.isin(valid_chrom)]
    # create return table
    # SNP only. TO DO: ADD MNP and INDEL Support
    keep = mut['type'] == 'SNP'
    # chrom != X, Y
    eigen = mut[keep].copy()
    if eigen.shape[0] == 0:
        logger.warning('No mutations left in Eigen adjustment')
        return None
    if coding:
        logger.info('Retrieving Eigen Coding Scores')
        # name for version 1.1. A single file for coding.
        name = 'Eigen_hg19_coding_annot_04092016.tab.bgz'
        file_path = os.path.join(path, name)
        assert os.path.isfile(file_path), 'Cannot find eigen coding in {}'.format(file_path)
        # open one eigen
        tb = tabix.open(file_path)
        # row apply
        func = lambda x: query_eigen_SNP(tb, x[0], x[1], x[2], x[4], x[5])
        eigen['fscore'] = eigen.apply(func, axis=1)
    else:
        logger.info('Retrieving Eigen Non-Coding Scores')
        # One file per chrom for non-coding. 22 in total.
        file_dict = {str(i) : 'Eigen_hg19_noncoding_annot_chr{}.tab.bgz'.format(i) for i in range(1, 23)}
        # check files, must be 22 True
        file_names = np.array(list(file_dict.values()))
        check_file = np.array([os.path.isfile(f) for f in file_names])
        assert np.sum(check_file) == 22, 'Cannot find eigen noncoding in {}'.format(", ".join(file_names[~check_file]))
        # open 22 eigen files
        file_dict = {k: tabix.open(v) for k, v in file_dict.items()}
        # row apply
        func = lambda x: query_eigen_SNP(file_dict[str(x[0])], x[0], x[1], x[2], x[4], x[5])
        eigen['fscore'] = eigen.apply(func, axis=1)
    return eigen


def query_eigen_SNP(tb, chrom, start, end, ref, alt):
    ''' Find eigen score for a SNP
    '''
    # logger.debug((chrom, start, end, ref, alt))
    res = tb.query(str(chrom), start, end)
    for i in res: # iter through records
        # check ref
        assert ref == i[2], 'Reference allele in mutation table does not match eigen records'
        if alt == i[3]:
            return np.float(i[len(i)-3]) # Eigen-phred, last but three col
    # No score found. (silent for coding)
    return np.nan


def query_eigen_MNP():
    ''' Find eigen score for a MNP
    '''
    pass


def query_eigen_indel():
    ''' Find eigen score for a INDEL
    '''
    pass

def get_cadd():
    ''' Get CADD scores with tabix
    '''
    pass


def func_adj(res, mut, method, path_eigen, is_coding, cutoff=85):
    ''' Main wrapper for functional adjustment
    '''
    support_method = ['eigen']
    assert method in support_method, 'Invalid functional score method. Must be chosen from {}'.format(support_method)
    if method == 'eigen':
        eigen = get_eigen(mut, path_eigen, is_coding)
        # mean eigen score per bin
        ## For coding, mean score of non-silent mutations
        fscore = eigen.dropna().pivot_table(index='binID', values='fscore', aggfunc=np.mean)
    res['fscore'] = fscore
    # fill bin without fscore with 0
    res.fscore.fillna(0, inplace=True)
    threshold = np.percentile(res.fscore, cutoff)
    # adjusted mutation rate
    ## for bin w/o fscore, MuAdj will be inf.
    res['MuAdj'] = res.Mu * threshold / res.fscore
    # Pval: binom(nMutSample, length, MuAdj)
    res['Pval'] = [binom_test(x, n, p, 'greater') if p<1 else 1 for x, n, p in zip(res.nMutSample, res.Length, res.MuAdj)]
    res['Qval'] = multipletests(res.Pval, method='fdr_bh')[1]
    return res
