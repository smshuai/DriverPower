''' Functional adjustment module for DriverPower
'''

import os
import sys
import tabix
import logging
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import binom_test
from driverpower.model import do_binom_test
# create  logger
logger = logging.getLogger('FUNC ADJ')


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
        file_dict = {str(i) : os.path.join(path, 'Eigen_hg19_noncoding_annot_chr{}.tab.bgz'.format(i)) for i in range(1, 23)}
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

def get_cadd(mut, path='/u/sshuai/sshuai/func_score/cadd/v1.3'):
    ''' Get CADD scores with tabix
    '''
    # make chrom string
    mut['chrom'] = mut.chrom.astype(str)
    # create return table
    # SNP only. TO DO: ADD MNP and INDEL Support
    keep = mut['type'].isin(['SNP','DEL','INS'])
    cadd = mut[keep].copy()
    if cadd.shape[0] == 0:
        logger.warning('No mutations left in CADD adjustment')
        return None
    logger.info('Retrieving CADD SNP Scores')
    # name for version 1.3.
    # A single file for SNP.
    # Use pre-computed PCAWG indel scores
    snp   = 'whole_genome_SNVs.tsv.gz'
    indel = 'PCAWG.INDELS.CADD.v1.3.tsv.gz'
    snp_path = os.path.join(path, snp)
    indel_path = os.path.join(path, indel)
    assert os.path.isfile(snp_path), 'Cannot find CADD SNP scores in {}'.format(snp_path)
    assert os.path.isfile(indel_path), 'Cannot find CADD PCAWG indel scores in {}'.format(indel_path)
    # open one CADD
    tb_snp = tabix.open(snp_path)
    tb_indel = tabix.open(indel_path)
    # row apply
    func = lambda x: query_cadd(x[3], tb_snp, tb_indel, x[0], x[1], x[2], x[4], x[5])
    cadd['fscore'] = cadd.apply(func, axis=1)
    return cadd


def query_cadd(categ, tb_snp, tb_indel, chrom, start, end, ref, alt, phred=True):
    ''' Find CADD score for a mutation based on type
    '''
    if categ == 'SNP':
        return query_cadd_SNP(tb_snp, chrom, start, end, ref, alt, phred)
    elif categ in ['INS', 'DEL']:
        return query_cadd_indel(tb_indel, chrom, start, end, ref, alt, phred)
    else:
        logger.warning('No CADD score because mutation {}:{}-{}({}>{}) is not an indel or a SNP'.format(chrom,start,end,ref,alt))
        return np.nan


def query_cadd_SNP(tb, chrom, start, end, ref, alt, phred=True):
    ''' Find CADD score for a SNP
    '''
    if phred:
        idx = 5 # use CADD phred score
    else:
        idx = 4 # use CADD raw score
    # logger.debug((chrom, start, end, ref, alt))
    res = tb.query(str(chrom), start, end)
    for i in res: # iter through records
        # check ref
        assert ref == i[2], 'Reference allele in mutation table does not match CADD records'
        if alt == i[3]:
            return np.float(i[idx]) # CADD raw or phred
    # No score found
    return np.nan


def query_cadd_indel(tb, chrom, start, end, ref, alt, phred=True):
    ''' Find CADD score for an indel
    '''
    idx = 6 if phred else 5
    res = tb.query(str(chrom), start, end)
    for i in res:
        if ref == i[3] and alt == i[4]:
            return np.float(i[idx])
    # No score found
    return np.nan

def load_func_scores(func_conf_path, methods):
    ''' Load the configure file for functional scores
    Args:
        func_conf_path - str, path to the conf file
        methods - list, name of methods will be used
    '''
    conf = pd.read_csv(func_conf_path, header=0)
    support_methods = set([i.upper() for i in conf.name.unique()]) # all in upper case
    methods = set([i.upper() for i in methods])
    assert support_methods.issuperset(methods), \
            'Cannot find configuration for the following method(s): {}'.format(', '.join(methods-support_methods))


def func_adj(res, mut, method, dir_func, is_coding, cutoff=85):
    ''' Main wrapper for functional adjustment
    '''
    support_method = ['eigen', 'cadd']
    assert method in support_method, 'Invalid functional score method. Must be chosen from {}'.format(support_method)
    if method == 'eigen':
        dir_eigen = os.path.join(os.path.expanduser(dir_func), 'eigen', 'v1.1')
        eigen = get_eigen(mut, dir_eigen, is_coding)
        # mean eigen score per bin
        ## For coding, mean score of non-silent and silent POINT mutations
        fscore = eigen.fillna(0).pivot_table(index='binID', values='fscore', aggfunc=np.mean)
        res['nEigenAll'] = eigen.pivot_table(index='binID', values='fscore', aggfunc=len)
        res['nEigen'] = eigen.dropna().pivot_table(index='binID', values='fscore', aggfunc=len)
    elif method == 'cadd':
        dir_cadd = os.path.join(os.path.expanduser(dir_func), 'cadd', 'v1.3')
        cadd = get_cadd(mut, dir_cadd)
        fscore = cadd.fillna(0).pivot_table(index='binID', values='fscore', aggfunc=np.mean)
        res['nCaddAll'] = cadd.pivot_table(index='binID', values='fscore', aggfunc=len)
        res['nCadd'] = cadd.dropna().pivot_table(index='binID', values='fscore', aggfunc=len)
    res['fscore'] = fscore
    # fill bin without fscore with 0
    res.fscore.fillna(0, inplace=True)
    threshold = np.percentile(res.fscore, cutoff)
    logger.info('Use functional score {} ({} percentile) to adjust mutation rate'.format(threshold, cutoff))
    # adjusted mutation rate
    ## for bin w/o fscore, MuAdj will be inf.
    res['MuAdj'] = res.Mu * threshold / res.fscore
    # Pval: binom(nMutSample, length, MuAdj)
    res['Pval'] = [binom_test(x, n, p, 'greater') if p<1 else 1 for x, n, p in zip(res.nMutSample, res.Length, res.MuAdj)]
    res['Qval'] = multipletests(res.Pval, method='fdr_bh')[1]
    res.sort_values('Pval', inplace=True)
    return res
##
# v0.5.0 Detect
##
def func_adj_new(result, ftuple, N, use_gmean):
    ''' Perform functional adjustment.
    Args:
        result - pd.DF, having BGMR, func_scores, nSample, nMut, length
        ftuple - tuple, (func_name, func_thresh, func_cut)
        N - int, number of donors
    Return:
        padj - np.array, adjusted pvals
        qadj - np.array, adjusted qvals
    '''
    # check func_name
    if ftuple[0] not in result.columns:
        logger.error('Bin-level functional score name {} not found in result'.format(ftuple[0]))
        sys.exit(1)
    func_cut = ftuple[2] if ftuple[1] is None else result[ftuple[0]][result.nMut>0].fillna(0).quantile(ftuple[1]/100)
    logger.info('Cutoff - {} - {}'.format(ftuple[0], func_cut))
    muadj = result['BGMR'] * func_cut / result[ftuple[0]]
    padj, qadj = do_binom_test(result, N, muadj, use_gmean)
    return padj, qadj
