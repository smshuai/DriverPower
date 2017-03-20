''' Helper functions
'''

def get_filter(ct, cg, len_threshold=500, recur_threshold=2, return_recur=False, return_tab=False):
    ''' Obtain filter based on length and recurrence.
    Args:
        ct - pd.DF. Count table, 4 cols, no index
        cg - pd.DF. Coverage table, indexed by binID
        len_threshold   - int. Bins with effective length < len_threshold will be discarded
        recur_threshold - int. Bins mutated in < recur_threshold samples will be discarded
        return_recur    - bool. Whether or not return the recur
        return_tab      - bool. Whether or not return the filter table.
    Return:
        keep  - np.array. Index of bins passed all filters
        recur - pd.Series. Indexed by binID. Recurrence of all bins.
        filter_tab - pd.DF. Indexed by binID. Two columns, cg and recur
    '''
    # ct pivot by binID and sid
    ct_byb_bys = ct.pivot_table(values='ct', index='binID', columns='sid', aggfunc=sum).fillna(0)
    # Sum of cg per binID
    cg_sum = cg.sum(axis=1).copy()
    nbin = cg_sum.shape[0]
    logger.info('{} bins before filter'.format(nbin))
    # Recurrence table
    recur = (ct_byb_bys > 0).sum(axis=1)
    filter_tab = pd.concat([cg_sum, recur], axis=1)
    filter_tab.columns = ['cg', 'recur']
    filter_tab.fillna(0, inplace=True) # recur can be 0
    assert np.array_equal(filter_tab.index, cg.index), "Filter index does not match CG index"
    # 1. at least mutated in recur_threshold samples
    keep1 = np.where(filter_tab.recur >= recur_threshold)[0] # index of bins pass
    logger.info('{} ({:.2f}%) bins have mutations in at least {} samples'.format(keep1.shape[0], keep1.shape[0]/nbin*100, recur_threshold))
    # 2. at least len_threshold bp
    keep2 = np.where(filter_tab.cg >= len_threshold)[0] # index of bins pass
    logger.info('{} ({:.2f}%) bins have effective length >= {} bp'.format(keep2.shape[0], keep2.shape[0]/nbin*100, len_threshold))
    # intersection of keep1 and keep2
    keep = np.intersect1d(keep1, keep2)
    logger.info('{} ({:.2f}%) bins passed all filters'.format(keep.shape[0], keep.shape[0]/nbin*100))
    if return_recur:
        return keep, filter_tab.recur
    if return_tab:
        return keep, filter_tab
    return keep # index

def get_gmean(y, recur):
    ''' Use binomial response y and recur to produce a new gmean response.
    '''
    logger.info('Use geometric mean as response')
    ynew = np.zeros(y.shape, dtype=int)
    ynew[:,0] = np.sqrt(recur * y[:,0])
    ynew[:,1] = y.sum(1) - ynew[:,0]
    return ynew

