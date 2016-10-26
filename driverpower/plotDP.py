import numpy as np
import scipy as sp

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns
from statsmodels.graphics.api import abline_plot
sns.set(style="white")

def plot_corr(pred, obs):
    ''' Plot predicted value vs. observed value
    '''
    rval = sp.stats.pearsonr(pred, obs)[0]
    fig, ax = plt.subplots()
    ax.scatter(pred, obs)
    abline_plot(ax=ax, slope=1, intercept=0)
    ax.set_title('R = {}'.format(rval))
    ax.set_ylabel('Observed values')
    ax.set_xlabel('Predictd values')
    ax.set_xlim([min(min(pred), min(obs)), max(max(pred), max(obs))])
    ax.set_ylim([min(min(pred), min(obs)), max(max(pred), max(obs))])
    return fig, ax

def bar_plot(heights, labels):
    ''' Plot a barplot with vertical labels
    '''
    fig, ax = plt.subplots()
    n = len(heights)
    ax.bar(np.arange(n), heights)
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_xticklabels(labels, rotation='vertical')
    return fig, ax

def pval_qqplot(pvals, labels=None):
    ''' Plot a pval qqplot and label top points when provided.
    Args:
        pvals - List or np.array of pvals
        labels - List or np.array of top point labels. Default None.
    '''
    obs_pval = -np.log10(np.sort(pvals))
    exp_pval = -np.log10( np.arange(1, len(obs_pval)+1)/float(len(obs_pval+1)))
    fig, ax = plt.subplots()
    ax.scatter(exp_pval, obs_pval)
    abline_plot(ax=ax, slope=1, intercept=0)
    ax.set_ylabel('Observed p-values')
    ax.set_xlabel('Expected p-values')
    if labels:
        ntop = len(labels)
        for g, x, y in zip(labels, exp_pval[:ntop], obs_pval[:ntop]):
            ax.annotate(g, xy=(x+0.1, y+0.1))
    return fig, ax