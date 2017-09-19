"""
Make training elements from test elements.

Usage:
    $ python ./make_training_elements.py element.bed bed12 50 exclude.bed


Data requirements:
1. element.bed
    * bed4 or bed12
    * no header
    * e.g. bed4:
        chr1	565900	566050	element_1
        chr1	566760	566910	element_2

Output:
    training_elements.bed

"""

import sys
import pandas as pd
import numpy as np
import seaborn as sns
from pybedtools import BedTool, chromsizes


HG19 = chromsizes('hg19')
CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
          'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
BIN_SIZE = 1000000  # size of window to sampling from


def useBED4(bed4, whitelist, nsample=500):
    '''

    Args:
        bed4:
        whitelist:
        nsample:

    Returns:

    '''
    global HG19, CHROMS, BIN_SIZE
    bed4_df = bed4.to_dataframe()
    # calculate raw sizes
    bed4_df['size'] = bed4_df['end'] - bed4_df['start']
    test_size = bed4_df.pivot_table(index='name', values='size', aggfunc='sum')
    # remove <100bp elements
    test_size = test_size[test_size['size']>=100]
    first = True
    # iterate through chroms
    np.random.seed(97)
    for chrom in CHROMS:
        print('Processing chromosome ' + chrom)
        chrom_size = HG19[chrom][1]  # get chrom size based on hg19
        num_bin = int(np.ceil(chrom_size / BIN_SIZE))
        for bid in range(num_bin):
            # sample nsample start coordinates
            hi = (bid+1)*BIN_SIZE if (bid+1)*BIN_SIZE < chrom_size else chrom_size
            starts = np.random.randint(low=(bid)*BIN_SIZE, high=hi, size=nsample)
            # sample nsample sizes
            sizes = np.random.choice(test_size['size'], nsample, replace=False)
            ends = starts + sizes
            ends[ends>chrom_size] = chrom_size  # capped at chrom.size
            if first:
                train_elements = pd.DataFrame({'chrom': chrom, 'start': starts, 'end': ends})
                first = False
            else:
                tmp = pd.DataFrame({'chrom': chrom, 'start': starts, 'end': ends})
                train_elements = pd.concat([train_elements, tmp])
    # QC elements
    train_elements = train_elements.loc[:, ('chrom', 'start', 'end')]
    train_bed = BedTool.from_dataframe(train_elements)
    # annotate training elements with test elements and callable regions
    train_bed = train_bed.annotate(files=[bed4.fn, whitelist.fn])
    train_elements = train_bed.to_dataframe(names=['chrom', 'start', 'end', 'pct_test', 'pct_callable'])
    # keep <20% test and >80% callable
    keep = np.logical_and(train_elements.pct_test<.2, train_elements.pct_callable>.8)
    train_elements = train_elements[keep]
    # add name and remove duplicates
    train_elements['name'] = train_elements.apply(lambda x: x.chrom + ':' + str(x.start) + '-' + str(x.end), axis=1)
    train_elements.drop_duplicates('name', inplace=True)
    train_elements.sort_values(['chrom', 'start'], inplace=True)
    train_bed = BedTool.from_dataframe(train_elements.loc[:, ('chrom', 'start', 'end', 'name')])
    # summary & plot
    tot_cov = train_bed.total_coverage()
    tot_genome_size = np.sum([HG19[i][1] for i in CHROMS])
    print('Total number of elements made: ' + str(train_elements.shape[0]))
    print('Total percent of genome (hg19) span: ' + str(tot_cov/tot_genome_size*100))
    ax = sns.kdeplot(np.log10(test_size['size']), shade=True, label="input", )
    ax = sns.kdeplot(np.log10(train_elements['end']-train_elements['start']), shade=True, label="output")
    ax.set_xlabel('log10 length of elements (bp)')
    ax.set_ylabel('Density')
    ax.get_figure().savefig('./length_distribution.pdf')
    # save to disk
    train_bed.moveto('./training_elements.bed')


def useBED12(bed12, whitelist, nsample=500):
    '''

    Args:
        bed12:
        callable:

    Returns:

    '''
    # convert bed12 to bed6
    bed6 = bed12.bed6()
    useBED4(bed6, whitelist, nsample)


if __name__ == '__main__':
    assert len(sys.argv) == 5, 'Usage: python ./make_training_elements.py TEST_ELEMENTS FORMAT FOLD WHITELIST'
    test = BedTool(sys.argv[1])
    bedFlag = sys.argv[2]
    fold = int(sys.argv[3])
    whitelist = BedTool(sys.argv[4])
    if bedFlag == 'bed4':
        useBED4(test, whitelist, fold)
    elif bedFlag == 'bed12':
        useBED12(test, whitelist, fold)
    else:
        sys.stderr.write('Unknown input BED format: {}. Must be "bed12" or "bed4"\n'.format(bedFlag))
        sys.exit(1)

