#!/usr/bin/env python
'''This script is used to generate mutation table and count table for DriverPower.

'''
import sys
import pandas as pd
import numpy as np
import tempfile
from subprocess import Popen, PIPE


def main():
    # mutation could be in MAF or MAF-like format **without** header
    mut_path = sys.argv[1]
    # BED no header. Columns should be chrom, start, end, binID...
    bed_path = sys.argv[2]
    # callable
    callable_path = sys.argv[3]
    # 0-index of chrom, pos (start, 1-indexed), ref, alt, sample_id
    # ex. of Broad sim.: '0,1,4,5,3'
    col_idx = [int(i) for i in sys.argv[4].split(',')]
    # two output files with header
    out_mut = sys.argv[5]
    out_ct = sys.argv[6]
    # process mut
    mut = pd.read_table(mut_path, sep='\t', header=None, usecols=col_idx,
                        dtype={col_idx[0]: str, col_idx[1]: np.int})
    print("Input {} mutations".format(mut.shape[0]))
    # Extract based on col_idx
    mut = mut.loc[:,col_idx]
    mut.columns = ['chrom', 'pos', 'ref', 'alt', 'sid']
    # Add type column
    mut['type'] = 'Unknown'
    len_ref = mut.ref.apply(len)
    len_alt = mut.alt.apply(len)
    # SNP len_ref == len_alt == 1 and no "-"
    # DNP len_ref == len_alt == 2 and no "-"
    # TNP len_ref == len_alt == 3 and no "-"
    # ONP len_ref == len_alt >= 4 and no "-"
    eq_len = len_ref == len_alt
    mut.loc[np.logical_and(eq_len, len_ref == 1), 'type'] = 'SNP'
    mut.loc[np.logical_and(eq_len, len_ref == 2), 'type'] = 'DNP'
    mut.loc[np.logical_and(eq_len, len_ref == 3), 'type'] = 'TNP'
    mut.loc[np.logical_and(eq_len, len_ref >= 4), 'type'] = 'ONP'
    # DEL and INS
    mut.loc[mut.alt == '-', 'type'] = 'DEL'
    mut.loc[mut.ref == '-', 'type'] = 'INS'
    # Add start column
    mut['start'] = mut['pos'] - 1
    # Add end column
    # SNP: end = pos
    # DNP: end = pos+1
    # TNP: end = pos+2
    # ONP: end = pos+len(ref)-1
    # INS: end = pos+1
    # DEL: end = pos+len(ref)-1
    mut['end'] = mut.pos
    mut.loc[mut.type == 'DNP', 'end'] = mut.loc[mut.type == 'DNP', 'pos'] + 1
    mut.loc[mut.type == 'TNP', 'end'] = mut.loc[mut.type == 'TNP', 'pos'] + 2
    mut.loc[mut.type == 'ONP', 'end'] = \
            mut.loc[mut.type == 'ONP', 'pos'] + len_ref[mut.type == 'ONP'] - 1
    mut.loc[mut.type == 'DEL', 'end'] = \
            mut.loc[mut.type == 'DEL', 'pos'] + len_ref[mut.type == 'DEL'] - 1
    mut.loc[mut.type == 'INS', 'end'] = mut.loc[mut.type == 'INS', 'pos'] + 1
    # convert mutation to bed
    mut = mut.loc[:, ('chrom', 'start', 'end', 'type', 'ref', 'alt', 'sid')]
    # save mut
    tmp_mut = tempfile.NamedTemporaryFile('w')
    mut.to_csv(tmp_mut, sep='\t', header=False, index=False) # could be slow for large file
    # intersect with bed
    bed = pd.read_table(bed_path, sep='\t', header=None)
    bed = bed.loc[:, (0,1,2,3)]
    bed.columns = ('chrom', 'start', 'end', 'binID')
    # 'chr1' to '1'
    if 'chr' in bed.chrom[1]:
        bed.chrom = bed.chrom.apply(lambda x: x[3:])
    # save bed to tmp file
    tmp_bed = tempfile.NamedTemporaryFile('w')
    bed.to_csv(tmp_bed, sep='\t', header=False, index=False)
    # run bedtools with subprocess
    intersect1 = Popen(['bedtools', 'intersect', '-a', tmp_mut.name, '-b', tmp_bed.name,
                        '-wa', '-wb'], stdout=PIPE)
    intersect2 = Popen(['bedtools', 'intersect', '-a', 'stdin', '-b', callable_path],
                       stdin=intersect1.stdout, stdout=PIPE)
    tmp_mut.close()
    tmp_bed.close()
    # pivot table
    mut = pd.read_table(intersect2.stdout, sep='\t', header=None)
    # Add a header
    print("Output {} mutations that intersect bins and are callable".format(mut.shape[0]))
    mut.columns = ('chrom', 'start', 'end', 'type', 'ref', 'alt', 'sid',
                  'chrom_bin', 'start_bin', 'end_bin', 'binID')
    ct = mut.pivot_table(index=['binID', 'sid'], values=['chrom'], aggfunc=len)
    ct.columns = ['ct']
    # save outputs
    mut.to_csv(out_mut, sep='\t', index=False)
    ct.to_csv(out_ct, sep='\t')
    print('Done!')

if __name__ == '__main__':
    main()
