''' This script is used to process PanCan.wig.
PanCan.wig:
    $bigWigToWig PanCan.bigwig PanCan.wig
    Values in bigWig are number of pairs (out of 1111) covered according to these criteria: http://www.nature.com/nbt/journal/v31/n3/extref/nbt.2514-S1.pdf, with 14/8 used as cutoffs for tumor/normal (respectively), and 19 as the normal cutoff at a dbSNP site. Up to 40 patients were selected per cohort to calculate this track.

Args:
    cutoff - Percentage of samples have coverage

Usage:
    python process_pancan_wig.py < PanCan.wig
'''
import fileinput

cutoff = 0.4
threshold = 1111 * cutoff

for line in fileinput.input():
    if 'fixedStep' in line: # info row
        print(line.strip())
    else:
        print(1 if int(line.strip()) >= threshold else 0)

