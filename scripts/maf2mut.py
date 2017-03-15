#!/usr/bin/env python

import sys

header = sys.stdin.readline().strip().split()
chrom_ix = header.index('Chromosome')
start_ix = header.index('Start_Position')
end_ix = header.index('End_Position')
class_ix = header.index('Variant_Classification')
type_ix = header.index('Variant_Type')
sid_ix = header.index('Tumor_Sample_Barcode')
ref_ix = header.index('Reference_Allele')
alt_ix = header.index('Tumor_Seq_Allele2')

# output
print("\t".join(['chrom', 'start', 'end', 'Variant_Classification', 'type', 'ref', 'alt', 'sid']))

for line in sys.stdin:
    mut = line.strip().split('\t')
    start = int(mut[start_ix]) - 1
    print("\t".join([mut[chrom_ix], str(start), mut[end_ix], mut[class_ix],
                     mut[type_ix], mut[ref_ix], mut[alt_ix], mut[sid_ix]]))
