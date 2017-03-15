#!/usr/bin/env python
''' Convert VCFs to mutation table in DriverPower
Usage:
    vcfs2mut.py VCF_LIST > MUT
Args:
    VCF_LIST - a two-columns tsv file. The first column is donor ID, and the second column is .vcf path
Note:
    VCFs are one per donor.
'''

import sys

def process_vcf(path, sid):
    valid_dna = set('ATCG,')
    Ncomma = 0
    with open(path) as vcf:
        for mut in vcf:
            if mut[0] == '#': continue
            chrom, pos, mutid, ref, alt = mut.split('\t')[:5]
            pos = int(pos)
            # check ref and alt
            if not (set(ref).issubset(valid_dna) and set(alt).issubset(valid_dna)):
                sys.stderr.write('WARNING: other letters rather than ATCG in ref or alt for variant {}'.format(mut))
                continue
            if ',' in alt:
                # alt is like G,T
                Ncomma += 1
                alt = alt.split(',')[0]
            # remove 'chr'
            chrom = chrom.replace('chr', '')
            if len(ref) == len(alt) == 1:
                # SNP
                start = pos - 1
                end = pos
                mut_type = 'SNP'
            elif len(ref) == 1 and len(alt) > 1:
                # INS
                # from: chr1 16079825 C CA
                # to: chr1 16079824 16079826 - A
                ref = '-'
                alt = alt[1:]
                start = pos - 1
                end = pos + 1
                mut_type = 'INS'
            elif len(alt) == 1 and len(ref) > 1:
                # DEL
                # from: chr1 2258532 GAA G
                # to: 1 2258532 2258534 AA - 
                alt = '-'
                ref = ref[1:]
                start = pos
                end = pos+len(ref)
                mut_type = 'DEL'
            else:
                sys.stderr.write('ERROR: unrecognized variant {}'.format(mut))
                sys.exit(1)
            dp_mut = [chrom, str(start), str(end), mut_type, ref, alt, sid]
            print("\t".join(dp_mut))
    return Ncomma

if __name__ == '__main__':
    with open(sys.argv[1]) as vcfs:
        print("\t".join(['chrom', 'start', 'end', 'type', 'ref', 'alt', 'sid']))
        Ncomma = 0
        for row in vcfs:
            if row[0] == '#': continue
            sid, path = row.strip().split('\t')
            Ncomma += process_vcf(path, sid)
        sys.stderr.write('Num comma in alt = {}'.format(Ncomma))

