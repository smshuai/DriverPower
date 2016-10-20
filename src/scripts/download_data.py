import synapseclient
from subprocess import run
'''
Resources page:
    https://wiki.oicr.on.ca/display/PANCANCER/PCAWG-2%2C5%2C9%2C14+final+analysis+set+resources
Python ENV:
    source activate py35
'''

syn = synapseclient.Synapse()
syn.login('shimin.shuai@mail.utoronto.ca','shuai930123')

# cds bed
cds_bed_id = 'syn7180178'
cds_bed = syn.get(cds_bed_id)
run(['mv', cds_bed.path, '../annotation/'])

# MAF October_2016_whitelist_2583.snv_mnv_indel.maf.gz
maf_id = 'syn7364923'
maf = syn.get(maf_id)
run(['mv', maf.path, '../data/'])

# MuTect_coverage
cv_id = 'syn7289561'
cv = syn.get(cv_id)
run(['mv', cv.path, '../data/'])


