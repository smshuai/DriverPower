''' Convert coverage wig to bed. fixedStep wig
'''
import fileinput
import sys

def eprint(*args, **kwargs):
    ''' Print to stderr
    '''
    print(*args, file=sys.stderr, **kwargs)

is_first = False # new bin or not
N = 0 # number of bases
for line in fileinput.input():
    if 'fixedStep' in line: # info line
        # 'fixedStep chrom=1 start=10001 step=1 span=1'
        # eprint(line.strip())
        if 'chrom' in locals():
            print('\t'.join([chrom, str(start), str(end), value]))
        info = [i.split('=') for i in line.strip().split(' ')]
        assert info[1][0] == 'chrom', 'chrom not found'
        chrom = info[1][1]
        assert info[2][0] == 'start', 'start not found'
        start = int(info[2][1]) - 1 # from 1-based to 0-based 
        assert info[3][0] == 'step', 'step not found'
        step = int(info[3][1])
        assert step == 1, 'step is not 1'
        assert info[4][0] == 'span', 'span not found'
        span = int(info[4][1])
        assert span == 1, 'span is not 1'
        is_first = True
    else:
        if N % 50000000 == 0:
            eprint('{}M bases processed'.format(int(N/1000000)))
        N += 1
        if is_first:
            is_first = False
            value = line.strip() # current value
            end = start + 1 # current end coordinate
        else:
            if line.strip() == value:
                # keep start, move end
                end += 1
            else:
                # print one bed record to stdout
                print('\t'.join([chrom, str(start), str(end), value]))
                # move start, end and change value
                start = end
                end += 1
                value = line.strip()
# print last element
print('\t'.join([chrom, str(start), str(end), value]))
eprint('Done! In total, {} bases processed.'.format(N))
