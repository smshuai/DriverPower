'''Extract SVs for an element from PCAWG SV calls
Use bedpe.gz file
'''
import sys
import gzip


def is_overlap(a, b):
    ''' a, b are intervals (start, end)
    True if a overlaped b
    '''
    return (min(a[1], b[1]) - max(a[0], b[0])) > 0


def main():
    chrom = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    tumor = sys.argv[4]
    file = sys.argv[5]
    a = (start, end)
    with gzip.open(file, 'rt') as f:
        for line in f:
            row = line.strip().split("\t")
            # compare the first part
            if chrom == row[0]:
                if is_overlap(a, (int(row[1]), int(row[2]))):
                    print(tumor + "\t" + line.strip())
            elif chrom == row[3]:
                # second part of bedpe
                if is_overlap(a, (int(row[4]), int(row[5]))):
                    print(tumor + "\t" + line.strip())


if __name__ == '__main__':
    main()
