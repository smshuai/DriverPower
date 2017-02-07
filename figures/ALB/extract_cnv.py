'''Extract CNV for an element from PCAWG CNV calls
'''
import sys


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
    with open(file) as f:
        for line in f:
            row = line.strip().split("\t")
            if chrom == row[0]:
                if is_overlap(a, (int(row[1]), int(row[2]))):
                    # find the CNV for element
                    print(tumor + "\t" + line.strip())


if __name__ == '__main__':
    main()
