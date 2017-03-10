# DriverPower
version: 0.5.0dev

## Introduction

Cancer driver mutations are genomic variants that confer selective growth advantages to tumors, which are rare comparing to the huge amount of passenger mutations. Detecting genomic elements harboring driver mutations is an important yet challenging task, especially for non-coding cancer drivers. DriverPower is a combined burden and functional impact test for coding and non-coding cancer driver elements.

## Installation

DriverPower requires Python >= 3.5 and some other packages. If you don't have Python 3.5 or higher, we highly recommend you to install Python with [Anaconda](https://www.continuum.io/downloads). To install DriverPower, you can either download and install it from [Python Package Index (PyPI)](https://pypi.python.org/pypi/DriverPower/) by typing the following command in your terminal
```bash
$ pip install driverpower
```
or clone this git repository first and install with
```bash
$ git clone https://github.com/smshuai/DriverPower.git
$ cd DriverPower && pip install .
```

To view available sub-commands and help of DriverPower, type
```
$ driverpower -h
usage: driverpower [-h] [-v] {preprocess,select,model} ...

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print the version of DriverPower

The DriverPower sub-commands include:
  {preprocess,select,model}
    preprocess          Load and preprocess data
    select              Run feature selection on preprocessed data
    model               Find driver bins with preprocessed training and test
                        data
```

To view help for a particular sub-commands like `select`, type
```
$ driverpower select -h
usage: driverpower select [-h] -d PATH_DATA [--scaler {robust,standard,none}]
                          [--sampling SAMPLING] [-o OUT]

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -d PATH_DATA, --data PATH_DATA
                        Path to the preprocessed training set (HDF5)

optional parameters:
  --scaler {robust,standard,none}
                        robust or standard (default: robust). Scaler used to
                        scale the data
  --sampling SAMPLING   Number > 0 (default: 1). Sampling the data based on
                        the provided value. Value in (0,1] is used as a
                        fraction. Value > 1 is used as the number of data
                        points.
  -o OUT, --output OUT  Path to the output file (default:
                        ./feature_select.tsv)
```

## A quick example

Suppose we would like to identify potential drivers from protein coding genes and promoter of genes for Liver-HCC.
First, specify test data information in `test_files.tsv`.
```
name element feature  func_names
CDS "./cds.bed.gz"  "./example.Liver-HCC.cds.features.gz" "cadd,dann,eigen_coding"
promoter  "./promCore.bed.gz" "./example.Liver-HCC.promCore.features.gz"  "cadd,dann,eigen_non_coding"
```

```bash
driverpower detect \
  -variant ./example.Liver-HCC.mut.gz \ # observed mutations
  -trainH5 ./trainingData.h5 \ # training data in HDF5 format
  -testFile ./test_files.tsv \ # Name and location of test element sets and features
  -callable ./calllable.bed.gz # whitelist genome regions

```
## Input data requirements

### Variant table

Variant table (-variant) should be a tab-delimited text file **with** header. Variant table records mutations (SNPs, MNPs and indels)  and each row in the table corresponds to one mutation. This table can be derived from MAF or VCF files. A minimum of six columns are required:

1. `chrom`: Chromosome of the mutation in [1-22, X, Y]. Corresponds to `Chromosome` in [MAF Specification](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification).
2. `start`: 0-indexed start coordinate of the mutation. Corresponds to `Start_Position - 1` in MAF.
3. `end`: 1-indexed end coordinate of the mutation. Corresponds to `End_Position` in MAF.
4. `ref`: Reference allele of the mutation. Corresponds to `Reference_Allele` in MAF.
5. `alt`: Mutated allele of the mutation. Corresponds to `Tumor_Seq_Allele2` in MAF.
6. `sid`: Sample identifier.
7. `type`: [OPTIONAL] Type of the mutation, in [SNP, DNP, TNP, ONP, INS, DEL]. Corresponds to `Variant_Type` in MAF. If `type` is not provided, DriverPower will deduce `type` from `ref`and `alt`.
8. One or multiple columns of functional impact scores: [OPTIONAL]. One functional impact measurement per column. If no functional impact score is provided, DriverPower can retrieve scores based on a configuration file.

Example:

| chrom | start   | end     | type | ref | alt   | sid     | CADD  | DANN  | ... |
|-------|---------|---------|------|-----|-------|---------|-------|-------|-----|
| 17    | 7577604 | 7577606 | INS  | -   | AACCT | DO36801 | 6.834 | NA    | ... |
| 17    | 7578405 | 7578406 | SNP  | C   | T     | DO7990  | 4.992 | 0.999 | ... |

### Test files table
Test files table (-testFile) is a CSV file

### Element table

Element table (-element) should be a 4-column BED file (no header). Four columns are in order of chromosome, start, end and binID. Chromosome names should not contain 'chr' and in [1-22, X, Y]. If you have a [BED12](https://genome.ucsc.edu/FAQ/FAQformat#format1) file, an one-liner conversion with [BedTools](http://bedtools.readthedocs.io/en/latest/) could be

```
bedtools bed12tobed6 -i in.bed12 | cut -f1,2,3,4 > out.bed4
```

Example:

(First three coding blocks of *TP53*)

|    |         |         |      |
|----|---------|---------|------|
| 17 | 7565256 | 7565332 | TP53 |
| 17 | 7569523 | 7569562 | TP53 |
| 17 | 7572926 | 7573008 | TP53 |

### Feature table

Feature table is also a tab-delimited text file **with** header. The first column in the table should be `binID`. `binID` in column 1 must be unique. Column 2 to the last column should be features used in mutation rate prediction with one feature per column.

Example:

| binID | GERP    | E128-H3K27ac | ... |
|-------|---------|--------------|-----|
| TP53  | 4.80287 | 1.19475      | ... |
| KRAS  | 3.56563 | 2.53435      | ... |

## Parameters


## LICENSE
DriverPower is distributed under the terms of the [GNU General Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.txt).

## Change Log
- 2017/03/02: Release version 0.4.0, which is used in the PCAWG driver analysis.

## TODO
For v0.5.0:
- New 'detect' module that uses test BEDs and mutations directly, and use multiple scores
- Use configure file to locate functional scores
- DANN scores (and FunSeq2 scores maybe)
- New way of calculating functional impact of elements (average of best of each sample)

Future plans:
- New prediction algorithms (GBM)
