# DriverPower
version: 0.3.3

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

To see usage and help of DriverPower, type
```bash
$ driverpower -h
```
## Input data requirements

DriverPower requires four input files: mutation table, count table, feature (covariates) table and effective length (coverage) table. All tables should be in TSV format **with** header. Compressed TSV files (*.gzip, *.bz2, *.zip, *.xz) are also acceptable.

### Mutation table

Mutation table should be a tab-delimited text file **with** header. Mutation table records mutations (SNPs, MNPs and indels) in test bins and each row in the table corresponds to one mutation. This table can be derived from MAF or VCF files. Only eight columns are required for DriverPower and extra columns in the table will be ignored:

1. `chrom`: Chromosome of the mutation in [1-22, X, Y]. Corresponds to `Chromosome` in [MAF Specification](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification).
2. `start`: 0-indexed start coordinate of the mutation. Corresponds to `Start_Position - 1` in MAF.
3. `end`: 1-indexed end coordinate of the mutation. Corresponds to `End_Position` in MAF.
4. `type`: Type of the mutation, in [SNP, DNP, TNP, ONP, INS, DEL]. Corresponds to `Variant_Type` in MAF.
5. `ref`: Reference allele of the mutation. Corresponds to `Reference_Allele` in MAF.
6. `alt`: Mutated allele of the mutation. Corresponds to `Tumor_Seq_Allele2` in MAF.
7. `sid`: Sample identifier.
8. `binID`: Bin identifier. 

Example:

| chrom | start   | end     | type | ref | alt   | sid     | binID |
|-------|---------|---------|------|-----|-------|---------|-------|
| 17    | 7577604 | 7577606 | INS  | -   | AACCT | DO36801 | TP53  |
| 17    | 7578405 | 7578406 | SNP  | C   | T     | DO7990  | TP53  |

### Count table

Mutation count table should be a tab-delimited text file (.tsv) **with** header. Mutation count table can be derived from somatic variants files (such as VCF or MAF). The following 4 columns are required:

1. `binID`: The identifier of bins, such as gene or promoter names.
2. `sid`: The identifier of samples.
3. `categ`: The category of mutation count. For SNPs, `categ` should be one of the 96 categories defined in [COSMIC mutational signatures](http://cancer.sanger.ac.uk/cosmic/signatures). For MNPs and indels, `categ` should be 'other'
4. `ct`: The number of mutations in the dataset for this bin, sample and category.

Example:

| binID | sid     | categ | ct |
|-------|---------|-------|----|
| KRAS  | DO49481 | ACC>A | 1  |
| KRAS  | DO51525 | other | 1  |

Note that 'ACC>A' in the first line means ACC to AAC point mutation.

### Feature table

Feature table is also a tab-delimited text file **with** header. The first column in the table should be `binID`. `binID` in column 1 must be unique. Column 2 to the last column should be features used in mutation rate prediction with one feature per column.

Example:

| binID | GERP    | E128-H3K27ac | ... |
|-------|---------|--------------|-----|
| KRAS  | 4.80287 | 1.19475      | ... |
| KRAS  | 3.56563 | 2.53435      | ... |

### Effective length table

Effective length table is also a tab-delimited text file **with** header. The first column in the table should be `binID`. `binID` in column 1 must be unique. Column 2 to column 33 should be 32 possible triple nucleotide contexts: ACA, ACC, ACG, ACT, ATA, ATC, ATG, ATT, CCA, CCC, CCG, CCT, CTA, CTC, CTG, CTT, GCA, GCC, GCG, GCT, GTA, GTC, GTG, GTT, TCA, TCC, TCG, TCT, TTA, TTC, TTG, TTT.

Example:

| binID | ACA | ACC | ... |
|-------|-----|-----|-----|
| KRAS  | 37  | 12  | ... |
| TP53  | 47  | 42  | ... |

## Parameters

- len_threshold=500
- recur_threshold=2
- scaler=['robust', 'standard']
- feature_seletion=['lassocv', 'rndlasso', 'spearman']
- model=['glm', 'xboost', 'dnn']
- func_score=['eigen', 'cadd']
- func_cutoff=80

## LICENSE
DriverPower is distributed under the terms of the [GNU General Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.txt).



## TODO


- MNP and INDEL support for eigen
- 
