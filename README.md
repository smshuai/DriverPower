# DriverPower
version: 0.5

## Input data requirements

DriverPower requires three input files: mutation count table, feature (covariates) table and effective length (coverage) table.

### Mutation counts table

Mutation count table should be a tab-delimited text file (.tsv) **with** header. Mutation count table can be derived from somatic variants files (such as VCF or MAF). The following 4 columns are required (with header), additional columns will be truncated:

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