#!/usr/bin/env bash
## Please make sure that DriverPower/script/make_features is in PATH
## Prerequisites:
##   bedtools; bigwigaverageoverbed; R with bioconductor
## To install bedtools and bigwigaverageoverbed with conda, run the following line of code
# conda install ucsc-bigwigaverageoverbed bedtools
## To install bioconductor, see https://www.bioconductor.org/install/

if [[ -d ./tmp/ ]]; then rm -r ./tmp/; fi
if [[ -d ./output/ ]]; then rm -r ./output/; fi
if [[ -d ./database/ ]]; then rm -r ./database/; fi

mkdir -p ./output ./tmp/
cp ../../../data/callable.bed.gz ./tmp/

## Download some raw data files
wget --directory-prefix=database/fasta http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr1.fa.gz # FASTA hg19 (chr1 only)
gzip -d database/fasta/chr1.fa.gz
# BED-like format features
wget --directory-prefix=database/bed http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromSynth/wgEncodeOpenChromSynthPanisletsPk.bed.gz
wget --directory-prefix=database/bed http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/E003-DNase.macs2.narrowPeak.gz
# bigwig format features
wget --directory-prefix=database/bw http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E002-H3K4me3.pval.signal.bigwig

## STEP 1: make coverage table and nucleotide content table
bedtools bed12tobed6 -i ./example_element.bed12 > ./tmp/example_element.bed6
# INPUT: ./chr1.fa ./tmp/example_element.bed6 ./tmp/callable.bed.gz
# OUTPUT: ./output/example.cg ./output/example.totcg
generate_cg.sh ./database/fasta/chr1.fa ./tmp/example_element.bed6 ./tmp/callable.bed.gz ./tmp/example.cg ./tmp/example.totcg

## STEP 2: make nucleotide conent features
generate_nuc_covar.py ./tmp/example.cg ./output/example_nuc_features.tsv

## STEP 3: make bigwig features
# INPUT: ./database/bw/ ./example_element.bed12
# OUPUT: ./example_bw_features.tsv
bws2cv.sh ./database/bw/ ./example_element.bed12 ./output/example_bw_features.tsv

## STEP 4: make bed-like features
# INPUT: ./database/bed/ ./example_element.bed6
# OUPUT: ./output/example_bed_features.tsv
beds2cv.sh ./database/bed/ ./tmp/example_element.bed6 ./tmp/example.totcg ./output/example_bed_features.tsv

## STEP 5: combine all features
combine_cv.sh ./output/ ./example_features.tsv

## Check difference
echo "Difference between example_features.tsv and expected_features.tsv"
diff example_features.tsv expected_features.tsv