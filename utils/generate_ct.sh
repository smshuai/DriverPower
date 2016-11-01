#!/usr/bin/env bash
<<generate_ct.sh
This script is used to generate mutation count over bed for DriverPower.
    Usage:
        generate_ct.sh MAF BED CALLABLE OUTMUT OUTCT
    Args:
        MAF      - Input MAF, use PCAWG final release
        BED      - Input BED, use 4/5 columns
        CALLABLE - Input callable regions, use data/blacklist_region/callable50p.bed
        OUTMUT   - Output mutation table path 
        OUTCT    - Output count table path
    Output count table have the following columns:
        binID - ID of bins, could be non-unique. 4th column in BED.
        sid   - ID of samples (donors). 43rd column in MAF.     
        categ - For SNP, categ is the triple nucleotide context plus the to-base, such as 'ACT>A'.
                For indel/DNP/TNP/ONP, categ is 'other'  
        ct    - Mutation count for this bin, sample and categ. Only mutations inside CALLABLE will be considered.
generate_ct.sh

MAF=$1
BED=$2
CALLABLE=$3
OUTMUT=$4
OUTCT=$5

# generate random tmp dir
tmp_dir=./tmp.`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32`
if [[ ! -d $tmp_dir ]]
then
    # dir not exist
    mkdir $tmp_dir
else
    (>&2 echo "Warning: $tmp_dir exist. Overwriting.")
    rm -r $tmp_dir
    mkdir $tmp_dir
fi

# input number of mutations
nmut_in=($(wc -l $MAF))
(>&2 echo "INFO: Input number of mutations: ${nmut_in[0]}")

(>&2 echo "STEP: Converting MAF to BED")
cut -f2,3,4,7,8,10,16,43 $MAF | # cut chrom, start, end, variant_type, ref, alter, ref_context, Donor_ID
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3, $4, $5, $6, toupper(substr($7,10,3)), $8}' > $tmp_dir/maf.bed # convert to BED (0-based)
# ref_context is 21bp long. substr($7,10,3) gives 3 nucleotide context.

(>&2 echo "STEP: Filtering with callable regions")
bedtools intersect -a $tmp_dir/maf.bed -b $CALLABLE -wa > $tmp_dir/mut.callable.bed # create callable mutation table
rm $tmp_dir/maf.bed

# nmut after callable filter
nmut_callable=($(wc -l $tmp_dir/mut.callable.bed))
(>&2 echo "INFO: Callable number of mutations: ${nmut_callable[0]}")

# Create categ.
(>&2 echo "STEP: Creating category for mutations")
cat $tmp_dir/mut.callable.bed | 
awk 'BEGIN {OFS="\t"; a["T"]="A";a["A"]="T";a["C"]="G";a["G"]="C"}
{
    if ($4!="SNP") # not a snp
        print $0, "other";
    else if ($4=="SNP" && $5=="A") # from=A, categ is reverse complement triple > to
        # a[substr($7,3,1)] complement of the 3rd base in $7
        print $0, a[substr($7,3,1)]a[substr($7,2,1)]a[substr($7,1,1)]">"a[$6];
    else if ($4=="SNP" && $5=="C") # from=C,  categ is triple nucleotide > to
        print $0, $7">"$6;
    else if ($4=="SNP" && $5=="G") # from=G, categ is reverse complement triple > to
        print $0, a[substr($7,3,1)]a[substr($7,2,1)]a[substr($7,1,1)]">"a[$6];
    else if ($4=="SNP" && $5=="T") # from=T, same as C
        print $0, $7">"$6;
    else
    {
        print "ERROR: " $0 "\nFrom base for SNP is not A, C, G, or T" > "/dev/stderr";
        exit 1
    } 

}' > $tmp_dir/mut.callable.with.categ.bed
rm $tmp_dir/mut.callable.bed

# Add binID
(>&2 echo "STEP: Intersecting BED")
cat $BED | awk '{ sub(/^chr/, "", $0); print }' | # remove 'chr' in chrom, if any.
bedtools intersect -a $tmp_dir/mut.callable.with.categ.bed -b stdin -wa -wb > $tmp_dir/mut.in.bin.bed

nmut_bin=($(wc -l $tmp_dir/mut.in.bin.bed))
(>&2 echo "INFO: Number of mutations intersected BED: ${nmut_bin[0]}")

# Pivot table with python pandas
(>&2 echo "STEP: Creating count table")
pivot_ct.py $tmp_dir/mut.in.bin.bed $tmp_dir/ct.tsv


# add header to mut
head -1 ~/current/data/PDAC/exon/PDAC.exon.mut.tsv > $OUTMUT
cat $tmp_dir/mut.in.bin.bed >> $OUTMUT
mv $tmp_dir/ct.tsv $OUTCT

rm -r $tmp_dir

(>&2 echo "STEP: Finished")
exit 0


