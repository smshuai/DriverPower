#!/usr/bin/env bash
# run from src/

# Change chrom 23 and 24 in PanCan50p.bed to X and Y
cat ../data/blacklist_region/PanCan50p.bed |
awk 'BEGIN{OFS="\t"} {if ($1==23) $1="X"; else if ($1==24) $1="Y"; print}' > tmp0.bed
awk 'BEGIN{SUM=0} {SUM += $3 - $2} END{ printf "Input %.0f bp\n", SUM }' tmp0.bed

# Use pancan.wig
cat tmp0.bed | awk '{if ($4==1) print }' > tmp1.bed
awk 'BEGIN{SUM=0} {SUM += $3 - $2} END{ printf "Remove PanCan50p, left %.0f bp\n", SUM }' tmp1.bed

# Use hg19_ucsc_gap.bed
cat ../data/blacklist_region/hg19_ucsc_gap.bed | awk '{ sub(/^chr/, "", $0); print }' | bedtools subtract -a tmp1.bed -b stdin > tmp2.bed
awk 'BEGIN{SUM=0} {SUM += $3 - $2} END{ printf "Remove UCSC gaps, left %.0f bp\n", SUM }' tmp2.bed

# Use wgEncodeDacMapabilityConsensusExcludable.bed
cat ../data/blacklist_region/wgEncodeDacMapabilityConsensusExcludable.bed | awk '{ sub(/^chr/, "", $0); print }' | bedtools subtract -a tmp2.bed -b stdin > tmp3.bed
awk 'BEGIN{SUM=0} {SUM += $3 - $2} END{ printf "Remove ENCODE DAC, left %.0f bp\n", SUM }' tmp3.bed

# Use ../data/blacklist_region/wgEncodeDukeMapabilityRegionsExcludable.bed
cat ../data/blacklist_region/wgEncodeDukeMapabilityRegionsExcludable.bed | awk '{ sub(/^chr/, "", $0); print }' | bedtools subtract -a tmp3.bed -b stdin > tmp4.bed
awk 'BEGIN{SUM=0} {SUM += $3 - $2} END{ printf "Remove ENCODE Duke, left %.0f bp\n", SUM }' tmp4.bed

# flatten bed
sort -k1,1 -k2,2n tmp4.bed | bedtools merge -i stdin > tmp5.bed
mv tmp5.bed ../data/blacklist_region/callable50p.bed
rm tmp0.bed tmp1.bed tmp2.bed tmp3.bed tmp4.bed
