#!/usr/bin/env bash
cnv_dir="../consensus.20170119.somatic.cna/"
cat ../ALB/Liver-HCC.tumor_barcode.txt | while read line
do
	arr=($line)
	barcode=${arr[0]}
	tumor=${arr[1]}
	fname="${barcode}.consensus.20170119.somatic.cna.txt"
	if [[ -f $cnv_dir/$fname ]]; then
		python ../ALB/extract_cnv.py \
		9 104182842 104198062 $tumor $cnv_dir/$fname
	else
		(>&2 echo "$fname not found in $cnv_dir")
	fi
done