#!/usr/bin/env bash
cnv_dir="cnv_Liver-HCC/"
cat ./Liver-HCC.tumor_barcode.txt | while read line
do
	arr=($line)
	barcode=${arr[0]}
	tumor=${arr[1]}
	fname="${barcode}.consensus.20170119.somatic.cna.txt"
	if [[ -f $cnv_dir/$fname ]]; then
		python extract_cnv.py \
		4 74262830 74286838 $tumor $cnv_dir/$fname
	else
		(>&2 echo "$fname not found in $cnv_dir")
	fi
done