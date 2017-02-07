#!/usr/bin/env bash
sv_dir="sv/"
cat ./Liver-HCC.tumor_barcode.txt | while read line
do
	arr=($line)
	barcode=${arr[0]}
	tumor=${arr[1]}
	fname="${barcode}.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz"
	if [[ -f $sv_dir/$fname ]]; then
		# 500kb
		python extract_sv.py \
		4 73762830 74786838 $tumor $sv_dir/$fname
	else
		(>&2 echo "$fname not found in $sv_dir")
	fi
done