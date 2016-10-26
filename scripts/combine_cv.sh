#!/usr/bin/env bash

cv_dir=$1
out_path=$2

cvs="$cv_dir/*.tsv"

first=true

for cv in $cvs
do
    if [[ ! -f $cv ]]
    then
        echo "ERROR: $cv not found. Aborting"
        exit 1
    else
        echo "$cv"
    fi
    
    # sort cv by binID
    (head -n 1 $cv && tail -n +2 $cv | sort -k1,1) > cv.sorted
    cv="cv.sorted"

    if [ "$first" == true ]
    then
        cp $cv $out_path
        first=false
        continue
    fi

    # compare binIDs (col1)
    difference=`diff <(cut -f1 $out_path) <(cut -f1 $cv)`
    if [ "$difference" == '' ] # identical
    then
        paste $out_path <(cut -f2- $cv) > tmp.combine
        mv tmp.combine $out_path
    else
        echo "ERROR: binID in $cv is not the same as $out_path. Aborting"
        exit 1
    fi
done

rm cv.sorted
exit 0
