#!/usr/bin/env bash

driverpower detect \
    --variant    "./data/example.mut.gz" \
    --trainH5    "./data/example.train.h5" \
    --testFile   "./test_list.tsv" \
    --callable   "./annot/callable.bed.gz" \
    --cohortName "example.CNS-Medullo" \
    --outDir     "./output/"
