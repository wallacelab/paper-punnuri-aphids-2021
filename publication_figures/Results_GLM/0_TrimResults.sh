#! /bin/bash

# This code snippet was used to trim the original TASSEL GLM results files down to only the necessary columns to save space

for infile in *.csv; do
    cat  $infile  | cut -f1-4,6-7 -d',' > tmp.csv
    mv -f tmp.csv $infile
    gzip $infile
    break
done
