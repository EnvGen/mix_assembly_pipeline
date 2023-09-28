#!/bin/bash -l

touch annotation/rfam/list_rRNA_genes_in_mix.txt

for filename in annotation/rfam/*_mix_RFAM.tblout; do
    n=$(grep -c 'rRNA' $filename )
    if [ "$n" -gt "0" ]; then
            grep 'rRNA' $filename | awk '{ print $3 }' | sort | uniq >> annotation/rfam/list_rRNA_genes_in_mix.txt #1
    fi
done
