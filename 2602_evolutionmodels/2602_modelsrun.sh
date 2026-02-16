#!/bin/bash

set -eouv pipefail

fasta="2602_htlv_391.fasta"
threads=$(nproc)

while IFS= read -r line; do
    # Read the fasta and partitions file paths from the line
    partitions=$(echo "$line")
    name="$(basename "$partitions" .txt)"
    echo $name
    echo $partitions
    iqtree3 -s $fasta \
    -T $threads \
    -m MFP+MERGE \
    -p $partitions \
    --prefix $name \
    -B 1000
done < "models.txt"