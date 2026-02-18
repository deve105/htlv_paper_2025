#!/bin/bash

set -eou pipefail

fasta="2602_htlv_391.fasta"

while IFS= read -r line || [[ -n "$line" ]]; do
    # Read the fasta and partitions file paths from the line
    partitions=$(echo "$line")
    name="$(basename "$partitions" .txt)"
    echo $name
    echo $partitions
    iqtree3 -s $fasta \
    -T AUTO \
    -m MFP \
    -p $partitions \
    --prefix $name \
    -B 1000
done < $1