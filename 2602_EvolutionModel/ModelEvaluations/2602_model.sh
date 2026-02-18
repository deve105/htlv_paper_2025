#!/bin/bash

set -eouv pipefail

# IQ-TREE 3.0
# Position 1: fasta file
# Position 2: partitions file
# Position 3: output prefix (optional, default: iqtree2_output)
# Based on: https://iqtree.github.io/doc/iqtree-doc.pdf

partitions=$1
name=$2

iqtree3 -s $fasta \
    -T $threads \
    -m MFP+MERGE \
    -p $partitions \
    --prefix $name
    -B 1000
#-mset raxml 
