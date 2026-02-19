#!/bin/bash

set -eou pipefail

fasta="2602_htlv_391.fasta"
model=$1

raxml-ng --all \
    --threads auto{16} \
    --seed 12345 \
    --bs-trees autoMRE{5000} \
    --msa-format FASTA \
    --data-type DNA \
    --prefix 2602_${model}_ \
    --msa "$fasta" \
    --model ${model} \
    --tree pars{20},rand{20} \
    --bs-cutoff 0.03 \
    --blmin 1e-6 \
    --brlen scaled \
    --bs-metric fbp,tbe \
    --outgroup L02534.1