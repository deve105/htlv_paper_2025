#!/bin/bash

set -eou pipefail

fasta="seqkit_norecomb_1201.fasta"

raxml-ng --all \
    --threads auto{12} \
    --seed 12345 \
    --bs-trees autoMRE{1000} \
    --msa-format FASTA \
    --data-type DNA \
    --prefix 2602_model_ \
    --msa "$fasta" \
    --model partition.txt \
    --tree pars{20},rand{20} \
    --bs-cutoff 0.03 \
    --blmin 1e-6 \
    --outgroup L02534.1,KX905203.1,KX905202.1,KF242506.1,KF242505.1,JX891479.1,JX891478.1