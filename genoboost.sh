#!/bin/bash
#
# GenoBoost training

set -eux

# output directory
dir_wgt="./test/result/1kg_n10000/train/"
# prefix of plink1 file
file_plink="./test/data/1kg_n10000/genot"
# covariate file
file_cov="./test/data/1kg_n10000/genot.cov"

export RUST_BACKTRACE=full
cargo run --release --bin boosting_rust -- \
    train \
    --dir "$dir_wgt" \
    --file_plink "$file_plink" \
    --file_cov "$file_cov" \
    --iter 100 \
    --learning_rate 0.1 \
    --clip_sample_weight "top0.1" \
    --prune_snv 0.1
