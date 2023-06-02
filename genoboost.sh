#!/bin/bash
#
# GenoBoost training and score
# requires `cargo`

set -eux


# output directory of training
dir_wgt="./test/result/1kg_n10000/train/"
# output directory of score
dir_score="./test/result/1kg_n10000/score/"
# prefix of plink1 file
file_plink="./test/data/1kg_n10000/genot"
# covariate file
file_cov="./test/data/1kg_n10000/genot.cov"
# learning rate parameters
learning_rates="0.1 0.5"

# compile
export RUST_BACKTRACE=full
cargo build --release -p boosting_rust
cp ./target/release/boosting_rust ./genoboost

# train
./genoboost train \
    --dir "$dir_wgt" \
    --file_plink "$file_plink" \
    --file_cov "$file_cov" \
    --learning_rates $learning_rates \
    --iter 100 \
    --clip_sample_weight "top0.1" \
    --prune_snv 0.1

# score
./genoboost score \
    --dir_score "$dir_score" \
    --iters 10 30 50 100 \
    --file_plink "$file_plink" \
    --file_cov "$file_cov" \
    --dir_wgt "$dir_wgt" \
    --learning_rates $learning_rates
