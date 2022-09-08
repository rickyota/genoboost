#!/bin/bash
#
# GenoBoost training

set -eux

file_plink="./test/data/1kg_n10000/genot"
file_cov="./test/data/1kg_n10000/genot.cov"
dir_wgt="./test/result/1kg_n10000/train/"

# delete wgt
fwgt="${dir_wgt}wgt/boosting.wgt"
echo "delete wgt file: ${fwgt}"
rm -f "$fwgt"

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


