#!/bin/bash
#
# GenoBoost training

set -eux

file_plink="./test/data/toy1/genot"
file_cov="./test/data/toy1/genot.cov"
dir_wgt="./test/result/toy1/train/"

n_iter=20

# delete wgt
fwgt="${dir_wgt}wgt/boosting.wgt"
echo "delete wgt file: ${fwgt}"
rm -f "$fwgt"

export RUST_BACKTRACE=full

cargo run --release --bin boosting_rust -- \
    train \
    --dir "$dir_wgt" \
    --file_plink "$file_plink" \
    --iter $n_iter \
    --file_cov "$file_cov"
