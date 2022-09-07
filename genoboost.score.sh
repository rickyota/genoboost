#!/bin/bash
#
# GenoBoost score

set -eux

file_plink="./test/data/toy1/genot"
file_cov="./test/data/toy1/genot.cov"
dir_wgt="./test/result/toy1/train/"
dir_score="./test/result/toy1/score/"

ns_iter="1 2 3 5 10 100"

export RUST_BACKTRACE=full

cargo run --release --bin boosting_rust -- \
    score \
    --dir_score "$dir_score" \
    --file_plink "$file_plink" \
    --iters $ns_iter \
    --file_cov "$file_cov" \
    --dir_wgt "$dir_wgt"
