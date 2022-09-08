#!/bin/bash
#
# GenoBoost score

set -eux

# output directory
dir_score="./test/result/1kg_n10000/score/"
# directory of training
dir_wgt="./test/result/1kg_n10000/train/"
# prefix of plink1 file
file_plink="./test/data/1kg_n10000/genot"
# covariate file
file_cov="./test/data/1kg_n10000/genot.cov"

ns_iter="1 2 3 5 10 50"

export RUST_BACKTRACE=full

cargo run --release --bin boosting_rust -- \
    score \
    --dir_score "$dir_score" \
    --file_plink "$file_plink" \
    --iters $ns_iter \
    --file_cov "$file_cov" \
    --dir_wgt "$dir_wgt"
