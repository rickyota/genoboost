#!/bin/bash
#
# GenoBoost training and score
# requires `cargo`
# 
# occasionally does not stop 

set -eux

# output directory of training
dir_wgt="./result/train/"
# output directory of score
dir_score="./result/score/"
# prefix of plink1 file
file_plink="./example/genot"
# covariate file
file_cov="./example/genot.cov"

# compile
export RUST_BACKTRACE=full
cargo build --manifest-path ./projects_rust/Cargo.toml --release --bin genoboost
cp ./projects_rust/target/release/genoboost ./genoboost

# train
./genoboost train \
    --dir "$dir_wgt" \
    --file-genot "$file_plink" \
    --file-phe "$file_cov" \
    --cov age,sex \
    --major-a2-train \
    --verbose

# score
./genoboost score \
    --dir-score "$dir_score" \
    --dir-wgt "$dir_wgt" \
    --file-genot "$file_plink" \
    --file-phe "$file_cov" \
    --cov age,sex
