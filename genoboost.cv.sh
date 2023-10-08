#!/bin/bash
#
# GenoBoost 5-fold cross-validation
# requires `conda` and `cargo`

set -eux

# output directory
dir="./result/cross_validation/"
# prefix of plink1 file
file_plink="./test/data/1kg_maf0.1_m1k/genot"
# covariate file
file_cov="./test/data/1kg_maf0.1_m1k/genot.cov"

# compile
export RUST_BACKTRACE=full
cargo build --manifest-path ./projects_rust/Cargo.toml --release --bin genoboost
cp ./projects_rust/target/release/genoboost ./genoboost

# train
./genoboost train \
    --dir "${dir}/train" \
    --file-genot "$file_plink" \
    --file-phe "$file_cov" \
    --cov age,sex \
    --cross-validation 5 \
    --major_a2_train

# score
./genoboost score \
    --dir-score "${dir}/score" \
    --dir-wgt "${dir}/train" \
    --file-genot "$file_plink" \
    --file-phe "$file_cov" \
    --cov age,sex \
    --cross-validation 5
