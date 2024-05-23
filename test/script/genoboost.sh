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
file_plink="./test/data/1kg_maf0.1_m1k/genot"
# covariate file
file_cov="./test/data/1kg_maf0.1_m1k/genot.cov"

# compile
export RUST_BACKTRACE=full
cargo build --manifest-path ./projects_rust/Cargo.toml --release --bin genoboost
cp ./projects_rust/target/release/genoboost ./genoboost

# train
./genoboost train \
    --dir "$dir_wgt" \
    --file-genot "$file_plink" \
    --file-phe "$file_cov" \
	--learning-rates "0.5 0.2" \
    --cov age,sex \
    --major-a2-train \
    --seed 51


	#--verbose \

# score
./genoboost score \
	--verbose \
    --dir-score "$dir_score" \
    --dir-wgt "$dir_wgt" \
    --file-genot "$file_plink" \
    --file-phe "$file_cov" \
    --cov age,sex
