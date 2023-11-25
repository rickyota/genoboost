#!/bin/bash
#
# GenoBoost training and score

set -eux

# output directory of training
dir_wgt="./result/train/"
# output directory of score
dir_score="./result/score/"
# prefix of plink1 file
file_plink="./example/genot"
# covariate file
file_cov="./example/genot.cov"

function genoboost-docker() {
    docker run -it -v "$(pwd)":/opt/ --env RUST_BACKTRACE=full rickyota/genoboost:latest "$@"
}

# train
genoboost-docker train \
    --dir "$dir_wgt" \
    --file-genot "$file_plink" \
    --file-phe "$file_cov" \
    --cov age,sex \
    --major-a2-train \
    --seed 51

# score
genoboost-docker score \
    --dir-score "$dir_score" \
    --dir-wgt "$dir_wgt" \
    --file-genot "$file_plink" \
    --file-phe "$file_cov" \
    --cov age,sex
