#!/bin/bash
#
# GenoBoost training and score

set -eux

# mount to this path
dir_data="/work/data/"
dir_result="/work/result/"
# output directory of training
dir_wgt="${dir_result}train/"
# output directory of score
dir_score="${dir_result}score/"
# prefix of plink1 file
file_plink="${dir_data}genot"
# covariate file
file_cov="${dir_data}genot.cov"
# learning rate parameters
learning_rates="0.1 0.5"

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
