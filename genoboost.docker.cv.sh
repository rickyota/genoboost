#!/bin/bash
#
# GenoBoost 5-fold cross-validation

set -eux

# mount to this path
dir_data="/work/data/"
dir_result="/work/result/"
# output directory
dir="${dir_result}"
# prefix of plink1 file
file_plink="${dir_data}genot"
# covariate file
file_cov="${dir_data}genot.cov"
# learning rate parameters
learning_rates="0.1 0.5"

# output directory of samples
dir_sample="${dir}samples/"
# output directory of cross-validation
dir_cv="${dir}cross_validation/"

# create cv dataset
mkdir -p "$dir_sample"
eval "$(conda shell.bash hook)"
conda activate genoboost
python -m projects.genetics_py.src.dataset \
    --cross_validation \
    --cross_validation_n 5 \
    --dout "${dir_sample}" \
    --fplink "$file_plink"

# train
for cvi in {0..4}; do
    dir_wgt_cv="${dir_cv}tr.cv${cvi}/"
    fin_sample="${dir_sample}tr.cv${cvi}.samples"
    ./genoboost train \
        --dir "$dir_wgt_cv" \
        --file_plink "$file_plink" \
        --file_cov "$file_cov" \
        --file_sample "$fin_sample" \
        --learning_rates $learning_rates \
        --iter 100 \
        --clip_sample_weight "top0.1" \
        --prune_snv 0.1
done

# score
for cvi in {0..4}; do
    dir_wgt_cv="${dir_cv}tr.cv${cvi}/"
    dir_score_cv="${dir_cv}va.cv${cvi}/"
    fin_sample="${dir_sample}va.cv${cvi}.samples"
    ./genoboost score \
        --dir_score "$dir_score_cv" \
        --iters 10 30 50 100 \
        --file_plink "$file_plink" \
        --file_cov "$file_cov" \
        --file_sample "$fin_sample" \
        --dir_wgt "$dir_wgt_cv" \
        --learning_rates $learning_rates

    dir_score_cv="${dir_cv}ts.cv${cvi}/"
    fin_sample="${dir_sample}test.samples"
    ./genoboost score \
        --dir_score "$dir_score_cv" \
        --iters 10 30 50 100 \
        --file_plink "$file_plink" \
        --file_cov "$file_cov" \
        --file_sample "$fin_sample" \
        --dir_wgt "$dir_wgt_cv" \
        --learning_rates $learning_rates
done
