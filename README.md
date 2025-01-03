# GenoBoost

[![GenoBoost](https://github.com/rickyota/genoboost/actions/workflows/genoboost.yml/badge.svg)](https://github.com/rickyota/genoboost/actions/workflows/genoboost.yml)
[![Release](https://github.com/rickyota/genoboost/actions/workflows/publish.yml/badge.svg)](https://github.com/rickyota/genoboost/actions/workflows/publish.yml)
[![Build](https://github.com/rickyota/genoboost/actions/workflows/build.yml/badge.svg)](https://github.com/rickyota/genoboost/actions/workflows/build.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11257869.svg)](https://doi.org/10.5281/zenodo.11257869)


## <a name="started"></a>Getting Started

Download a program for your computer from [here](https://github.com/rickyota/genoboost/releases), and run

```bash
$ genoboost train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --verbose
```

## Table of Contents

- [GenoBoost](#genoboost)
  - [Getting Started](#getting-started)
  - [Table of Contents](#table-of-contents)
  - [News](#news)
  - [Introduction](#introduction)
  - [Users' Guide](#users-guide)
    - [Installation](#installation)
      - [Plink1 Input](#plink1-input)
      - [Plink2 Input (compile)](#plink2-input-compile)
      - [Plink2 Input (docker)](#plink2-input-docker)
    - [Train GenoBoost Model](#train-genoboost-model)
      - [Simplest Usage](#simplest-usage)
      - [Without Validation](#without-validation)
      - [Input Plink2](#input-plink2)
      - [Cross-validation](#cross-validation)
      - [Options for Training](#options-for-training)
    - [Calculate Sample Scores](#calculate-sample-scores)
      - [Simplest Usage](#simplest-usage-1)
      - [Without Validation](#without-validation-1)
      - [Input Plink2](#input-plink2-1)
      - [Cross-validation](#cross-validation-1)
      - [Options for Score](#options-for-score)
  - [Advanced Guide](#advanced-guide)
    - [Advanced Installation](#advanced-installation)
      - [Docker](#docker)
      - [Singularity](#singularity)
    - [Computational Time](#computational-time)

## <a name="news"></a>News

- [v1.2.0](https://github.com/rickyota/genoboost/releases/tag/v1.2.0) (Dec 28, 2024)
    - Clean code and remove some requirements for plink2.
- [v1.1.0](https://github.com/rickyota/genoboost/releases/tag/v1.1.0) (May 23, 2024)
    - Clean code.
- [v1.0.8](https://github.com/rickyota/genoboost/releases/tag/v1.0.8) (Nov 25, 2023)
    - Initial version.
    - Tested on Rocky Linux 8.9 and MacOS 14.3.1.


## <a name="introduction"></a>Introduction

GenoBoost is a polygenic score method to capture additive and non-additive genetic inheritance effects.
So far, most polygenic score methods use the additive model, which exploits an effect size ($\alpha$) and a constant value for each SNV.
GenoBoost exploits three SNV scores ($s_0, s_1, s_2$) for each SNV, corresponding to SNV score for major homozygotes, heterozygotes, and homozygotes. Additive GenoBoost uses linear scores for the SNV scores ($s_2-s_1=s_1-s_0$), and non-additive GenoBoost uses three independent scores to model general non-additive effects.

<img src='readme/img/genoboost-score.png' width=300>

To fit the model, boosting method is used. GenoBoost iteratively selects the most associated SNV with the phenotype and adds the variant to the polygenic score function. When calculating the association, the effects already in polygenic score function are excluded to avoid selecting SNVs with duplicated effects.
After the fitting, the SNV scores for selected SNVs are written in the output folder.

There are two parameters: iteration numbers and learning rates. The default values for learning rates are (0.05, 0.1, 0.2, 0.5).

Covariate effects are fitted only once before starting fitting SNVs using multi-variable logistic regression.

## <a name="user-guide"></a>Users' Guide

For now, the input genotype format is allowed for plink1 or plink2 only.

### <a name="install"></a>Installation

#### <a name="install-plink1"></a>Plink1 Input

If you want to input plink1, download a compiled program for Linux (tested on Rocky Linux<=8.9), macOS (tested on <=14.3.1), and Windows (not tested) from [here][release]. This should take less than 1 minute.

#### <a name="install-plink2-compile"></a>Plink2 Input (compile)

If you want to input plink2 genotype file, you can compile program by yourself as below or [use docker or singularity](#advanced-installation). You can use plink1 format as well.

First, install `rust` as instructed [here][rust-install] if not installed. Then,

```bash
git clone https://github.com/rickyota/genoboost.git
cd genoboost
cargo build --manifest-path ./projects_rust/Cargo.toml --release --bin genoboost
cp ./projects_rust/target/release/genoboost ./genoboost
```

and you can use `genoboost` program. This should take less than 5 minutes.

Using arm architecture, including Macbook M1 and M2 chips, will stop or slow down the software due to the unavailability of SIMD.
I plan to deal with it in the future.

### <a name="train"></a>Train GenoBoost Model

GenoBoost returns the SNV weights file with $s_0, s_1, s_2$ for each SNV in one line.

<img src='readme/img/wgt.png' width=800>

#### <a name="train-simple"></a>Simplest Usage

You can run GenoBoost at least with plink1 genotype files and, in most cases, a covariates file.

See `./example/` for reference of file format. For example, the covariates file should be tab-delimited containing `iid` header corresponding to the sample id in plink1 fam file.

<img src='readme/img/cov.png' width=300>

With the minimum options, GenoBoost produces SNV weights list with the best parameter.
SNV weights list is computed from randomly extracted 80% training samples, and the best parameter is determined in the remaining 20% validation samples. You can control how to split the samples with a random seed.
Write the column name to be used in covariates file after `--cov`.
It is important that the major allele is set to a2 (alternative allele) by `--major-a2-train`since $s_2$ is winsorized. This option is unnecessary if the major allele is already set as the reference allele in genotype file.

```bash
$ ./genoboost train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --seed 55
```

This test code should take less than 2 minutes.

#### <a name="train-train-only"></a>Without Validation

If you want to treat all samples as a training dataset, use `--train-only` option. GenoBoost produces SNV weights each for learning rate. Use `--iter-snv` or `--iter` to control the maximum number of SNVs or iterations for training.

```bash
$ ./genoboost train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --train-only \
    --iter-snv 10000
```

#### <a name="train-plink2"></a>Input Plink2

If you use plink2 genotype file (`.pgen`, `.psam` and `.pvar` or `.pvar.zst`), use `--genot-format plink2` or `--genot-format plink2-vzs`.

If the phenotype is accompanied by covariates in the phenotype file, use `--phe` for the phenotype name. If phenotypes and covariates are in plink2 psam file, do not use `--file-phe`.

Control/case format should be `0/1` or `1/2`.

```bash
$ ./genoboost train \
    --dir ./result \
    --file-genot ./example/genot2 \
    --genot-format plink2-vzs \
    --file-phe ./example/genot2.phe \
    --phe PHENO1 \
    --cov age,sex \
    --major-a2-train \
    --seed 55
```

#### <a name="train-cv"></a>Cross-validation

If you want to run k-fold cross-validation, use `--cross-validation [k]`. GenoBoost will split the samples into k chunks and run k times training with one of the chunks being the validation samples. You can control how to split the samples with a random seed.

```bash
$ ./genoboost train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --cross-validation 5 \
    --seed 55
```

#### <a name="train-option"></a>Options for Training

`--dir <DIR>` : Directory to output.

`--file-genot <FILE>`: Prefix of a plink1 or plink2 file (.bed, .fam, .bim or .pgen, .psam, .pvar/.pvar.zst should exist).

`--genot-format [FORMAT]`: {`plink`, `plink2`, `plink2-vzs`}. Genotype format. Default is `plink`.

`--file-phe [FILE]`: Phenotype or covariates file.

`--phe [NAME]`: Phenotype name specified in --file-phe or psam file.

`--cov [NAMES]`: Covariates names in comma-delimited format. ex. `age,sex,PC1-PC10`.

`--file-sample [FILE]`: Sample file for training. One line for one sample id.

`--file-sample-val [FILE]`: Sample file for validation.

`--file-snv [FILE]`: Snv file for training. One line for one SNV id.

`--major-a2-train`: Set major allele as a2 (alternative allele) in training dataset.

`--iter-snv [NUMBER]`, `--iter [NUMBER]`: Maximum number of SNVs or iterations for training.

`--learning-rates [NUMBERS]`: Learning rates in space-delimited format. Default value is `"0.5 0.2 0.1 0.05"`.

`--cross-validation [NUMBER]`: Number of cross-validations.

`--seed [NUMBER]`: Random seed to control sample split.

`--train-only`: Run without validation.

`--verbose`: Let GenoBoost speak more!

### <a name="score"></a>Calculate Sample Scores

GenoBoost returns a polygenic score for each sample. GenoBoost outputs scores without covariates (`score.tsv`) and with covariates (`score.withcov.tsv`).

<img src='readme/img/score.png' width=200>

#### <a name="score-simple"></a>Simplest Usage

With the minimum options, GenoBoost will calculate sample scores from SNV weights with the best parameters determined in the validation dataset.

```bash
$ ./genoboost score \
    --dir-score ./result_score \
    --dir-wgt ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex
```

#### <a name="score-train-only"></a>Without Validation

If you did not use the validation dataset in the training phase, GenoBoost will output sample scores for all parameters. You have to specify the number of SNVs in `--iters`.

```bash
$ ./genoboost score \
    --dir-score ./result_score \
    --dir-wgt ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --train-only \
    --iters "10 20 50"
```

#### <a name="score-plink2"></a>Input Plink2

Use `--genot-format`, `--file-phe` etc. for plink2 as shown in [training phase](#train-plink2).

```bash
$ ./genoboost score \
    --dir ./result \
    --file-genot ./example/genot2 \
    --genot-format plink2-vzs \
    --file-phe ./example/genot2.phe \
    --cov age,sex
```

#### <a name="score-cv"></a>Cross-validation

If you used cross-validation in the training phase, use `--cross-validation [k]`.

```bash
$ ./genoboost score \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --cross-validation 5
```

#### <a name="score-option"></a>Options for Score

`--dir <DIR>`: Directory to output score files.

`--dir-wgt [DIR]`: The same directory specified for training.

`--file-wgt [FILE]`: Use this specific SNV weight file.

`--file-genot <FILE>`: Prefix of a plink1 or plink2 file (`.bed`, `.fam`, `.bim` or `.pgen`, `.psam`, `.pvar/.pvar.zst` should exist).

`--genot-format [FORMAT]`: {`plink`, `plink2`, `plink2-vzs`}. Genotype format. Default is `plink`.

`--file-phe [FILE]`: Covariates file.

`--cov [NAMES]`: Covariates names in comma-delimited format. ex. `age,sex,PC1-PC10`.

`--file-sample [FILE]`: Sample file for calculating scores. One line for one sample id.

`--iters [NUMBERS]`: Number of SNVs used as a parameter.

`--use-iter`: Also output sample score with the number of iterations as a parameter in addition to the number of SNVs.

`--learning-rates [NUMBERS]`: Learning rates in space-delimited format. Default value is `"0.5 0.2 0.1 0.05"`.

`--cross-validation [NUMBER]`: Number of cross-validations.

`--seed [NUMBER]`: Random seed to control sample split.

`--train-only`: Run without validation.

`--verbose`: Let GenoBoost speak more!

## <a name="advanced-guide"></a>Advanced Guide

### <a name="advanced-installation"></a>Advanced Installation

Using docker or singularity is recommended.

#### <a name="docker"></a>Docker

```bash
$ docker pull rickyota/genoboost:latest \
$ docker run -it -v "$(pwd)":/opt/ rickyota/genoboost:latest \
    train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --verbose
```

#### <a name="singularity"></a>Singularity

```bash
$ singularity build genoboost.sif  ./docker/genoboost.def
$ singularity run --no-home --pwd /opt/ --bind "$(pwd)":/opt/ genoboost.sif \
    train \
    --dir ./result \
    --file-genot ./example/genot \
    --file-phe ./example/genot.cov \
    --cov age,sex \
    --major-a2-train \
    --verbose
```

### <a name="computational-time"></a>Computational Time

For ~216 thousand training samples and ~1.1 million SNVs, GenoBoost would take 10 hours to output weights for 10,000 unique SNVs.

[release]: https://github.com/rickyota/genoboost/releases
[rust-install]: https://www.rust-lang.org/tools/install
