# GenoBoost v0.3.0

[![GenoBoost](https://github.com/rickyota/genoboost/actions/workflows/genoboost.yml/badge.svg)](https://github.com/rickyota/genoboost/actions/workflows/genoboost.yml)
[![Release](https://github.com/rickyota/genoboost/actions/workflows/publish.yml/badge.svg)](https://github.com/rickyota/genoboost/actions/workflows/publish.yml)

GenoBoost is a polygenic score method to capture non-additive genetic inheritance effects.
So far, most polygenic score methods use the additive model.
GenoBoost exploits three SNV scores for each SNV to model general non-additive effects.

## Usage

Download a folder including program and toy example [here](https://github.com/rickyota/genoboost/releases).

See toy example for file format.

### Train GenoBoost model

```bash
$ genoboost train \
    --dir ./result \
    --file_plink ./example/genot \
    --file_cov ./example/genot.cov \
    --file_sample ./example/genot.train.sample
```

--dir: Directory to output.

--file_plink: Prefix of a plink1 file (.bed, .bim, .fam should exist).

--file_cov: Covariate file.

--file_sample: Sample file for training.

--file_phe: [optional] Phenotype file. If not set, phenotype in --file_plink is used.

--phe: [optional] Phenotype name indicated in --file_phe.

--boost_type: Genetic inheritance model to use. "nonadd" or "add". Default is "nonadd".


### Calculate sample scores

```bash
$ genoboost score \
    --dir_score ./result_score \
    --dir_wgt ./result \
    --file_plink ./example/genot \
    --file_cov ./example/genot.cov \
    --file_sample ./example/genot.train.sample
```

--dir_score: Directory to output.

--dir_wgt: Same directory indicated on training. 

--file_plink: Prefix of a plink1 file (.bed, .bim, .fam should exist).

--file_cov: Covariate file.

--file_sample: Sample file for training.



## Advanced Usage

### Docker

Using docker or singularity is recommended.

Run GenoBoost on a example dataset in `./test/data/1kg_n10000` (1000 samples x 10000 SNVs).

```bash
$ docker run -td \
    -v "$(pwd)/test/data/1kg_n10000":/work/data:ro -v "$(pwd)/result":/work/result \
    rickyota/genoboost:latest \
    bash ./genoboost.docker.cv.sh
```

or

```bash
$ singularity build geno.sif docker://rickyota/genoboost:latest
$ singularity exec \
    --bind "$(pwd)/test/data/1kg_n10000":/work/data,"$(pwd)/result":/work/result \
    --no-home --pwd /opt/genoboost geno.sif \
    bash ./genoboost.docker.cv.sh
```

Result files are now in `./result/` .

### Rust

Otherwise, you can directly run GenoBoost with `cargo` and `conda`.

```bash
bash ./genoboost.cv.sh
```
