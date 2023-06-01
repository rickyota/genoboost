# GenoBoost v0.2.0

[![GenoBoost](https://github.com/rickyota/genoboost/actions/workflows/genoboost.yml/badge.svg)](https://github.com/rickyota/genoboost/actions/workflows/genoboost.yml)

Polygenic score method for non-additive models.

## Usage

### Download program

Download program for [linux](https://github.com/rickyota/genoboost/releases).

Other versions are available [Release](https://github.com/rickyota/genoboost/releases).


### Docker (Advanced)

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

### Rust (Advanced)

Otherwise, you can directly run GenoBoost with `cargo` and `conda`.

```bash
bash ./genoboost.cv.sh
```
