#!/bin/bash

# DO NOT HAVE TO RUN THIS locally.
# assumed to be called by github workflows

# After `git tag -a v0.3.0 -m "v0.3.0" ` and push to github,
# then github workflows will automatically publish.

set -eu

# ex. genoboost-linux-musl
artifact_name="$1"
target="$2"

d_publish="./${artifact_name}/"


cargo build --manifest-path ./projects_rust/Cargo.toml \
    --release --target=${target} \
    --bin genoboost

mkdir -p ${d_publish}
cp ./projects_rust/target/${target}/release/genoboost ${d_publish}/genoboost

mkdir -p ${d_publish}/sample/
cp ./test/data/toy1/* ${d_publish}/sample/


if [[ ${target} == *"windows"* ]]; then
	tar -cvzf  ./${artifact_name}.zip ${d_publish}
else
	zip -r ./${artifact_name}.zip ${d_publish}
fi

