#!/bin/bash

# assumed to be called by github workflows

set -eu

d_publish="./${1}/"

bash ./build.static.sh
mkdir -p ${d_publish}
cp ./projects_rust/target/x86_64-unknown-linux-musl/release/genoboost ${d_publish}/genoboost

mkdir -p ${d_publish}/sample/
cp ./test/data/toy1/* ${d_publish}/sample/


zip -r ./genoboost.zip ${d_publish}

