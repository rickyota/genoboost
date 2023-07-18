#!/bin/bash

set -eu

# cargo clean

cargo build \
    --release \
    --target=x86_64-unknown-linux-musl \
    --manifest-path ./projects_rust/Cargo.toml \
    --no-default-features \
    --bin genoboost

#function rust_musl_builder(){
#    docker run --rm -it -v "$(pwd)":/home/rust/src messense/rust-musl-cross:x86_64-musl "$@"
#}
#
#rust_musl_builder cargo build \
#    --release \
#    --target=x86_64-unknown-linux-musl \
#    --manifest-path ./projects_rust/Cargo.toml \
#    --no-default-features \
#    --bin genoboost


