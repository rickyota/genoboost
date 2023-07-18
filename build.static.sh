#!/bin/bash

set -eu


cargo build --manifest-path ./projects_rust/Cargo.toml \
    --release --target=x86_64-unknown-linux-musl \
    --bin genoboost



#function rust_musl_builder(){
#    docker run --rm -it -v "$(pwd)":/home/rust/src messense/rust-musl-cross:x86_64-musl "$@"
#}
#
#rust_musl_builder cargo build --manifest-path ./projects_rust/Cargo.toml \
#    --release --target=x86_64-unknown-linux-musl \
#    --bin genoboost




