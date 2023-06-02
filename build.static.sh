#!/bin/bash

set -eu

cargo build --manifest-path ./projects_rust/Cargo.toml \
    --release --target=x86_64-unknown-linux-musl \
    --bin genoboost


