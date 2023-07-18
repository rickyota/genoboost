#!/bin/bash

set -eu

cargo build \
    --release \
    --target=x86_64-unknown-linux-musl \
    --manifest-path ./projects_rust/Cargo.toml \
    --no-default-features \
    --bin genoboost
