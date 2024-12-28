#!/bin/bash


source /bio/lmod/lmod/init/bash
module use /nfs/data06/ricky/app/.modulefiles
module purge
module load clang/16.0.2


cargo clean --manifest-path ./projects_rust/Cargo.toml --release

export RUST_BACKTRACE=full
cargo build --manifest-path ./projects_rust/Cargo.toml --release --bin genoboost 

