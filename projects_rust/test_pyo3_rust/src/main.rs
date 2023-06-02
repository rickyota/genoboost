//! test for apple M1 mac
//! output error on avx (able to run without Rosetta)
//! `cargo run --bin test_rust --target x86_64-apple-darwin`
//! 
//! run `conda activate geentics_pyo3"

use test_pyo3_rust;

fn main() {
    println!("Hello, world!");

    println!("test");
    test_pyo3_rust::test();

    println!("Done!");
}
