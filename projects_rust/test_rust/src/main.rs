//! test for apple M1 mac
//! output error on avx (able to run without Rosetta)
//! `cargo run --bin test_rust --target x86_64-apple-darwin`

use test_rust;

fn main() {
    println!("Hello, world!");

    println!("test");
    test_rust::test();

    println!("Done!");
}
