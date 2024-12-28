fn main() {
    // arg for rustc compile
    // This is same as
    // $ export RUSTFLAGS='-C link-args=-fopenmp'
    // Better use build.rs to make build.sh simple
    //
    println!("cargo:rustc-link-arg=-fopenmp");
    // seems not necessary
    // [ref](https://github.com/rust-lang/cc-rs/issues/266)
    //println!("cargo:rustc-link-lib=gomp");
}
