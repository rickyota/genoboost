use walkdir::WalkDir;

fn main() {
    let dirs = ["./lib/pgenlib/"];
    // exclude src/unuse
    //let dirs = [
    //    "./lib/pgenlib/src/",
    //    "./lib/pgenlib/src/include/",
    //    "./lib/pgenlib/src/simde/",
    //    "./lib/pgenlib/src/simde/x86",
    //    "./lib/pgenlib/src/simde/x86/avx512",
    //];

    // all files with .cc or .cpp
    let cpps: Vec<String> = dirs
        .map(|dir| {
            WalkDir::new(dir)
                .into_iter()
                .map(|x| x.unwrap().path().display().to_string())
        })
        .into_iter()
        .flatten()
        .filter(|x| x.ends_with(".cpp") || x.ends_with(".cc"))
        .collect();

    // cannot print in build.rs
    //println!("cpps: {:?}", cpps);

    // how to add openmp?
    // https://users.rust-lang.org/t/binding-openmp-c-function/40196/4

    // should be .rs with #[cxx::bridge]
    cxx_build::bridge("src/lib.rs")
        .files(&cpps)
        .flag_if_supported("-fopenmp")
        .compile("pgenlib-bridge");

    //.flag_if_supported("-static")

    // seems not necessary
    //.flag_if_supported("-std=c++")
    //.flag_if_supported("-std=c++11")

    // seems not necessary
    //.flag_if_supported("-lgomp")

    // arg for rustc compile
    // This is same as
    // $ export RUSTFLAGS='-C link-args=-fopenmp'
    // Better use build.rs to make build.sh simple
    //
    println!("cargo:rustc-link-arg=-fopenmp");
    // seems not necessary
    // [ref](https://github.com/rust-lang/cc-rs/issues/266)
    //println!("cargo:rustc-link-lib=gomp");

    println!("cargo:rerun-if-changed=src/lib.rs");
    // all cpp and h
    dirs.map(|dir| {
        WalkDir::new(dir)
            .into_iter()
            .map(|x| x.unwrap().path().display().to_string())
    })
    .into_iter()
    .flatten()
    .filter(|x| {
        x.ends_with(".cpp") || x.ends_with(".cc") || x.ends_with(".hpp") || x.ends_with(".h")
    })
    .for_each(|x| println!("cargo:rerun-if-changed={:?}", x));
}
