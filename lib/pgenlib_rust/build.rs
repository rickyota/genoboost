// https://rendered-obsolete.github.io/2018/09/30/rust-ffi-ci.html

// later: https://rust-lang.github.io/rust-bindgen/non-system-libraries.html
// if you have .o : https://rust-lang.github.io/rust-bindgen/tutorial-3.html

// unnecessary?
//extern crate bindgen;
//use bindgen;

use std::env;
use std::path::PathBuf;

use bindgen::CargoCallbacks;
use cmake;

fn main() {

    // cmake
    let dst = cmake::build("src/pgenlib");
    println!("cargo:rustc-link-search=native={}", dst.display());
    println!("cargo:rustc-link-lib=dylib=stdc++");
    // ok?
    // should be same as project name in CMakeLists.txt
    println!("cargo:rustc-link-lib=static=pgenlib");
    
    println!("cargo:rustc-link-arg=-fopenmp");
    println!("cargo:rustc-link-lib=gomp");



    // binding
    let bindings = bindgen::Builder::default()
        .header("src/pgenlib/pgenlibr_wrapc.h")
        .parse_callbacks(Box::new(CargoCallbacks))
        .clang_arg("-xc++")
        .clang_arg("-std=c++11") 
        .allowlist_file("src/pgenlib/pgenlibr_wrapc.h")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings");
}

