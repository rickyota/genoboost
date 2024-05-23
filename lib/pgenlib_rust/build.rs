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
    
    // [ref](https://stackoverflow.com/questions/50642574/how-can-i-specify-linker-flags-arguments-in-a-build-script)
    println!("cargo:rustc-link-arg=-fopenmp");

    println!("cargo:rustc-link-lib=gomp");
    // added 230823
    // [ref](https://github.com/rust-or/highs-sys/blob/master/build.rs)
    //println!("cargo:rustc-link-lib=dylib=gomp");



    // binding
    let bindings = bindgen::Builder::default()
        .header("src/pgenlib/pgenlibr_wrapc.hpp")
        .parse_callbacks(Box::new(CargoCallbacks))
        .clang_arg("-xc++")
        .clang_arg("-std=c++11") 
        .clang_arg("-stdlib=libc++") 
        .allowlist_file("src/pgenlib/pgenlibr_wrapc.hpp")
        .generate()
        .expect("Unable to generate bindings");
    //let bindings = bindgen::Builder::default()
    //    .header("src/pgenlib/pgenlibr_wrapc.h")
    //    .parse_callbacks(Box::new(CargoCallbacks))
    //    .clang_arg("-xc++")
    //    .clang_arg("-std=c++11") 
    //    .allowlist_file("src/pgenlib/pgenlibr_wrapc.h")
    //    .generate()
    //    .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings");
}

