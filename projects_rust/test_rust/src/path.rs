use std::ffi::OsStr;
use std::ffi::OsString;
use std::fs::{self, File};
use std::path::{Path, PathBuf};

pub fn f<T: AsRef<Path>>(fin: T) {
    match File::open(fin) {
        Ok(_) => true,
        Err(_) => false,
    };
}

pub fn g_osstr(f: &Path, ext: impl AsRef<OsStr>) -> PathBuf {
    let mut fnew: OsString = f.into();
    fnew.push(".");
    fnew.push(ext.as_ref());
    fnew.into()
}

// this is fine
pub fn g_path(f: &Path, ext: impl AsRef<Path>) -> PathBuf {
    let mut fnew: OsString = f.into();
    fnew.push(".");
    fnew.push(ext.as_ref());
    fnew.into()
}

pub fn test() {
    /*
    let mut str = String::from("./test_dir/");
    //let mut str = String::from("./test_dir");
    // work for both
    str = str + "/";

    let mut p = PathBuf::from(str);
    //let mut p = PathBuf::from("./test_dir/");
    println!("path {:?}", p);
    p = p.canonicalize().unwrap();
    println!("canonical path {:?}", p);
    fs::create_dir_all(&p).unwrap();

    p.push("test.wgt");
    println!("wgt {:?}", p);

    f(p);
     */

    let p = PathBuf::from("./test_dir/abc");
    let p_ext = g_osstr(&p, "csv");
    println!("p_osstr {:?}", p_ext);

    let p = PathBuf::from("./test_dir/abc");
    let p_ext = g_osstr(&p, "csv");
    println!("p_path {:?}", p_ext);

    let p = PathBuf::from("./test_dir/abc");
    let p2 = PathBuf::from("csv");
    let p_ext = g_osstr(&p, &p2);
    println!("p2_path {:?}", p_ext);

    // if PathBuf->String, Path->str, then
    // &PathBuf => &str  : yes
    let p = PathBuf::from("./test_dir/abc");
    let p2: &Path = &p;
    println!("p2 {:?}", p2);
    // to bring back
    let p = PathBuf::from("./test_dir/abc");
    let p2: &Path = &p;
}
