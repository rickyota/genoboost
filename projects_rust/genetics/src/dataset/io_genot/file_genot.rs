use std::path::PathBuf;

use crate::GenotFormat;

// TODO
#[allow(dead_code)]
#[derive(Clone)]
pub struct FileGenot {
    fgenot: PathBuf,
    format: GenotFormat,
    /*
    // only necessary once when loading first time
    sample_header: bool,
    sample_cols: Vec<String>,
    snv_header: bool,
    snv_cols: Vec<String>,
    */
}

impl FileGenot {
    //pub fn new(fgenot: PathBuf,format:GenotFormat)->Self{
    //}
}
