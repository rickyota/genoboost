use std::fs::File;
use std::path::Path;

pub fn exist_file(fin: &Path) -> bool {
    fin.exists()
    //Path::new(fin).exists()
}

pub fn check_exist_file(fin: &Path) {
    if !exist_file(fin) {
        panic!("File does not exist: {:?}", fin);
    }
}

/// This contains exist_file()
pub fn able_open_file(fin: &Path) -> bool {
    File::open(fin).is_ok()
}

pub fn file_size(file: &Path) -> Option<usize> {
    File::open(file)
        .ok()
        .map(|f| f.metadata().unwrap().len() as usize)
}

pub fn check_open_file(fin: &Path) {
    // This is included in open error, but written to clarify error message
    if !exist_file(fin) {
        panic!("File does not exist: {:?}", fin);
    }
    if let Err(_) = File::open(fin) {
        panic!("Cannot open file: {:?}", fin);
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    //use crate::{io_genot, GenotFormat};

    use super::*;

    #[test]
    #[should_panic]
    fn test_check_exist_file() {
        let fin = PathBuf::from("../../test/data/toy1/does_not_exist.csv");
        check_exist_file(&fin);
    }

}