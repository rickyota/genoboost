use std::fs::File;
use std::fs::OpenOptions;
use std::io::BufWriter;
use std::path::Path;

pub fn exist_file(fin: &Path) -> bool {
    fin.exists()
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

pub fn is_nonzero(file: &Path) -> bool {
    file_size(file).unwrap() > 0
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

pub fn bufwriter(file: &Path) -> BufWriter<File> {
    let f = match File::create(&file) {
        Ok(f) => f,
        Err(e) => panic!(
            "Cannot create file {:?}, possibly directory does not exist: {}",
            &file, e
        ),
    };
    BufWriter::new(f)
}

pub fn bufwriter_append(file: &Path) -> BufWriter<File> {
    let f = match OpenOptions::new()
        .write(true)
        .append(true)
        .create(true)
        .open(file)
    {
        Ok(f) => f,
        Err(e) => panic!("Error opening file {:?}: {}", file, e),
    };
    BufWriter::new(f)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    #[should_panic]
    fn test_check_exist_file() {
        let fin = PathBuf::from("../../test/data/toy1/does_not_exist.csv");
        check_exist_file(&fin);
    }
}
