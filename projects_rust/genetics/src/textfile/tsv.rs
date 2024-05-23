//! use tsv instaed of text module
//! use "csv" crate
//!
//!
//!TODO: you can directly load to class [ref](https://docs.rs/csv/latest/csv/tutorial/index.html#reading-with-serde)
//!
//!

/* use crate::textfile;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path; */

/* // test
//fn load_first_content<T: std::io::Read>(buf: &mut BufReader<T>) -> String {
fn load_first_content<T: std::io::Read>(buf: T) -> String {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(buf);
    rdr.records().next().unwrap().unwrap()[0].to_owned()
} */

#[cfg(test)]
mod tests {
    //use super::*;
    /*     use crate::{io_genot, GenotFormat};
    use std::path::PathBuf; */

    /*     #[test]
    fn test_load_first_content() {
        let content = "abc\tdef\n1\trs5\n";
        //let mut buf=BufReader::new(content.as_bytes());
        let content = content.as_bytes();
        let content = load_first_content(content);
        assert_eq!(content, "abc");
    } */
}
