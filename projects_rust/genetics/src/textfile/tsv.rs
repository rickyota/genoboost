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

    /*
    #[test]
    fn test_compute_num_line() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        //let fin_bim = fin.clone() + ".bim";
        let m = compute_num_line(&fin_bim, None);
        assert_eq!(m, 3);

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        //let fin_fam = fin.clone() + ".fam";
        let n = compute_num_line(&fin_fam, None);
        assert_eq!(n, 5);
    }

    #[test]
    fn test_compute_num_line_compress() {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2Vzs;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        let m = compute_num_line(&fin_bim, Some(String::from("zst")));
        assert_eq!(m, 4);

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        let n = compute_num_line(&fin_fam, None);
        assert_eq!(n, 11);
    }

    #[test]
    fn test_compute_num_column() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        //let fin_bim = fin.clone() + ".bim";
        let col: usize = compute_num_column(&fin_bim, None);
        assert_eq!(col, 6);

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        //let fin_fam = fin.clone() + ".fam";
        let col = compute_num_column(&fin_fam, None);
        assert_eq!(col, 6);
    }

    #[test]
    fn test_compute_num_column_compress() {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2Vzs;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        let m = compute_num_column(&fin_bim, Some(String::from("zst")));
        assert_eq!(m, 5);

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        let n = compute_num_column(&fin_fam, None);
        assert_eq!(n, 4);
    }
    */
}
