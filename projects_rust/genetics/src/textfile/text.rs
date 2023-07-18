//! if input is fin: &Path, then return Option ??
//!
//! When you use <R: std::io::Read>, you cannot load twice in the function,
//! so explicitly use &[u8]
//!

use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::time::Instant;

pub fn isin_header<R: std::io::Read>(col: &str, buf: R) -> bool {
    let header = load_first_line(buf);
    isin_string(col, &header)
}

/// col: "ID"
/// line: "#ID\tABC"
pub fn isin_string(col: &str, line: &str) -> bool {
    line.contains(&col)
}

fn load_first_line<R: std::io::Read>(buf: R) -> String {
    let buf = BufReader::new(buf);
    // unwrap() twice since lines() returns Option<Result<>>
    buf.lines().next().unwrap().unwrap()

    /*
    let mut buf = BufReader::new(buf);
    let mut line = String::new();
    //let mut line = String::with_capacity(512);
    //let mut buf = BufReader::new(File::open(fin).unwrap());
    buf.read_line(&mut line).unwrap();
    // remove last \n
    line = line.trim_end().to_string();
    println!("{}", line);
    line
    */

    /*
    let mut buf_a = BufReader::new(buf);
    for line in buf_a.lines() {
        println!("{}", line.unwrap());
    }
    println!("loop end2");
    panic!();

    */
}

/* // TODO: <P: AsRef<Path>>
//pub fn read_file_to_end(fin: &Path, compress: Option<&str>) -> Option<Vec<u8>> {
pub fn read_file_to_end(fin: &Path, compress: Option<&str>) -> Vec<u8> {
    let mut buf: Vec<u8> = vec![];

    match compress {
        None => {
            // -> you cannot load stdin twice [ref](https://stackoverflow.com/questions/17958571/reading-input-from-stdcin-twice)
            // -> You have to load only once at first and use them
            let _ = File::open(fin)
                .unwrap_or_else(|_| panic!("Cannot open file {:?}", fin))
                .read_to_end(&mut buf)
                .unwrap_or_else(|_| panic!("Cannot read file {:?}", fin));
        }
        Some(compress) => match compress {
            "zst" => {
                buf = zstd::stream::decode_all(
                    File::open(fin).unwrap_or_else(|_| panic!("Cannot open file {:?}", fin)),
                )
                .unwrap_or_else(|_| panic!("Cannot decode zstd file {:?}", fin))
            }
            _ => panic!("Unsupported compress type {}", compress),
        },
    };

    buf
} */

pub fn read_file_to_end(fin: &Path, compress: Option<&str>) -> Result<Vec<u8>, Box<dyn Error>> {
    log::debug!("read file: {:?}", fin);

    let mut buf: Vec<u8> = vec![];
    match compress {
        None => {
            // -> you cannot load stdin twice [ref](https://stackoverflow.com/questions/17958571/reading-input-from-stdcin-twice)
            // -> You have to load only once at first and use them
            let _ = File::open(fin)?.read_to_end(&mut buf)?;
        }
        Some(compress) => match compress {
            "zst" => {
                buf = zstd::stream::decode_all(File::open(fin)?)?;
            }
            _ => panic!("Unsupported compress type {}", compress),
        },
    };

    Ok(buf)
}

/*
pub fn read_file_to_end(fin: &Path, compress: Option<String>) -> Vec<u8> {
    if !textfile::able_open_file(fin) {
        panic!("Cannot open {:?}", fin);
    }
    let mut buf: Vec<u8> = vec![];
    match compress {
        None => {
            File::open(fin).unwrap().read_to_end(&mut buf).unwrap();
        }
        Some(compress) => match &*compress {
            "zst" => {
                buf = zstd::stream::decode_all(File::open(fin).unwrap()).unwrap();
            }
            _ => panic!(),
        },
    };

    buf
}
*/

// use this for .fam, .bim, whose delimiter is not necessarily '\t'
pub fn compute_num_line_text(fin: &Path, compress: Option<&str>) -> Result<usize, Box<dyn Error>> {
    let buf = read_file_to_end(fin, compress)?;

    Ok(compute_num_line_buf(&buf[..]))
}

//fn compute_num_line_buf<T: std::io::Read>(buf: &mut BufReader<T>) -> usize {
pub fn compute_num_line_buf<T: std::io::Read>(buf: T) -> usize {
    let mut buf = BufReader::new(buf);

    let mut num_line: usize = 0;

    // one line of .bim and .fam should be smaller than 128.
    // otherwise automatically enlarge capacity
    let mut line = String::with_capacity(512);
    //let mut buf = BufReader::new(File::open(fin).unwrap());
    while buf.read_line(&mut line).unwrap() > 0 {
        //println!("{}", line);
        //log::debug!("{}", line);
        num_line += 1;
        line.clear(); // clear to reuse the buffer
    }
    //Some(num_line)
    num_line
}

/// read first line to get number of columns
/// this can be replaced by tsv::compute_num_column() but later
pub fn compute_num_column_text(
    fin: &Path,
    compress: Option<&str>,
) -> Result<usize, Box<dyn Error>> {
    let buf = read_file_to_end(fin, compress)?;

    Ok(compute_num_column_buf(&buf[..]))
}

//pub fn compute_num_column_buf<T: std::io::Read>(buf: &mut BufReader<T>) -> usize {
pub fn compute_num_column_buf<R: std::io::Read>(buf: R) -> usize {
    let buf = BufReader::new(buf);
    // unwrap() twice since lines() returns Option<Result<>>
    let line_first = buf.lines().next().unwrap().unwrap();
    let mut num_cols = 0;
    for _ in line_first.split_whitespace() {
        num_cols += 1;
    }
    //Some(num_cols)
    num_cols
}

// when you use buf: R,  cannot read twice
//pub fn load_table_buf<R: std::io::Read>(buf: R, header: bool) -> Vec<Vec<String>> {
pub fn load_table_buf(buf: &[u8], header: bool) -> Vec<Vec<String>> {
    //check_exist_file(fin);

    let num_col = compute_num_column_buf(buf);
    log::debug!("col: {}", num_col);
    let num_line = compute_num_line_buf(buf);
    log::debug!("line: {}", num_line);

    let mut valss: Vec<Vec<String>> = Vec::with_capacity(num_col);
    //let m: Vec<String>;
    for _ in 0..num_col {
        valss.push(Vec::with_capacity(num_line));
    }

    // one line of .bim and .fam should be smaller than 128.
    let mut line = String::with_capacity(512);
    let mut buf = BufReader::new(buf);
    //let mut buf = BufReader::new(File::open(buf).unwrap());

    if header {
        buf.read_line(&mut line).unwrap();
        line.clear();
    }

    let mut line_i = 1;
    while buf.read_line(&mut line).unwrap() > 0 {
        //for vals in valss.iter_mut() {
        //let content_i = &mut valss[line_i];
        //buf.read_line(&mut line).unwrap();
        line_i += 1;

        //log::debug!("line {}", line);

        // TODO: check len == num_col
        //let mut words = line.split_whitespace();
        //for (i, word) in .enumerate() {
        let mut col_line = 0;
        for (i, word) in line.split_whitespace().enumerate() {
            //log::debug!("i, word {},{}", i, word);
            valss[i].push(word.to_owned());
            col_line = i;
        }
        assert_eq!(
            col_line + 1,
            num_col,
            "Line {} does not have the same columns number as header: actual {} vs expected {}.",
            line_i,
            col_line + 1,
            num_col
        );

        //for (word, val) in line.split_whitespace().zip(vals.iter_mut()) {
        //    *val = word.to_owned();
        //}
        line.clear(); // clear to reuse the buffer
    }

    /*
    // do not use `cs` since not sure about delimiter
    let file = File::open(fin).unwrap();
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b" ")
        .has_headers(false)
        .from_reader(file);
    //let mut rdr = csv::Reader::from_reader(file);

    for result in rdr.records() {
        let line = result.unwrap();
        log::debug!("{:?}", line);
    }
    */

    valss
}

/// return Vec<Vec<String>> (num_columns x num_lines)
/// -> when openning file, option is natural. do unwrap() in that func.
pub fn load_table(fin: &Path, header: bool) -> Vec<Vec<String>> {
    let buf = read_file_to_end(fin, None).unwrap();

    load_table_buf(&buf[..], header)
}

// TODO: remove load_table and rename this load_table
pub fn load_table_compress(fin: &Path, header: bool, compress: Option<&str>) -> Vec<Vec<String>> {
    let buf = read_file_to_end(fin, compress).unwrap();

    load_table_buf(&buf[..], header)
}

pub fn load_table_header(fin: &Path) -> Vec<String> {
    //let num_col = compute_num_column(fin);
    //let mut vals: Vec<String> = Vec::with_capacity(num_col);

    //check_exist_file(fin);

    // one line of .bim and .fam should be smaller than 128.
    let mut line = String::with_capacity(512);
    let mut buf = BufReader::new(File::open(fin).unwrap());

    buf.read_line(&mut line).unwrap();

    let vals: Vec<String> = line
        .split_whitespace()
        .map(|word| word.to_owned())
        .collect();

    vals
}

pub fn load_table_header_buf<R: std::io::Read>(buf: R) -> Vec<String> {
    //let num_col = compute_num_column(fin);
    //let mut vals: Vec<String> = Vec::with_capacity(num_col);

    /*     if !able_open_file(fin) {
        return None;
    } */
    //check_exist_file(fin);

    // one line of .bim and .fam should be smaller than 128.
    let mut line = String::with_capacity(512);
    //let mut buf = BufReader::new(File::open(fin).unwrap());
    let mut buf = BufReader::new(buf);

    buf.read_line(&mut line).unwrap();

    let vals: Vec<String> = line
        .split_whitespace()
        .map(|word| word.to_owned())
        .collect();

    vals
}

pub fn coli_of_header_buf(buf: &[u8], col: &str) -> Option<usize> {
    let header = load_table_header_buf(buf);

    for (coli, col_header) in header.iter().enumerate() {
        if col == col_header {
            return Some(coli);
        }
    }
    None
}

pub fn load_table_col_by_header_buf(buf: &[u8], col: &str) -> Vec<String> {
    let cols = vec![col];
    load_table_cols_by_header_buf(buf, &cols)[0].clone()
}

/// orders are aligned to cols
/// if header is duplicated, the first column is used.
pub fn load_table_cols_by_header_buf(buf: &[u8], cols: &[&str]) -> Vec<Vec<String>> {
    let header = load_table_header_buf(buf);

    let mut cols_i: Vec<usize> = vec![];
    for col_use in cols.iter() {
        for (coli, col) in header.iter().enumerate() {
            if col_use == col {
                cols_i.push(coli);
                break;
            }
        }
    }

    load_table_cols_buf(buf, &cols_i, true)
}

/// orders are aligned to cols
pub fn load_table_cols_by_header(fin: &Path, cols: &[&str]) -> Vec<Vec<String>> {
    let header = load_table_header(fin);

    let mut cols_i: Vec<usize> = vec![];
    for col_use in cols.iter() {
        for (coli, col) in header.iter().enumerate() {
            if col_use == col {
                cols_i.push(coli);
                break;
            }
        }
    }

    load_table_cols(fin, &cols_i, true)
}

/// orders are aligned to cols
pub fn load_table_cols_buf(buf: &[u8], cols: &[usize], header: bool) -> Vec<Vec<String>> {
    //let buf = read_file_to_end(fin, None).unwrap();

    let num_col = compute_num_column_buf(buf);
    //let num_col = compute_num_column_buf(&buf[..]);
    //let num_col = compute_num_column_text(fin, None).unwrap();
    //let num_col = compute_num_column_text(fin, None).unwrap();
    let num_col_use = cols.len();

    if *cols.iter().max().unwrap() >= num_col {
        panic!(
            "Indicated cols exceed number of columns: {} vs  {}",
            *cols.iter().max().unwrap(),
            num_col
        );
    }

    // col in file -> col in cols
    let mut colmap: Vec<Option<usize>> = Vec::with_capacity(num_col);
    for i in 0..num_col {
        let mut colv = None;
        for (coli, col) in cols.iter().enumerate() {
            if *col == i {
                colv = Some(coli);
                break;
            }
        }
        colmap.push(colv);
    }

    let num_line = compute_num_line_buf(buf);
    //let num_line = compute_num_line_buf(&buf[..]);

    let mut valss: Vec<Vec<String>> = Vec::with_capacity(num_col_use);
    for _ in 0..num_col_use {
        valss.push(Vec::with_capacity(num_line));
    }

    // once one line exceed 512, it would automatically extended.
    let mut line = String::with_capacity(512);
    let mut bufr = BufReader::new(buf);
    //let mut bufr = BufReader::new(&buf[..]);
    //let mut buf = BufReader::new(File::open(fin).unwrap());

    if header {
        bufr.read_line(&mut line).unwrap();
        line.clear();
    }

    while bufr.read_line(&mut line).unwrap() > 0 {
        for (i, word) in line.split_whitespace().enumerate() {
            if let Some(col) = colmap[i] {
                valss[col].push(word.to_owned());
            }
        }
        line.clear(); // clear to reuse the buffer
    }

    if header {
        assert_eq!(valss[0].len(), num_line - 1);
    } else {
        assert_eq!(valss[0].len(), num_line);
    }

    valss
}

// AVOID loading from fin several times
// when fin is stdin, you cannot access twice
// first load to buf and treat them later
/// return Vec<Vec<String>> (cols x num_lines)
pub fn load_table_cols(fin: &Path, cols: &[usize], header: bool) -> Vec<Vec<String>> {
    let buf = read_file_to_end(fin, None).unwrap();

    load_table_cols_buf(&buf[..], cols, header)

    /*     let num_col = compute_num_column_buf(&buf[..]);
    //let num_col = compute_num_column_text(fin, None).unwrap();
    //let num_col = compute_num_column_text(fin, None).unwrap();
    let num_col_use = cols.len();

    if *cols.iter().max().unwrap() >= num_col {
        panic!(
            "Indicated cols exceed number of columns: {} vs  {}",
            *cols.iter().max().unwrap(),
            num_col
        );
    }

    // col in file -> col in cols
    let mut colmap: Vec<Option<usize>> = Vec::with_capacity(num_col);
    for i in 0..num_col {
        let mut colv = None;
        for (coli, col) in cols.iter().enumerate() {
            if *col == i {
                colv = Some(coli);
                break;
            }
        }
        colmap.push(colv);
    }

    let num_line = compute_num_line_buf(&buf[..]);

    let mut valss: Vec<Vec<String>> = Vec::with_capacity(num_col_use);
    for _ in 0..num_col_use {
        valss.push(Vec::with_capacity(num_line));
    }

    // once one line exceed 512, it would automatically extended.
    let mut line = String::with_capacity(512);
    let mut bufr = BufReader::new(&buf[..]);
    //let mut buf = BufReader::new(File::open(fin).unwrap());

    if header {
        bufr.read_line(&mut line).unwrap();
        line.clear();
    }

    while bufr.read_line(&mut line).unwrap() > 0 {
        for (i, word) in line.split_whitespace().enumerate() {
            if let Some(col) = colmap[i] {
                valss[col].push(word.to_owned());
            }
        }
        line.clear(); // clear to reuse the buffer
    }

    if header {
        assert_eq!(valss[0].len(), num_line - 1);
    } else {
        assert_eq!(valss[0].len(), num_line);
    }

    Some(valss) */
}

/// return Vec<String> (num_lines)
pub fn load_table_col(fin: &Path, col: usize, header: bool) -> Vec<String> {
    let num_line = compute_num_line_text(fin, None).unwrap();
    log::debug!("num_line: {}", num_line);

    let mut vals: Vec<String> = Vec::with_capacity(num_line);
    //log::debug!("vals len,cap; {},{}", vals.len(), vals.capacity());

    // once one line exceed 512, it would automatically extended.
    let mut line = String::with_capacity(512);
    let mut buf = BufReader::new(File::open(fin).unwrap());
    //log::debug!("ab");

    if header {
        buf.read_line(&mut line).unwrap();
        line.clear();
    }

    //let mut line_i: usize = 0;
    while buf.read_line(&mut line).unwrap() > 0 {
        //log::debug!("line {}", line);
        for (i, word) in line.split_whitespace().enumerate() {
            //log::debug!("word {}", word);
            if i == col {
                //log::debug!("word to push {}", word);
                vals.push(word.to_owned());
                break;
            }
        }
        //line_i += 1;
        line.clear(); // clear to reuse the buffer
    }

    if header {
        assert_eq!(vals.len(), num_line - 1);
    } else {
        assert_eq!(vals.len(), num_line);
    }

    vals
}

// as fast as load_byte2 but complecated
#[allow(dead_code)]
fn load_byte1(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    //let bed_size = calculate_bed_size(n, m);
    let v_size = crate::io_genot::calculate_bed_size_genotype(m, n);
    log::debug!("{}", v_size);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }
    let mut reader = BufReader::new(File::open(fin).unwrap());
    let n: usize = 100_000_000;
    //let n: usize = 1_000_000;
    let mut buf: Vec<u8> = Vec::with_capacity(n);

    // first for 3
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    log::debug!("{:?}", buf);

    unsafe {
        buf.set_len(n);
    }
    let mut loop_times = 0;
    // if we know the size is 9, then
    for i in 0..((v_size - 1) / n) {
        reader.read_exact(&mut buf).unwrap();
        // https://stackoverflow.com/questions/28219231/how-to-idiomatically-copy-a-slice
        v[n * i..n * (i + 1)].copy_from_slice(&buf);
        loop_times += 1;
    }
    log::debug!("loop times{}", loop_times);
    // for remaining
    let n_remain: usize = v_size - (v_size - 1) / n * n;
    unsafe {
        buf.set_len(n_remain);
    }
    reader.read_exact(&mut buf).unwrap();
    v[(v_size - 1) / n * n..v_size].copy_from_slice(&buf);
    v
}

// fast and simple
/// Use this function to load bytes.
#[allow(dead_code)]
fn load_byte2(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    // you can use "seek"
    // use this!!
    //let bed_size = calculate_bed_size(n, m);
    let v_size = crate::io_genot::calculate_bed_size_genotype(m, n);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }
    log::debug!("v_size: {}", v_size);

    let mut reader = BufReader::new(File::open(fin).unwrap());
    let n: usize = 100_000_000;
    //let n: usize = 1_000_000;
    let mut buf: Vec<u8> = Vec::with_capacity(n);

    // first for 3
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    log::debug!("{:?}", buf);

    unsafe {
        buf.set_len(n);
    }

    let mut loop_times = 0;
    let mut i_ptr: usize = 0;
    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            n => {
                log::debug!("i_ptr,n:{},{}", i_ptr, n);
                let buf = &buf[..n];
                v[i_ptr..i_ptr + n].copy_from_slice(&buf);
                i_ptr = i_ptr + n;
            }
        }
        loop_times += 1;
    }
    log::debug!("loop times{}", loop_times);
    v
}

#[allow(dead_code)]
fn load_byte3(fin: &str, n: usize, m: usize) -> Vec<u8> {
    let mut file = File::open(fin).unwrap();
    // this is ok but not sure fast or slow
    let mut buf: Vec<u8> = Vec::new();
    // meaningless, lastly cap is 32
    //let n = 9;
    //let mut buf: Vec<u8> = Vec::with_capacity(n + 2);
    // this wont work
    //unsafe {
    //    buf.set_len(n);
    //}
    log::debug!("len: {}", buf.len());
    log::debug!("cap: {}", buf.capacity());
    let _ = file.read_to_end(&mut buf).unwrap();
    log::debug!("len: {}", buf.len());
    // larger than cap
    log::debug!("cap: {}", buf.capacity());
    log::debug!("buf[0]: {}", buf[0]);
    //log::debug!("{:?}", buf);

    //let v_size = calculate_bed_genotype_size(n, m);
    //tmp
    let v_size = crate::io_genot::calculate_bed_size(n, m);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }

    (&mut v).copy_from_slice(&buf);
    v
}

pub fn run_byte(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    // fast but complecated
    log::debug!("way 1");
    let istart = Instant::now();
    let v = load_byte1(fin, n, m);
    log::debug!("load_byte1: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);

    // as fast as way 1
    // also easy to write
    log::debug!("way 2");
    let istart = Instant::now();
    let v = load_byte2(fin, n, m);
    log::debug!("load_byte2: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);

    /*
    // seems slow -> could be because file size was larger than mem size?
    log::debug!("way 3");
    let istart = Instant::now();
    let v = load_byte3(fin, n, m);
    log::debug!("load_byte3: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);
    */

    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{io_genot, GenotFormat};
    use std::path::PathBuf;

    #[test]
    fn test_isin_string() {
        let line = String::from("#ID\tABC\tdef");
        let col = String::from("ID");
        assert_eq!(isin_string(&col, &line), true);

        let col = String::from("IDA");
        assert_eq!(isin_string(&col, &line), false);
    }

    #[test]
    fn test_load_first_line() {
        let content = "abc\tdef\n1\trs5\n";
        //let mut buf=BufReader::new(content.as_bytes());
        let content = content.as_bytes();
        let content = load_first_line(content);
        assert_eq!(content, "abc\tdef");
    }

    #[test]
    fn test_read_file_to_end() {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2Vzs;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        let buf = read_file_to_end(&fin_bim, Some("zst")).unwrap();

        let buf = BufReader::new(&buf[..]);

        assert_eq!(
            buf.lines().next().unwrap().unwrap(),
            "#CHROM\tPOS\tID\tREF\tALT"
        )
    }

    #[test]
    #[should_panic]
    fn test_read_file_to_end_panic() {
        let fin = PathBuf::from("notexist.txt");
        let buf = read_file_to_end(&fin, None).unwrap();
    }

    #[test]
    fn test_compute_num_line() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        //let fin_bim = fin.clone() + ".bim";
        let m = compute_num_line_text(&fin_bim, None).unwrap();
        assert_eq!(m, 3);

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        //let fin_fam = fin.clone() + ".fam";
        let n = compute_num_line_text(&fin_fam, None).unwrap();
        assert_eq!(n, 5);
    }

    #[test]
    fn test_compute_num_line_compress() {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2Vzs;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        let m = compute_num_line_text(&fin_bim, Some("zst")).unwrap();
        assert_eq!(m, 4);

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        let n = compute_num_line_text(&fin_fam, None).unwrap();
        assert_eq!(n, 11);
    }

    #[test]
    fn test_compute_num_column() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        //let fin_bim = fin.clone() + ".bim";
        let col: usize = compute_num_column_text(&fin_bim, None).unwrap();
        assert_eq!(col, 6);

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        //let fin_fam = fin.clone() + ".fam";
        let col = compute_num_column_text(&fin_fam, None).unwrap();
        assert_eq!(col, 6);
    }

    #[test]
    fn test_compute_num_column_compress() {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2Vzs;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        let m = compute_num_column_text(&fin_bim, Some("zst")).unwrap();
        assert_eq!(m, 5);

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        let n = compute_num_column_text(&fin_fam, None).unwrap();
        assert_eq!(n, 4);
    }

    #[test]
    fn test_load_tsv() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        // let fin_bim = fin.clone() + ".bim";
        let vs = load_table(&fin_bim, false);
        assert_eq!(vs.len(), compute_num_column_text(&fin_bim, None).unwrap());
        assert_eq!(vs[0].len(), compute_num_line_text(&fin_bim, None).unwrap());
        assert_eq!(vs[0][0], "1");
        assert_eq!(vs[1][0], "rs1");
        assert_eq!(vs[5][0], "C");
        assert_eq!(vs[5][2], "C");
    }

    /// test with header
    #[test]
    fn test_load_tsv_2() {
        let fin_cov = PathBuf::from("../../test/data/toy1/genot.cov");

        let vs = load_table(&fin_cov, true);
        assert_eq!(vs.len(), compute_num_column_text(&fin_cov, None).unwrap());
        assert_eq!(
            vs[0].len(),
            compute_num_line_text(&fin_cov, None).unwrap() - 1
        );
        assert_eq!(vs[0][0], "1");
        assert_eq!(vs[1][0], "1");
        assert_eq!(vs[2][0], "1");
        assert_eq!(vs[3][0], "58");
        assert_eq!(vs[0][4], "5");
    }

    #[test]
    fn test_load_tsv_col() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        //let fin_fam = fin.clone() + ".fam";
        let v = load_table_col(&fin_fam, 4, false);
        assert_eq!(v.len(), compute_num_line_text(&fin_fam, None).unwrap());
        assert_eq!(v[0], "1");
        assert_eq!(v[4], "2");
    }

    #[test]
    fn test_load_byte2() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;

        let fin_fam = io_genot::fname_plinks_sample(&fin, gfmt, None);
        //let fin_fam = fin.clone() + ".fam";
        let n: usize = compute_num_line_text(&fin_fam, None).unwrap();

        let fin_bim = io_genot::fname_plinks_snv(&fin, gfmt, None);
        //let fin_bim = fin.clone() + ".bim";
        let m: usize = compute_num_line_text(&fin_bim, None).unwrap();

        let fin_bed = io_genot::fname_plinks_genot(&fin, gfmt, None);
        //let fin_bed = fin.clone() + ".bed";
        load_byte2(&fin_bed, n, m);
    }
}
