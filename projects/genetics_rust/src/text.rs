//! TODO: when fin does not exist, raise error or return Option<usize>
//! later might be better
// TODO: PathBuf -> T: AsRef<Path>

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::time::Instant;

/*
fn bufwrite(fout: &str) -> BufWriter<File> {
    let file = match File::create(&fout) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {}",
            &fout
        ),
    };
    BufWriter::new(file)
    //let mut writer = BufWriter::new(file);
}
 */

//pub fn exist_file_pathbuf(fin: &Path) -> bool {
//    fin.exists()
//}

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
    match File::open(fin) {
        Ok(_) => true,
        Err(_) => false,
    }
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

// or return Option<usize> when fin does not exist
pub fn compute_num_line(fin: &Path) -> Option<usize> {
    if !able_open_file(fin) {
        return None;
    }

    let mut num_line: usize = 0;

    // one line of .bim and .fam should be smaller than 128.
    // otherwise automatically enlarge capacity
    let mut line = String::with_capacity(512);
    let mut buf = BufReader::new(File::open(fin).unwrap());
    while buf.read_line(&mut line).unwrap() > 0 {
        //println!("{}", line);
        num_line += 1;
        line.clear(); // clear to reuse the buffer
    }
    Some(num_line)
}

/// read first line to get number of columns
pub fn compute_num_column(fin: &Path) -> Option<usize> {
    if !able_open_file(fin) {
        return None;
    }

    let buf = BufReader::new(File::open(fin).unwrap());
    // why shoud `.unwrap()` twice??
    // -> lines() returns Option<Result<>>
    let line_first = buf.lines().next().unwrap().unwrap();
    let mut num_cols = 0;
    for _ in line_first.split_whitespace() {
        num_cols += 1;
    }
    Some(num_cols)
}

// TODO: why this shape?
/// return Vec<Vec<String>> (num_columns x num_lines)
pub fn load_table(fin: &Path, header: bool) -> Option<Vec<Vec<String>>> {
    if !able_open_file(fin) {
        return None;
    }
    //check_exist_file(fin);

    let num_col = compute_num_column(fin).unwrap();
    println!("col: {}", num_col);
    let num_line = compute_num_line(fin).unwrap();
    println!("line: {}", num_line);

    let mut valss: Vec<Vec<String>> = Vec::with_capacity(num_col);
    //let m: Vec<String>;
    for _ in 0..num_col {
        valss.push(Vec::with_capacity(num_line));
    }

    // one line of .bim and .fam should be smaller than 128.
    let mut line = String::with_capacity(512);
    let mut buf = BufReader::new(File::open(fin).unwrap());

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

        //println!("line {}", line);

        // TODO: check len == num_col
        //let mut words = line.split_whitespace();
        //for (i, word) in .enumerate() {
        let mut col_line = 0;
        for (i, word) in line.split_whitespace().enumerate() {
            //println!("i, word {},{}", i, word);
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
        println!("{:?}", line);
    }
    */

    Some(valss)
}

pub fn load_table_header(fin: &Path) -> Option<Vec<String>> {
    //let num_col = compute_num_column(fin);
    //let mut vals: Vec<String> = Vec::with_capacity(num_col);

    if !able_open_file(fin) {
        return None;
    }
    //check_exist_file(fin);

    // one line of .bim and .fam should be smaller than 128.
    let mut line = String::with_capacity(512);
    let mut buf = BufReader::new(File::open(fin).unwrap());

    buf.read_line(&mut line).unwrap();

    let vals: Vec<String> = line
        .split_whitespace()
        .map(|word| word.to_owned())
        .collect();

    Some(vals)
}

/// return Vec<Vec<String>> (cols x num_lines)
pub fn load_table_cols(fin: &Path, cols: &[usize], header: bool) -> Option<Vec<Vec<String>>> {
    if !able_open_file(fin) {
        return None;
    }
    //check_exist_file(fin);

    let num_col = compute_num_column(fin).unwrap();
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

    let num_line = compute_num_line(fin).unwrap();

    let mut valss: Vec<Vec<String>> = Vec::with_capacity(num_col_use);
    for _ in 0..num_col_use {
        valss.push(Vec::with_capacity(num_line));
    }

    // once one line exceed 512, it would automatically extended.
    let mut line = String::with_capacity(512);
    let mut buf = BufReader::new(File::open(fin).unwrap());

    if header {
        buf.read_line(&mut line).unwrap();
        line.clear();
    }

    while buf.read_line(&mut line).unwrap() > 0 {
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

    Some(valss)
}

/// return Vec<String> (num_lines)
pub fn load_table_col(fin: &Path, col: usize, header: bool) -> Option<Vec<String>> {
    if !able_open_file(fin) {
        return None;
    }

    let num_line = compute_num_line(fin).unwrap();
    println!("num_line: {}", num_line);

    let mut vals: Vec<String> = Vec::with_capacity(num_line);
    //println!("vals len,cap; {},{}", vals.len(), vals.capacity());

    // once one line exceed 512, it would automatically extended.
    let mut line = String::with_capacity(512);
    let mut buf = BufReader::new(File::open(fin).unwrap());
    //println!("ab");

    if header {
        buf.read_line(&mut line).unwrap();
        line.clear();
    }

    //let mut line_i: usize = 0;
    while buf.read_line(&mut line).unwrap() > 0 {
        //println!("line {}", line);
        for (i, word) in line.split_whitespace().enumerate() {
            //println!("word {}", word);
            if i == col {
                //println!("word to push {}", word);
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

    Some(vals)
}

// as fast as load_byte2 but complecated
#[allow(dead_code)]
fn load_byte1(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    //let bed_size = calculate_bed_size(n, m);
    let v_size = super::plink::calculate_bed_size_genotype(m, n);
    println!("{}", v_size);
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
    println!("{:?}", buf);

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
    println!("loop times{}", loop_times);
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
    let v_size = super::plink::calculate_bed_size_genotype(m, n);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }
    println!("v_size: {}", v_size);

    let mut reader = BufReader::new(File::open(fin).unwrap());
    let n: usize = 100_000_000;
    //let n: usize = 1_000_000;
    let mut buf: Vec<u8> = Vec::with_capacity(n);

    // first for 3
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    println!("{:?}", buf);

    unsafe {
        buf.set_len(n);
    }

    let mut loop_times = 0;
    let mut i_ptr: usize = 0;
    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            n => {
                println!("i_ptr,n:{},{}", i_ptr, n);
                let buf = &buf[..n];
                v[i_ptr..i_ptr + n].copy_from_slice(&buf);
                i_ptr = i_ptr + n;
            }
        }
        loop_times += 1;
    }
    println!("loop times{}", loop_times);
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
    println!("len: {}", buf.len());
    println!("cap: {}", buf.capacity());
    let _ = file.read_to_end(&mut buf).unwrap();
    println!("len: {}", buf.len());
    // larger than cap
    println!("cap: {}", buf.capacity());
    println!("buf[0]: {}", buf[0]);
    //println!("{:?}", buf);

    //let v_size = calculate_bed_genotype_size(n, m);
    //tmp
    let v_size = super::plink::calculate_bed_size(n, m);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }

    (&mut v).copy_from_slice(&buf);
    v
}

pub fn run_byte(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    // fast but complecated
    println!("way 1");
    let istart = Instant::now();
    let v = load_byte1(fin, n, m);
    println!("load_byte1: {:?}", Instant::now().duration_since(istart));
    println!("{}", v[0]);

    // as fast as way 1
    // also easy to write
    println!("way 2");
    let istart = Instant::now();
    let v = load_byte2(fin, n, m);
    println!("load_byte2: {:?}", Instant::now().duration_since(istart));
    println!("{}", v[0]);

    /*
    // seems slow -> could be because file size was larger than mem size?
    println!("way 3");
    let istart = Instant::now();
    let v = load_byte3(fin, n, m);
    println!("load_byte3: {:?}", Instant::now().duration_since(istart));
    println!("{}", v[0]);
    */

    v
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use crate::plink;

    use super::*;

    #[test]
    #[should_panic]
    fn test_check_exist_file() {
        let fin = PathBuf::from("../../test/data/toy1/does_not_exist.csv");
        check_exist_file(&fin);
    }

    #[test]
    fn test_compute_num_column() {
        let fin = PathBuf::from("../../test/data/toy1/genot");

        let fin_bim = plink::fname_bim(&fin, None);
        //let fin_bim = fin.clone() + ".bim";
        let m = compute_num_line(&fin_bim).unwrap();
        assert_eq!(m, 3);

        let fin_fam = plink::fname_fam(&fin, None);
        //let fin_fam = fin.clone() + ".fam";
        let n = compute_num_line(&fin_fam).unwrap();
        assert_eq!(n, 5);
    }

    #[test]
    fn test_compute_num_line() {
        let fin = PathBuf::from("../../test/data/toy1/genot");

        let fin_bim = plink::fname_bim(&fin, None);
        //let fin_bim = fin.clone() + ".bim";
        let col: usize = compute_num_column(&fin_bim).unwrap();
        assert_eq!(col, 6);

        let fin_fam = plink::fname_fam(&fin, None);
        //let fin_fam = fin.clone() + ".fam";
        let col = compute_num_column(&fin_fam).unwrap();
        assert_eq!(col, 6);
    }

    #[test]
    fn test_load_tsv() {
        let fin = PathBuf::from("../../test/data/toy1/genot");

        let fin_bim = plink::fname_bim(&fin, None);
        // let fin_bim = fin.clone() + ".bim";
        let vs = load_table(&fin_bim, false).unwrap();
        assert_eq!(vs.len(), compute_num_column(&fin_bim).unwrap());
        assert_eq!(vs[0].len(), compute_num_line(&fin_bim).unwrap());
        assert_eq!(vs[0][0], "1");
        assert_eq!(vs[1][0], "rs1");
        assert_eq!(vs[5][0], "C");
        assert_eq!(vs[5][2], "C");
    }

    /// test with header
    #[test]
    fn test_load_tsv_2() {
        let fin_cov = PathBuf::from("../../test/data/toy1/genot.cov");

        let vs = load_table(&fin_cov, true).unwrap();
        assert_eq!(vs.len(), compute_num_column(&fin_cov).unwrap());
        assert_eq!(vs[0].len(), compute_num_line(&fin_cov).unwrap() - 1);
        assert_eq!(vs[0][0], "1");
        assert_eq!(vs[1][0], "1");
        assert_eq!(vs[2][0], "1");
        assert_eq!(vs[3][0], "58");
        assert_eq!(vs[0][4], "5");
    }

    #[test]
    fn test_load_tsv_col() {
        let fin = PathBuf::from("../../test/data/toy1/genot");

        let fin_fam = plink::fname_fam(&fin, None);
        //let fin_fam = fin.clone() + ".fam";
        let v = load_table_col(&fin_fam, 4, false).unwrap();
        assert_eq!(v.len(), compute_num_line(&fin_fam).unwrap());
        assert_eq!(v[0], "1");
        assert_eq!(v[4], "2");
    }

    #[test]
    fn test_load_byte2() {
        let fin = PathBuf::from("../../test/data/toy1/genot");

        let fin_fam = plink::fname_fam(&fin, None);
        //let fin_fam = fin.clone() + ".fam";
        let n: usize = compute_num_line(&fin_fam).unwrap();

        let fin_bim = plink::fname_bim(&fin, None);
        //let fin_bim = fin.clone() + ".bim";
        let m: usize = compute_num_line(&fin_bim).unwrap();

        let fin_bed = plink::fname_bed(&fin, None);
        //let fin_bed = fin.clone() + ".bed";
        load_byte2(&fin_bed, n, m);
    }
}
