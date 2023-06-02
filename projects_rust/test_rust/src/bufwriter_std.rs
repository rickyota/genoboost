//! Test reader, writer
//!  https://stackoverflow.com/questions/28370126/how-can-i-test-stdin-and-stdout
//!
//!

use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::io::Cursor;
use std::io::{prelude::*, BufReader, Read, SeekFrom};
use std::io::{BufRead, BufWriter, Write};

fn file_rwc() -> File {
    let fout = "./test.wgt";

    match OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open(&fout)
    {
        Ok(file) => file,
        Err(_) => panic!("file does not exist: {}", &fout),
    }
}

fn file_rac() -> File {
    let fout = "./test.wgt";

    match OpenOptions::new()
        .read(true)
        .append(true)
        .create(true)
        .open(&fout)
    {
        Ok(file) => file,
        Err(_) => panic!("file does not exist: {}", &fout),
    }
}

fn bufwrite_file() -> BufWriter<File> {
    let fout = "./test.wgt";

    let file = match File::create(&fout) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {}",
            &fout
        ),
    };

    let mut writer = BufWriter::new(file);
    let str = "new";
    writer.write(str.as_bytes()).unwrap();
    writer
}

fn bufwrite_file_append() -> BufWriter<File> {
    /*
    let fout = "./test.wgt";

    let file = match OpenOptions::new().append(true).open(&fout) {
        Ok(file) => file,
        Err(_) => panic!("file does not exist: {}", &fout),
    };
    */

    let file = file_rac();

    let mut writer = BufWriter::new(file);
    let str = "append";
    writer.write(str.as_bytes()).unwrap();
    writer
}

fn bufread_from() -> BufReader<File> {
    let file = file_rwc();

    let buf = BufReader::new(file);

    buf
}

fn bufread_from_append() -> BufReader<File> {
    let file = file_rac();

    let buf = BufReader::new(file);

    buf
}

fn bufwrite_stdout() -> BufWriter<io::Stdout> {
    let mut writer = BufWriter::new(io::stdout());
    writer
}

fn write_both<W: std::io::Write>(writer: &mut BufWriter<W>, str: &str) {
    //let str = "def\n";
    writer.write(str.as_bytes()).unwrap();
}

fn write_both2<W: Write>(writer: &mut W, str: &str) {
    //let str = "def\n";
    writer.write(str.as_bytes()).unwrap();
}

// &[u8] doesn't have Seek
//fn read_both2<R: BufRead>(mut reader: R, s: &mut [u8]) {
//fn read_both2<R: BufRead + std::io::Seek>(reader: &mut R, s: &mut [u8]) {
fn read_both2<R: BufRead>(reader: &mut R, s: &mut [u8]) {
    //reader.seek(SeekFrom::Start(3)).unwrap();
    //reader.consume(1);
    reader.consume(1);
    reader.read(s).unwrap();
    println!("read 1 {:?}", s);
    reader.read(s).unwrap();
    println!("read 2 {:?}", s);
}

fn read_both2_seek<R: BufRead + Seek>(reader: &mut R, s: &mut [u8]) {
    reader.seek(SeekFrom::Start(3)).unwrap();
    //reader.consume(1);
    //reader.consume(1);
    reader.read(s).unwrap();
    println!("read 1 {:?}", s);
    reader.read(s).unwrap();
    println!("read 2 {:?}", s);
}

//  https://stackoverflow.com/questions/28370126/how-can-i-test-stdin-and-stdout
fn prompt<R, W>(mut reader: R, mut writer: W, question: &str) -> String
where
    R: BufRead,
    W: Write,
{
    write!(&mut writer, "{}", question).expect("Unable to write");
    let mut s = String::new();
    reader.read_line(&mut s).expect("Unable to read");
    s
}

pub fn test() {
    // this is necessary for append file
    /*
    {
        // this create ./test.wgt
        let mut writer = bufwrite_file();
        let str = "def\n";
        writer.write(str.as_bytes()).unwrap();
        write_both(&mut writer, "jkl\n");
        write_both2(&mut writer, "abcd\n");
    }
    */

    {
        // append
        println!("append");
        let mut writer_append = bufwrite_file_append();
        //let str = "appenddef\n";
        //writer.write(str.as_bytes()).unwrap();
        //write_both(&mut writer, "jkl\n");
        //write_both2(&mut writer, "abcd\n");
    }

    // is it able to use same File to read?
    let mut buf_read = bufread_from();
    //let mut buf_read = bufread_from_append();
    let mut line = String::new();
    buf_read.read_line(&mut line).unwrap();
    println!("bufread {:?}", line);

    /*
    // write
    let mut writer_stdout = bufwrite_stdout();
    let str = "ghi\n";
    writer_stdout.write(str.as_bytes()).unwrap();
    write_both(&mut writer_stdout, "mno\n");
    write_both2(&mut writer_stdout, "abc\n");

    // test for BufReader
    //let reader = vec![0u8, 1, 2];
    //let reader = b"I'm George";
    let reader = vec![0u8, 1, 2, 3, 4, 5, 6];
    let mut buf = vec![0u8; 3];
    // &[u8] implements BufRead so force casting
    read_both2(&mut reader.as_slice(), &mut buf);
    //read_both2(&mut (&reader[..]), &mut buf);
    println!("read bufreader {:?}", buf);

    let reader = vec![0u8, 1, 2, 3, 4, 5, 6];
    let mut buf = vec![0u8; 3];
    let mut cur = Cursor::new(reader.as_slice());
    read_both2_seek(&mut cur, &mut buf);
    println!("read bufreader {:?}", buf);

    /*     let input = b"I'm George";
    let mut output = Vec::new();
    // What is &v[..]?
    // &[u8] implements BufRead so force casting
    let answer = prompt(&input[..], &mut output, "Who goes there?");
    println!("answer {:?}", answer); */
    */
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_out() {
        let mut writer_stdout = bufwrite_stdout();
        write_both(&mut writer_stdout, "pqr\n");
    }
}
