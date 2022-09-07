use genetics::samples::prelude::*;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn create_dir(dout: &Path) {
    fs::create_dir_all(&dout).unwrap();
}
//pub fn create_dir(dout: &Path) {
//    fs::create_dir_all(&dout).unwrap();
//}

// TODO: make nocov, is_nsnvs into enum
// call with iter and score
pub fn fname_score_createdir(dout: &Path, iter: usize, nocov: bool, is_nsnv: bool) -> PathBuf {
    let f = fname_score(dout, iter, nocov, is_nsnv);
    let d = f.parent().unwrap().to_path_buf();
    create_dir(&d);
    f
}

pub fn fname_score(dout: &Path, iter: usize, nocov: bool, is_nsnv: bool) -> PathBuf {
    if is_nsnv {
        fname_score_nsnv(dout, iter, nocov)
    } else {
        fname_score_iter(dout, iter, nocov)
    }
}

pub fn fname_score_iter(dout: &Path, iter: usize, nocov: bool) -> PathBuf {
    let mut d = dout.to_owned();
    let dscore;
    if nocov {
        dscore = "score.nocov.iter";
        //fout.to_owned() + ".nocov.iter" + &iter.to_string() + ".score"
        //dout.to_owned() + ".score.nocov.iter" + &iter.to_string()
    } else {
        dscore = "score.iter";
        //dout.to_owned() + ".score.iter" + &iter.to_string()
    }

    d.push(dscore);
    let fname = "iter".to_string() + &iter.to_string() + ".score";
    d.push(fname);
    d
}

pub fn fname_score_nsnv(dout: &Path, nsnv: usize, nocov: bool) -> PathBuf {
    let mut d = dout.to_owned();
    let dscore;
    if nocov {
        dscore = "score.nocov.nsnv";
        //fout.to_owned() + ".nocov.iter" + &iter.to_string() + ".score"
        //dout.to_owned() + ".score.nocov.iter" + &iter.to_string()
    } else {
        dscore = "score.nsnv";
        //dout.to_owned() + ".score.iter" + &iter.to_string()
    }

    d.push(dscore);
    let fname = "nsnv".to_string() + &nsnv.to_string() + ".score";
    d.push(fname);
    d
    /*     if nocov {
        dout.to_owned() + ".score.nocov.nsnv" + &nsnv.to_string()
        //fout.to_owned() + ".nocov.n" + &nsnv.to_string() + ".score"
    } else {
        dout.to_owned() + ".score.nsnv" + &nsnv.to_string()
        //fout.to_owned() + ".n" + &nsnv.to_string() + ".score"
    } */
}

pub fn write_scores(dout: &Path, scores: &[f64], ys: &Phe, samples_id: &[(String, String)]) {
    let file = match File::create(&dout) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {:?}",
            &dout
        ),
    };

    let mut writer = BufWriter::new(file);
    // assume word count of one line is 30
    // no problem when it exceeds
    let capacity = 30 * scores.len();
    let mut str = String::with_capacity(capacity);
    str.push_str("fid\tiid\tphe\trs\n");

    for ni in 0..scores.len() {
        str.push_str(&samples_id[ni].0);
        str.push_str("\t");
        str.push_str(&samples_id[ni].1);
        str.push_str("\t");
        str.push_str(&(ys.get_unchecked(ni) as u8).to_string());
        str.push_str("\t");
        str.push_str(&format!("{:.5}\n", scores[ni]));
    }

    writer.write(str.as_bytes()).unwrap();

    //for ni in
    // use .concat(), .join()?
    // https://users.rust-lang.org/t/fast-string-concatenation/4425/5
    // -> .push_str seems fastest
    // -> could be because use with_capacity beforehand
}

#[cfg(test)]
mod tests {
    //use super::*;
}
