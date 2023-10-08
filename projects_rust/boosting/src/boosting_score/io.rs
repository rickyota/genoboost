use genetics::samples::prelude::*;
use genetics::score as gscore;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn get_dir_score_cv(dout: &Path, cvi: usize) -> PathBuf {
    let mut d = dout.to_owned();
    let dpara = String::from("cv-") + &cvi.to_string();
    d.push(dpara);
    d
}

fn create_dir(dout: &Path) {
    fs::create_dir_all(&dout).unwrap();
}
//pub fn create_dir(dout: &Path) {
//    fs::create_dir_all(&dout).unwrap();
//}

// TODO: make nocov, is_nsnvs into enum
// call with iter and score
pub fn fname_score_concat_createdir(
    dout: &Path,
    //iter: usize,
    nocov: bool,
    is_nsnv: bool,
    integrate: bool,
) -> PathBuf {
    let f = if integrate {
        fname_score_integrate(dout, nocov)
    } else {
        fname_score_concat(dout, nocov, is_nsnv)
        //fname_score(dout, iter, nocov, is_nsnv)
    };
    let d = f.parent().unwrap().to_path_buf();
    create_dir(&d);
    f
}

// // TODO: make nocov, is_nsnvs into enum
// // call with iter and score
// pub fn fname_score_createdir(
//     dout: &Path,
//     iter: usize,
//     nocov: bool,
//     is_nsnv: bool,
//     integrate: bool,
// ) -> PathBuf {
//     let f = if integrate {
//         fname_score_integrate(dout, nocov)
//     } else {
//         fname_score(dout, iter, nocov, is_nsnv)
//     };
//     let d = f.parent().unwrap().to_path_buf();
//     create_dir(&d);
//     f
// }

fn fname_score_concat(dout: &Path, nocov: bool, is_nsnv: bool) -> PathBuf {
    let para = if is_nsnv {
        "n".to_string()
    } else {
        "iter".to_string()
    };

    let mut d = dout.to_owned();
    let fname = if nocov {
        "boosting_".to_string() + &para + ".score"
    } else {
        "boosting_".to_string() + &para + ".scorecov"
    };
    d.push(fname);
    d
}

// fn fname_score(dout: &Path, iter: usize, nocov: bool, is_nsnv: bool) -> PathBuf {
//     if is_nsnv {
//         fname_score_nsnv(dout, iter, nocov)
//     } else {
//         fname_score_iter(dout, iter, nocov)
//     }
// }

fn fname_score_integrate(dout: &Path, nocov: bool) -> PathBuf {
    let mut d = dout.to_owned();
    let fname = if nocov {
        "boosting.score"
        //"score.tsv"
    } else {
        "boosting.scorecov"
        //"score.withcov.tsv"
    };
    d.push(fname);
    d
}

// fn fname_score_iter(dout: &Path, iter: usize, nocov: bool) -> PathBuf {
//     let mut d = dout.to_owned();
//     let dscore = if nocov {
//         "score.iter"
//     } else {
//         "score.withcov.iter"
//     };
//     d.push(dscore);
//     let fname = "iter-".to_string() + &iter.to_string() + ".score";
//     d.push(fname);
//     d
// }

// fn fname_score_nsnv(dout: &Path, nsnv: usize, nocov: bool) -> PathBuf {
//     let mut d = dout.to_owned();
//     let dscore = if nocov {
//         "score.nsnv"
//     } else {
//         "score.withcov.nsnv"
//     };
//     d.push(dscore);
//     let fname = "nsnv-".to_string() + &nsnv.to_string() + ".score";
//     d.push(fname);
//     d
// }

pub fn write_scores_concat(
    fout: &Path,
    score_paras: &[Vec<f64>],
    paras: &[String],
    //scores: &[f64],
    samples_id: &[String],
    is_nsnv: bool,
    integrate: bool,
) {
    if integrate {
        assert_eq!(score_paras.len(), 1);
        gscore::write_scores_nopheno(fout, &score_paras[0], samples_id)
    } else {
        if is_nsnv {
            gscore::write_scores_paras_nopheno(fout, score_paras, "n", paras, samples_id)
        } else {
            gscore::write_scores_paras_nopheno(fout, score_paras, "iter", paras, samples_id)
        }
    }
}

// // moved to genetics::score
// //pub fn write_scores(fout: &Path, scores: &[f64], phe: &Phe, samples_id: &[(String, String)]) {
// //pub fn write_scores(fout: &Path, scores: &[f64], phe: &Phe, samples_id: &[String]) {
// pub fn write_scores(fout: &Path, scores: &[f64], samples_id: &[String]) {
//     let file = match File::create(&fout) {
//         Ok(file) => file,
//         Err(_) => panic!(
//             "Cannot create file, possibly directory does not exist: {:?}",
//             &fout
//         ),
//     };

//     let mut writer = BufWriter::new(file);
//     // assume word count of one line is 30
//     // no problem when it exceeds
//     let capacity = 30 * scores.len();
//     let mut score_string = String::with_capacity(capacity);
//     score_string.push_str("iid\tscore\n");
//     //score_string.push_str("fid\tiid\tphe\trs\n");

//     for ni in 0..scores.len() {
//         score_string.push_str(&samples_id[ni]);
//         score_string.push_str("\t");
//         //score_string.push_str(&samples_id[ni]);
//         //score_string.push_str("\t");
//         //score_string.push_str(&(phe.get_unchecked(ni) as u8).to_string());
//         //score_string.push_str("\t");
//         score_string.push_str(&format!("{:.5}\n", scores[ni]));
//     }

//     writer.write(score_string.as_bytes()).unwrap();

//     //for ni in
//     // use .concat(), .join()?
//     // https://users.rust-lang.org/t/fast-string-concatenation/4425/5
//     // -> .push_str seems fastest
//     // -> could be because use with_capacity beforehand
// }

/* // moved to genetics::score
//pub fn write_scores(fout: &Path, scores: &[f64], phe: &Phe, samples_id: &[(String, String)]) {
pub fn write_scores(fout: &Path, scores: &[f64], phe: &Phe, samples_id: &[String]) {
    let file = match File::create(&fout) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {:?}",
            &fout
        ),
    };

    let mut writer = BufWriter::new(file);
    // assume word count of one line is 30
    // no problem when it exceeds
    let capacity = 30 * scores.len();
    let mut score_string = String::with_capacity(capacity);
    score_string.push_str("fid\tiid\tphe\trs\n");

    for ni in 0..scores.len() {
        score_string.push_str(&samples_id[ni]);
        score_string.push_str("\t");
        score_string.push_str(&samples_id[ni]);
        score_string.push_str("\t");
        score_string.push_str(&(phe.get_unchecked(ni) as u8).to_string());
        score_string.push_str("\t");
        score_string.push_str(&format!("{:.5}\n", scores[ni]));
    }

    writer.write(score_string.as_bytes()).unwrap();

    //for ni in
    // use .concat(), .join()?
    // https://users.rust-lang.org/t/fast-string-concatenation/4425/5
    // -> .push_str seems fastest
    // -> could be because use with_capacity beforehand
} */

#[cfg(test)]
mod tests {
    //use super::*;
}
