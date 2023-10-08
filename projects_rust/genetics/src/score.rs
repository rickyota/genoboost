use crate::genot::prelude::*;
use crate::wgt::{Coef, WgtKind};
use crate::Covs;
use crate::Dataset;
use crate::Wgt;
use crate::Wgts;
//use core::panic;
use std::hash::Hash;
use std::collections::HashSet;
//use rayon::iter::ParallelBridge;
//use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
//use crate::samples::phe::Phe;
use crate::samples::prelude::*;

// (s0, s1, s2)
fn scores_add_coef(scores: &mut [f64], score_wgt: (f64, f64, f64), genot_mi: &GenotSnvRef) {
    let (s0, s1, s2) = score_wgt;
    scores
        .iter_mut()
        .zip(genot_mi.iter())
        .for_each(|(score, val)| {
            let score_add = match val {
                0 => s0,
                1 => s1,
                2 => s2,
                // panic for NA
                //3 => sm,
                _ => {
                    panic!("Wrong genotype. Possibly NA for linear model.")
                }
            };
            *score += score_add;
        })
}

// DO NOT use rayon here
// since very slow since n is small
pub fn add_score(scores: &mut [f64], wgt: &Wgt, genot: &Genot, covs: Option<&Covs>) {
    //log::debug!("wgt {:?}", wgt);
    match wgt.wgt().kind() {
        WgtKind::Snv(_, _, mi) => {
            //log::debug!("mi {}", mi.unwrap());
            // This should never panic.
            match mi {
                None => return,
                Some(mi) => {
                    let genot_mi = genot.to_genot_snv(*mi);
                    match wgt.wgt().model().coef() {
                        Coef::Linear(alpha_ti) => {
                            let s2 = alpha_ti * 2.0;
                            let s1 = alpha_ti;
                            let s0 = 0.0;
                            scores_add_coef(scores, (s0, s1, s2), &genot_mi);
                            // not checked but this should work
                            //scores
                            //    .iter_mut()
                            //    .zip(genot_mi.iter())
                            //    .for_each(|(score, val)| {
                            //        let score_add = match val {
                            //            2 => s2,
                            //            1 => s1,
                            //            0 => s0,
                            //            // panic for NA
                            //            //3 => sm,
                            //            _ => {
                            //                panic!("Wrong genotype. Possibly NA for linear model.")
                            //            }
                            //        };
                            //        *score += score_add;
                            //    })
                        }
                        Coef::Score3(score_wgt) => {
                            //let s2 = alpha_ti * 2.0;
                            //let s1 = alpha_ti;
                            //let s0 = 0.0;
                            scores_add_coef(scores, score_wgt, &genot_mi);
                        }
                        _ => unimplemented!(),
                    }
                }
            }

            //let genot_mi = genot.to_genot_snv(mi.unwrap());
            //match wgt.wgt().model().coef() {
            //    Coef::Linear(alpha_ti) => {
            //        let s2 = alpha_ti * 2.0;
            //        let s1 = alpha_ti;
            //        let s0 = 0.0;
            //        // not checked but this should work
            //        scores
            //            .iter_mut()
            //            .zip(genot_mi.iter())
            //            .for_each(|(score, val)| {
            //                let score_add = match val {
            //                    2 => s2,
            //                    1 => s1,
            //                    0 => s0,
            //                    // panic for NA
            //                    //3 => sm,
            //                    _ => panic!("Wrong genotype. Possibly NA for linear model."),
            //                };
            //                *score += score_add;
            //            })
            //    }
            //    _ => unimplemented!(),
            //}
        }
        WgtKind::Cov(_) => {
            let cov_name = wgt.wgt().kind().cov().name();

            // TODO: clean const
            if cov_name == "const" {
                match wgt.wgt().model().coef() {
                    Coef::Linear(alpha) => {
                        scores.iter_mut().for_each(|score| {
                            *score += alpha;
                        });
                        //scores.par_iter_mut().for_each(|score| {
                        //    *score += alpha;
                        //});

                        //for score in scores.iter_mut() {
                        //    *score += alpha;
                        //}
                    }
                    Coef::Binary((const_t, alpha)) => {
                        scores.iter_mut().for_each(|score| {
                            *score += alpha + const_t;
                        });
                        //scores.par_iter_mut().for_each(|score| {
                        //    *score += alpha + const_t;
                        //});
                        //for score in scores.iter_mut() {
                        //    *score += alpha + const_t;
                        //}
                    }
                    _ => unimplemented!(),
                }
            } else {
                // cov_name=const will fail
                let cov_vals = covs.unwrap().vals_id(cov_name);
                //let cov_vals = covs.vals_id(cov_name);
                //let cov_vals = cov_wgt.var().vals();
                match wgt.wgt().model().coef() {
                    Coef::Linear(alpha) => {
                        scores
                            .iter_mut()
                            .zip(cov_vals.iter())
                            .for_each(|(score, val)| {
                                *score += alpha * val;
                            });
                        //scores
                        //    .iter_mut()
                        //    .zip(cov_vals.iter())
                        //    .par_bridge()
                        //    .for_each(|(score, val)| {
                        //        *score += alpha * val;
                        //    });

                        //for (score, val) in scores.iter_mut().zip(cov_vals.iter()) {
                        //    *score += alpha * val;
                        //}
                    }
                    _ => unimplemented!(),
                }
            }
        }
    }
}

#[allow(dead_code)]
pub fn score(fout_score: &Path, wgts: &Wgts, dataset: &Dataset, samples_id: &[String]) {
    let n = dataset.genot().n();
    let genot = dataset.genot();
    let phe = dataset.samples().phe_unwrap();

    let mut scores = vec![0.0f64; n];
    for (wgt_i, wgt) in wgts.wgts().iter().enumerate() {
        if (wgt_i > 0) & (wgt_i % 100_000 == 0) {
            log::debug!("{:?}-th SNV", wgt_i);
        }

        add_score(&mut scores, wgt, genot, None);
    }

    write_scores(&fout_score, &scores, phe, samples_id);
    //io::write_scores(&fout_iteration, &scores, phe, samples_id);
}

pub fn score_nowrite(wgts: &Wgts, dataset: &Dataset) -> Vec<f64> {
    let n = dataset.genot().n();
    let genot = dataset.genot();
    //let phe = dataset.samples().phe();

    let mut scores = vec![0.0f64; n];
    for (wgt_i, wgt) in wgts.wgts().iter().enumerate() {
        if (wgt_i > 0) & (wgt_i % 100_000 == 0) {
            log::debug!("{:?}-th SNV", wgt_i);
        }

        add_score(&mut scores, wgt, genot, None);
    }

    //write_scores(&fout_score, &scores, phe, samples_id);
    scores
}

/* pub fn score(fout_score: &Path, wgts: &Wgts, dataset: &Dataset, samples_id: &[(String, String)]) {
    let n = dataset.genot().n();
    let genot = dataset.genot();
    let phe = dataset.samples().phe();

    let mut scores = vec![0.0f64; n];
    for (wgt_i, wgt) in wgts.wgts().iter().enumerate() {
        if (wgt_i>0) & (wgt_i % 100_000 == 0 ){
            log::debug!("{:?}-th SNV", wgt_i);
        }

        add_score(&mut scores, wgt, genot, None);
    }

    write_scores(&fout_score, &scores, phe, samples_id);
    //io::write_scores(&fout_iteration, &scores, phe, samples_id);
} */

//fn write_scores(fout: &Path, scores: &[f64], phe: &Phe, samples_id: &[(String, String)]) {
#[allow(dead_code)]
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
    let mut str = String::with_capacity(capacity);
    str.push_str("fid\tiid\tphe\trs\n");

    for ni in 0..scores.len() {
        str.push_str(&samples_id[ni]);
        str.push_str("\t");
        str.push_str(&samples_id[ni]);
        str.push_str("\t");
        str.push_str(&(phe.get_unchecked(ni) as u8).to_string());
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

pub fn write_scores_nopheno(fout: &Path, scores: &[f64], samples_id: &[String]) {
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
    let mut str = String::with_capacity(capacity);
    //str.push_str("fid\tiid\tphe\trs\n");
    //str.push_str("fid\tiid\tscore\n");
    str.push_str("id\tscore\n");

    for ni in 0..scores.len() {
        str.push_str(&samples_id[ni]);
        //str.push_str("\t");
        //str.push_str(&samples_id[ni]);
        //str.push_str("\t");
        //str.push_str(&(phe.get_unchecked(ni) as u8).to_string());
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

fn is_unique<T>(iter: T) -> bool
where
    T: IntoIterator,
    T::Item: Eq + Hash,
{
    let mut uniq = HashSet::new();
    iter.into_iter().all(move |x| uniq.insert(x))
}

pub fn write_scores_paras_nopheno(
    fout: &Path,
    score_paras: &[Vec<f64>],
    concat_para: &str,
    paras: &[String],
    samples_id: &[String],
) {
    if score_paras.len() != paras.len() {
        panic!("score_paras.len() != paras.len()");
    }

    if !is_unique(paras) {
        panic!("paras should be unique.");
    }

    let file = match File::create(&fout) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {:?}",
            &fout
        ),
    };

    let sample_n = score_paras[0].len();

    let mut writer = BufWriter::new(file);
    // assume word count of one line is 30
    // no problem when it exceeds
    let capacity = 30 * score_paras.len() * sample_n;
    let mut str = String::with_capacity(capacity);

    let header: String = "id\t".to_string()
        + &paras
            .iter()
            .map(|x| "score_".to_string() + concat_para + "-" + x)
            .collect::<Vec<String>>()
            .join("\t")
        + "\n";

    str.push_str(&header);
    //str.push_str("id\tscore\n");

    for ni in 0..sample_n {
        str.push_str(&samples_id[ni]);
        //str.push_str("\t");
        for scores in score_paras {
            str.push_str(&format!("\t{:.5}", scores[ni]));
        }
        str.push_str("\n");
    }

    writer.write(str.as_bytes()).unwrap();

    //for ni in
    // use .concat(), .join()?
    // https://users.rust-lang.org/t/fast-string-concatenation/4425/5
    // -> .push_str seems fastest
    // -> could be because use with_capacity beforehand
}
