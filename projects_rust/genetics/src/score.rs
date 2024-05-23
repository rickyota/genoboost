mod calc;
mod score_struct;

//use crate::alloc;
use crate::genot::{self, prelude::*};
use crate::wgt::{Coef, WgtKind};
use crate::Wgts;
use crate::{vec, Covs};
use crate::{Dataset, WgtTrait};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub use score_struct::*;

// for rayon::reduce
pub fn calc_score_wgt<W: WgtTrait>(
    n: usize,
    wgt: &W,
    genot: Option<&Genot>,
    covs: Option<&Covs>,
    allow_nonexist_snv: bool,
    // TODO: create Option<MissingTo>
    missing_to_mode: bool,
    missing_to_mean: bool,
) -> SampleScore {
    let mut scores = SampleScore::new(n);

    add_score(
        &mut scores,
        wgt,
        genot,
        covs,
        allow_nonexist_snv,
        missing_to_mode,
        missing_to_mean,
    );

    scores
}

//pub fn calc_score_wgt<W: WgtTrait>(
//    n: usize,
//    wgt: &W,
//    genot: Option<&Genot>,
//    covs: Option<&Covs>,
//    allow_nonexist_snv: bool,
//) -> Vec<f64> {
//    //let mut scores = vec![0.0f64; n];
//    let mut scores = alloc::with_capacity_align::<f64>(n + 32);
//    scores.resize(n + 32, 0.0f64);
//    add_score(&mut scores, wgt, genot, covs, allow_nonexist_snv);
//    scores
//}

// TODO: mv to score_struct.rs
fn mode_to_score(mode: u8, scores: (f64, f64, f64)) -> f64 {
    match mode {
        0 => scores.0,
        1 => scores.1,
        2 => scores.2,
        _ => panic!("Wrong genotype."),
    }
}

// DO NOT use rayon here
// since very slow when n is small
// though nested rayon does not raise error
// TODO: add Option<Snvs> for maf
//
// should implement sth in Wgt, Model to make this simple
/// Missing values
/// - If no missing values in genot, none of maf, missing_to_mode, missing_to_mean is used
/// - If missing values exists but all Score4, the same.
///
/// For Score3, missing_to_mode is used.
/// For Linear, LinearConst, eigher missing_to_mode or missing_to_mean is used.
///
/// Does not check if missing_to_mode and missing_to_mean are both true.
pub fn add_score<W: WgtTrait>(
    scores: &mut SampleScore,
    //scores_pad: &mut [f64],
    wgt: &W,
    genot: Option<&Genot>,
    covs: Option<&Covs>,
    allow_nonexist_snv: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
) {
    match wgt.wgt().kind() {
        WgtKind::Snv(snv_wgt) => {
            let mi = snv_wgt.index();
            //log::debug!("mi {}", mi.unwrap());

            let maf = snv_wgt.maf();

            let mi = match mi {
                None => {
                    if allow_nonexist_snv {
                        return;
                    } else {
                        panic!("SNV does not exist.")
                    }
                }
                Some(mi) => mi,
            };

            let gsnv = genot.unwrap().to_genot_snv(mi);
            match wgt.wgt().model().coef() {
                Coef::Single(alpha) => {
                    scores.scores_mut().iter_mut().for_each(|score| {
                        *score += alpha;
                    });
                }
                Coef::Binary((const_ti, alpha_ti)) => {
                    let threshold = wgt.wgt().model().threshold().unwrap();
                    let score_add_high = const_ti + alpha_ti;
                    let score_add_low = const_ti - alpha_ti;

                    // using par_bridge() here is very slow
                    scores
                        .scores_mut()
                        .iter_mut()
                        .zip(gsnv.iter())
                        .for_each(|(score, val)| {
                            let score_add = if val == 3 {
                                // for Binary, missing is not allowed.
                                // (should be `const_ti`??)
                                panic!("Cannot use mising on Binary coef.");
                            } else {
                                match (val as f64) > threshold {
                                    true => score_add_high,
                                    false => score_add_low,
                                }
                            };
                            *score += score_add;
                        });
                }

                Coef::Linear(alpha_ti) => {
                    let s0 = 0.0;
                    let s1 = alpha_ti;
                    let s2 = 2.0 * alpha_ti;

                    let sm = if missing_to_mode {
                        let mode = match maf {
                            Some(maf) => genot::maf_to_mode(maf),
                            None=>panic!("MAF is missing but --missing-to-mode is used. Check column name or use --fill-missing-in-dataset to fill missing values in test dataset.")
                        };
                        //let mode = genot::maf_to_mode(maf.unwrap());
                        mode_to_score(mode, (s0, s1, s2))
                    } else if missing_to_mean {
                        match maf {
                            Some(maf) => 2.0f64 * maf * alpha_ti,
                            None=>panic!("MAF is missing but --missing-to-mean is used. Check column name or use --fill-missing-in-dataset to fill missing values in test dataset.")
                        }
                    } else {
                        f64::NAN
                    };

                    calc::scores_add_coef(scores, (s0, s1, s2, sm), &gsnv);
                    //calc::scores_add_coef_nomissing(scores, (s0, s1, s2), &gsnv);
                }

                Coef::LinearConst((const_ti, alpha_ti)) => {
                    let s0 = const_ti;
                    let s1 = const_ti + alpha_ti;
                    let s2 = const_ti + 2.0 * alpha_ti;

                    let sm = if missing_to_mode {
                        let mode = match maf {
                            Some(maf) => genot::maf_to_mode(maf),
                            None=>panic!("MAF is missing but --missing-to-mode is used. Check column name or use --fill-missing-in-dataset to fill missing values in test dataset.")
                        };
                        //let mode = genot::maf_to_mode(maf.unwrap());
                        mode_to_score(mode, (s0, s1, s2))
                    } else if missing_to_mean {
                        match maf {
                            Some(maf) => const_ti + 2.0 * maf * alpha_ti,
                            None=>panic!("MAF is missing but --missing-to-mean is used. Check column name or use --fill-missing-in-dataset to fill missing values in test dataset.")
                        }
                        //const_ti + 2.0 * maf.unwrap() * alpha_ti
                    } else {
                        f64::NAN
                    };

                    calc::scores_add_coef(scores, (s0, s1, s2, sm), &gsnv);
                    //calc::scores_add_coef_nomissing(scores, (s0, s1, s2), &gsnv);
                    // TMP
                    //log::debug!("score_wgt: {:?}", (s0, s1, s2, sm));
                }

                Coef::Score3((s0, s1, s2)) => {
                    let sm = if missing_to_mode {
                        let mode = match maf {
                            Some(maf) => genot::maf_to_mode(maf),
                            None=>panic!("MAF is missing but --missing-to-mode is used. Check column name or use --fill-missing-in-dataset to fill missing values in test dataset.")
                        };

                        // let mode = genot::maf_to_mode(maf.unwrap());
                        mode_to_score(mode, (s0, s1, s2))
                    } else {
                        // missing_to_mean or else
                        f64::NAN
                    };
                    calc::scores_add_coef(scores, (s0, s1, s2, sm), &gsnv);
                    //calc::scores_add_coef_nomissing(scores, (s0, s1, s2), &gsnv);
                }

                Coef::Score4((s0, s1, s2, sm)) => {
                    calc::scores_add_coef(scores, (s0, s1, s2, sm), &gsnv);
                }
                _ => unimplemented!(),
            }
        }
        WgtKind::SnvInteraction(snv_inter_wgt) => {
            let (mi_1, mi_2) = snv_inter_wgt.indexs();
            let mi_1 = match mi_1 {
                None => {
                    if allow_nonexist_snv {
                        return;
                    } else {
                        panic!("SNV does not exist.")
                    }
                }
                Some(mi_1) => mi_1,
            };
            let mi_2 = match mi_2 {
                None => {
                    if allow_nonexist_snv {
                        return;
                    } else {
                        panic!("SNV does not exist.")
                    }
                }
                Some(mi_2) => mi_2,
            };

            let (maf_1, maf_2) = snv_inter_wgt.mafs();
            let (maf_1, maf_2) = (maf_1.unwrap(), maf_2.unwrap());
            if maf_1.is_nan() || maf_2.is_nan() {
                panic!("maf is nan.");
            }

            let genot_mi_1 = genot.unwrap().to_genot_snv(mi_1);
            let genot_mi_2 = genot.unwrap().to_genot_snv(mi_2);

            match wgt.wgt().model().coef() {
                Coef::LinearConstInteraction((const_ti, alpha_ti)) => {
                    if missing_to_mode || missing_to_mean {
                        let score_wgt = interaction_genotype_missing(
                            maf_1,
                            maf_2,
                            missing_to_mode,
                            missing_to_mean,
                        )
                        .linearconst(const_ti, alpha_ti)
                        .to_4x4();
                        // nosimd for missing
                        calc::scores_add_coef_interaction_nosimd(
                            scores,
                            score_wgt,
                            &genot_mi_1,
                            &genot_mi_2,
                        );
                    } else {
                        // assume no missing in genot
                        let score_wgt = interaction_genotype(maf_1, maf_2)
                            .linearconst(const_ti, alpha_ti)
                            .to_tuple();
                        // TMP simd is no missing only
                        calc::scores_add_coef_interaction_no_missing(
                            scores,
                            score_wgt,
                            &genot_mi_1,
                            &genot_mi_2,
                        );
                    }
                }
                _ => unimplemented!(),
            }
        }
        WgtKind::Cov(_) => {
            let cov_name = wgt.wgt().kind().cov().name();

            // TODO: clean const
            if cov_name == "const" {
                match wgt.wgt().model().coef() {
                    Coef::Linear(alpha) => {
                        scores.scores_mut().iter_mut().for_each(|score| {
                            *score += alpha;
                        });
                    }
                    Coef::Binary((const_t, alpha)) => {
                        let score_add = alpha + const_t;
                        scores.scores_mut().iter_mut().for_each(|score| {
                            *score += score_add;
                        });
                    }
                    _ => unimplemented!(),
                }
            } else {
                // cov_name=const will fail
                let cov_vals = covs.unwrap().vals_id(cov_name);
                match wgt.wgt().model().coef() {
                    Coef::Linear(alpha) => {
                        scores
                            .scores_mut()
                            .iter_mut()
                            .zip(cov_vals.iter())
                            .for_each(|(score, val)| {
                                *score += alpha * val;
                            });
                    }
                    _ => unimplemented!(),
                }
            }
        }
    }
}

//pub fn add_score<W: WgtTrait>(
//    scores_pad: &mut [f64],
//    wgt: &W,
//    genot: Option<&Genot>,
//    covs: Option<&Covs>,
//    allow_nonexist_snv: bool,
//) {
//    //log::debug!("wgt {:?}", wgt);
//    match wgt.wgt().kind() {
//        WgtKind::Snv(_, _, mi) => {
//            //log::debug!("mi {}", mi.unwrap());
//
//            let mi = match mi {
//                None => {
//                    if allow_nonexist_snv {
//                        return;
//                    } else {
//                        panic!("SNV does not exist.")
//                    }
//                }
//                Some(mi) => *mi,
//            };
//
//            let gsnv = genot.unwrap().to_genot_snv(mi);
//            //let genot_mi = genot.unwrap().to_genot_snv(mi.unwrap());
//            match wgt.wgt().model().coef() {
//                Coef::Single(alpha) => {
//                    scores_pad.iter_mut().for_each(|score| {
//                        *score += alpha;
//                    });
//                }
//                Coef::Binary((const_ti, alpha_ti)) => {
//                    let threshold = wgt.wgt().model().threshold().unwrap();
//                    let score_add_high = const_ti + alpha_ti;
//                    let score_add_low = const_ti - alpha_ti;
//
//                    // using par_bridge() here is very slow
//                    scores_pad
//                        .iter_mut()
//                        .zip(gsnv.iter())
//                        .for_each(|(score, val)| {
//                            let score_add = if val == 3 {
//                                // for Binary, missing is not allowed.
//                                // (should be `const_ti`??)
//                                panic!("Cannot use mising on Binary coef.");
//                            } else {
//                                match (val as f64) > threshold {
//                                    true => score_add_high,
//                                    false => score_add_low,
//                                }
//                            };
//                            *score += score_add;
//                        });
//                }
//
//                Coef::Linear(alpha_ti) => {
//                    let s0 = 0.0;
//                    let s1 = alpha_ti;
//                    let s2 = 2.0 * alpha_ti;
//
//                    calc::scores_add_coef_nomissing(scores_pad, (s0, s1, s2), &gsnv);
//                    //scores
//                    //    .iter_mut()
//                    //    .zip(genot_mi.iter())
//                    //    .for_each(|(score, val)| {
//                    //        let score_add = match val {
//                    //            0 => s0,
//                    //            1 => s1,
//                    //            2 => s2,
//                    //            //3 => sm,
//                    //            _ => panic!("Sth wrong"),
//                    //        };
//                    //        *score += score_add;
//                    //    });
//                }
//
//                Coef::LinearConst((const_ti, alpha_ti)) => {
//                    let s0 = const_ti;
//                    let s1 = const_ti + alpha_ti;
//                    let s2 = const_ti + 2.0 * alpha_ti;
//
//                    calc::scores_add_coef_nomissing(scores_pad, (s0, s1, s2), &gsnv);
//                    //scores
//                    //    .iter_mut()
//                    //    .zip(genot_mi.iter())
//                    //    .for_each(|(score, val)| {
//                    //        let score_add = match val {
//                    //            0 => s0,
//                    //            1 => s1,
//                    //            2 => s2,
//                    //            //3 => sm,
//                    //            _ => panic!("Sth wrong"),
//                    //        };
//                    //        *score += score_add;
//                    //    });
//                }
//
//                Coef::Score3((s0, s1, s2)) => {
//                    calc::scores_add_coef_nomissing(scores_pad, (s0, s1, s2), &gsnv);
//                    //scores
//                    //    .iter_mut()
//                    //    .zip(genot_mi.iter())
//                    //    .for_each(|(score, val)| {
//                    //        let score_add = match val {
//                    //            0 => s0,
//                    //            1 => s1,
//                    //            2 => s2,
//                    //            //3 => sm,
//                    //            _ => panic!("Sth wrong"),
//                    //        };
//                    //        *score += score_add;
//                    //    });
//                }
//
//                Coef::Score4((s0, s1, s2, sm)) => {
//                    calc::scores_add_coef(scores_pad, (s0, s1, s2, sm), &gsnv);
//                    //scores
//                    //    .iter_mut()
//                    //    .zip(genot_mi.iter())
//                    //    .for_each(|(score, val)| {
//                    //        let score_add = match val {
//                    //            0 => s0,
//                    //            1 => s1,
//                    //            2 => s2,
//                    //            3 => sm,
//                    //            _ => panic!("Sth wrong"),
//                    //        };
//                    //        *score += score_add;
//                    //    });
//                }
//                _ => unimplemented!(),
//            }
//        }
//        WgtKind::SnvInteraction(_, mi_1, _, mi_2) => {
//            let mi_1 = match mi_1 {
//                None => {
//                    if allow_nonexist_snv {
//                        return;
//                    } else {
//                        panic!("SNV does not exist.")
//                    }
//                }
//                Some(mi_1) => *mi_1,
//            };
//            let mi_2 = match mi_2 {
//                None => {
//                    if allow_nonexist_snv {
//                        return;
//                    } else {
//                        panic!("SNV does not exist.")
//                    }
//                }
//                Some(mi_2) => *mi_2,
//            };
//
//            let genot_mi_1 = genot.unwrap().to_genot_snv(mi_1);
//            let genot_mi_2 = genot.unwrap().to_genot_snv(mi_2);
//
//            match wgt.wgt().model().coef() {
//                //Coef::LinearConst((const_ti, alpha_ti)) => {
//                Coef::LinearConstInteraction((const_ti, alpha_ti)) => {
//                    let s0 = const_ti;
//                    let s1 = const_ti + alpha_ti;
//                    let s2 = const_ti + 2.0 * alpha_ti;
//                    let s4 = const_ti + 4.0 * alpha_ti;
//
//                    calc::scores_add_coef_interaction(
//                        scores_pad,
//                        (s0, s1, s2, s4),
//                        &genot_mi_1,
//                        &genot_mi_2,
//                    );
//                }
//                _ => unimplemented!(),
//            }
//        }
//        WgtKind::Cov(_) => {
//            let cov_name = wgt.wgt().kind().cov().name();
//
//            // TODO: clean const
//            if cov_name == "const" {
//                match wgt.wgt().model().coef() {
//                    Coef::Linear(alpha) => {
//                        scores_pad.iter_mut().for_each(|score| {
//                            *score += alpha;
//                        });
//                    }
//                    Coef::Binary((const_t, alpha)) => {
//                        let score_add = alpha + const_t;
//                        scores_pad.iter_mut().for_each(|score| {
//                            *score += score_add;
//                        });
//                    }
//                    _ => unimplemented!(),
//                }
//            } else {
//                // cov_name=const will fail
//                let cov_vals = covs.unwrap().vals_id(cov_name);
//                //let cov_vals = covs.vals_id(cov_name);
//                //let cov_vals = cov_wgt.var().vals();
//                match wgt.wgt().model().coef() {
//                    Coef::Linear(alpha) => {
//                        scores_pad
//                            .iter_mut()
//                            .zip(cov_vals.iter())
//                            .for_each(|(score, val)| {
//                                *score += alpha * val;
//                            });
//                        //scores
//                        //    .iter_mut()
//                        //    .zip(cov_vals.iter())
//                        //    .par_bridge()
//                        //    .for_each(|(score, val)| {
//                        //        *score += alpha * val;
//                        //    });
//
//                        //for (score, val) in scores.iter_mut().zip(cov_vals.iter()) {
//                        //    *score += alpha * val;
//                        //}
//                    }
//                    _ => unimplemented!(),
//                }
//            }
//        }
//    }
//}

/*
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

// TODO:
// DO NOT use rayon here
// since very slow since n is small
// TODO: test
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
 */

//#[allow(dead_code)]
//pub fn score(fout_score: &Path, wgts: &Wgts, dataset: &Dataset, samples_id: &[String]) {
//    let n = dataset.genot().n();
//    let genot = dataset.genot();
//    let phe = dataset.samples().phe_unwrap();
//
//    let mut scores = vec![0.0f64; n];
//    for (wgt_i, wgt) in wgts.wgts().iter().enumerate() {
//        if (wgt_i > 0) && (wgt_i % 100_000 == 0) {
//            log::debug!("{:?}-th SNV", wgt_i);
//        }
//
//        add_score(&mut scores, wgt, genot, None);
//    }
//
//    write_scores(&fout_score, &scores, phe, samples_id);
//    //io::write_scores(&fout_iteration, &scores, phe, samples_id);
//}

/// Use rayon
pub fn calc_score_wgts(
    wgts: &Wgts,
    dataset: &Dataset,
    allow_nonexist_snv: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
) -> SampleScore {
    let n = dataset.genot().n();
    let genot = dataset.genot();

    let scores_sum = wgts
        .wgts()
        .par_iter()
        .map(|wgt| {
            calc_score_wgt(
                n,
                wgt,
                Some(genot),
                None,
                allow_nonexist_snv,
                missing_to_mode,
                missing_to_mean,
            )
        })
        .reduce(
            || SampleScore::new(n),
            |mut a, b| {
                a.add_scores(&b);
                a
                //a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                //    *x += y;
                //});
                //a
            },
        );

    //.map(|wgt| calc_score_wgt(n, wgt, Some(genot), None, true))

    scores_sum
}

//pub fn calc_score_wgts(wgts: &Wgts, dataset: &Dataset) -> Vec<f64> {
//    let n = dataset.genot().n();
//    let genot = dataset.genot();
//
//    let scores_sum = wgts
//        .wgts()
//        .par_iter()
//        .map(|wgt| calc_score_wgt(n, wgt, Some(genot), None, true))
//        .reduce(
//            || vec![0.0f64; n],
//            |mut a, b| {
//                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
//                    *x += y;
//                });
//                a
//            },
//        );
//
//    //let mut scores = vec![0.0f64; n];
//    //for (wgt_i, wgt) in wgts.wgts().iter().enumerate() {
//    //    if (wgt_i > 0) && (wgt_i % 100_000 == 0) {
//    //        log::debug!("{:?}-th SNV", wgt_i);
//    //    }
//
//    //    add_score(&mut scores, wgt, Some(genot), None, true);
//    //}
//
//    //write_scores(&fout_score, &scores, phe, samples_id);
//    scores_sum
//}

// no rayon
//pub fn calculate_score(wgts: &Wgts, dataset: &Dataset) -> Vec<f64> {
//    let n = dataset.genot().n();
//    let genot = dataset.genot();
//
//    let mut scores = vec![0.0f64; n];
//    for (wgt_i, wgt) in wgts.wgts().iter().enumerate() {
//        if (wgt_i > 0) && (wgt_i % 100_000 == 0) {
//            log::debug!("{:?}-th SNV", wgt_i);
//        }
//
//        add_score(&mut scores, wgt, Some(genot), None, true);
//        //add_score(&mut scores, wgt, genot, None);
//    }
//
//    //write_scores(&fout_score, &scores, phe, samples_id);
//    scores
//}

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

////fn write_scores(fout: &Path, scores: &[f64], phe: &Phe, samples_id: &[(String, String)]) {
//#[allow(dead_code)]
//pub fn write_scores(fout: &Path, scores: &[f64], phe: &Phe, samples_id: &[String]) {
//    let file = match File::create(&fout) {
//        Ok(file) => file,
//        Err(_) => panic!(
//            "Cannot create file, possibly directory does not exist: {:?}",
//            &fout
//        ),
//    };
//
//    let mut writer = BufWriter::new(file);
//    // assume word count of one line is 30
//    // no problem when it exceeds
//    let capacity = 30 * scores.len();
//    let mut str = String::with_capacity(capacity);
//    str.push_str("fid\tiid\tphe\trs\n");
//
//    for ni in 0..scores.len() {
//        str.push_str(&samples_id[ni]);
//        str.push_str("\t");
//        str.push_str(&samples_id[ni]);
//        str.push_str("\t");
//        str.push_str(&(phe.get_unchecked(ni) as u8).to_string());
//        str.push_str("\t");
//        str.push_str(&format!("{:.5}\n", scores[ni]));
//    }
//
//    writer.write(str.as_bytes()).unwrap();
//
//    //for ni in
//    // use .concat(), .join()?
//    // https://users.rust-lang.org/t/fast-string-concatenation/4425/5
//    // -> .push_str seems fastest
//    // -> could be because use with_capacity beforehand
//}

pub fn write_scores(fout: &Path, scores: &SampleScore, samples_id: &[String]) {
    write_scores_vec(fout, scores.scores(), samples_id);
}

pub fn write_scores_vec(fout: &Path, scores: &[f64], samples_id: &[String]) {
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

    //for ni in 0..scores.len() {
    for ni in 0..samples_id.len() {
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
    // -> probabily due to using with_capacity beforehand
}

// use vec::has_unique_elements
//fn is_unique<T>(iter: T) -> bool
//where
//    T: IntoIterator,
//    T::Item: Eq + Hash,
//{
//    let mut uniq = HashSet::new();
//    iter.into_iter().all(move |x| uniq.insert(x))
//}

pub fn write_scores_paras(
    fout: &Path,
    score_paras: &[SampleScore],
    concat_para: &str,
    paras: &[String],
    samples_id: &[String],
) {
    let scores_paras_vec = score_paras.iter().map(|x| x.scores()).collect::<Vec<_>>();
    write_scores_paras_vec(fout, &scores_paras_vec, concat_para, paras, samples_id);
}

pub fn write_scores_paras_vec(
    fout: &Path,
    score_paras: &[&[f64]],
    //score_paras: &[Vec<f64>],
    concat_para: &str,
    paras: &[String],
    samples_id: &[String],
) {
    if score_paras.len() != paras.len() {
        panic!("score_paras.len() != paras.len()");
    }

    //if !is_unique(paras) {
    if !vec::has_unique_elements(paras) {
        panic!("paras should be unique.");
    }

    let file = match File::create(&fout) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {:?}",
            &fout
        ),
    };

    let sample_n = samples_id.len();
    //let sample_n = score_paras[0].len();

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
