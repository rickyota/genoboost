use crate::WgtBoost;
use genetics::genot::prelude::*;
//use genetics::wgt::{Coef, WgtKind};
use genetics::Covs;
use genetics::SampleScore;

use genetics::score as gscore;

// DO NOT deprecate
// This can avoid bug of allow_nonexist_snv
pub fn add_score_training(
    scores: &mut SampleScore,
    wgt: &WgtBoost,
    genot: Option<&Genot>,
    covs: Option<&Covs>,
) {
    // missing is filled in training
    gscore::add_score(scores, wgt, genot, covs, false, false, false);

    //  no check for speed
    //scores.check_no_nan();

}
//pub fn add_score(scores_pad: &mut [f64], wgt: &WgtBoost, genot: Option<&Genot>, covs: Option<&Covs>) {
//    gscore::add_score(scores_pad, wgt, genot, covs, false);
//}

/*
// moved to genetics::score::add_score()
// TODO: use SIMD
// TODO: test
// TODO: merge to genetics::score()
//
// DO NOT use rayon here
// since very slow when n is small
// though nested rayon does not raise error
//
// should implement sth in Wgt, Model to make this simple
pub fn add_score(scores: &mut [f64], wgt: &WgtBoost, genot: Option<&Genot>, covs: Option<&Covs>) {
    //log::debug!("wgt {:?}", wgt);
    match wgt.wgt().kind() {
        WgtKind::Snv(_, _, mi) => {
            //log::debug!("mi {}", mi.unwrap());
            // This should never panic. => ?
            let genot_mi = genot.unwrap().to_genot_snv(mi.unwrap());
            match wgt.wgt().model().coef() {
                Coef::Single(alpha) => {
                    scores.iter_mut().for_each(|score| {
                        *score += alpha;
                    });
                }
                Coef::Binary((const_ti, alpha_ti)) => {
                    let threshold = wgt.wgt().model().threshold().unwrap();
                    let score_add_high = const_ti + alpha_ti;
                    let score_add_low = const_ti - alpha_ti;

                    // using par_bridge() here is very slow
                    scores
                        .iter_mut()
                        .zip(genot_mi.iter())
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

                Coef::LinearConst((const_ti, alpha_ti)) => {
                    let s2 = const_ti + 2.0 * alpha_ti;
                    let s1 = const_ti + alpha_ti;
                    let s0 = const_ti;

                    // TODO: extract fn scores_add_coef_nomissing()
                    scores
                        .iter_mut()
                        .zip(genot_mi.iter())
                        .for_each(|(score, val)| {
                            let score_add = match val {
                                2 => s2,
                                1 => s1,
                                0 => s0,
                                //3 => sm,
                                _ => panic!("Sth wrong"),
                            };
                            *score += score_add;
                        });
                }

                Coef::Score3((s0, s1, s2)) => {
                    scores
                        .iter_mut()
                        .zip(genot_mi.iter())
                        .for_each(|(score, val)| {
                            let score_add = match val {
                                2 => s2,
                                1 => s1,
                                0 => s0,
                                //3 => sm,
                                _ => panic!("Sth wrong"),
                            };
                            *score += score_add;
                        });
                }

                Coef::Score4((s0, s1, s2, sm)) => {
                    scores
                        .iter_mut()
                        .zip(genot_mi.iter())
                        .for_each(|(score, val)| {
                            let score_add = match val {
                                2 => s2,
                                1 => s1,
                                0 => s0,
                                3 => sm,
                                _ => panic!("Sth wrong"),
                            };
                            *score += score_add;
                        });
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
                        scores.iter_mut().for_each(|score| {
                            *score += alpha;
                        });
                    }
                    Coef::Binary((const_t, alpha)) => {
                        let score_add = alpha + const_t;
                        scores.iter_mut().for_each(|score| {
                            *score += score_add;
                        });
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
