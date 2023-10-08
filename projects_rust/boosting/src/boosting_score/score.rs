// TODO: DO NOT USE wgti as iter!!! -> these two should be separated!! should refer to wgt.iter()

use std::collections::HashSet;
use std::path::Path;

use crate::wgt_boosts::WgtBoosts;
use crate::WgtBoost;
use genetics::genot::prelude::*;
use genetics::samples::CovsTrait;
use genetics::wgt::{Coef, WgtKind};
use genetics::{vec, Covs, Dataset};

//use rayon::iter::ParallelBridge;
//use rayon::prelude::*;
//use rayon::iter::ParallelBridge;

// TODO: use SIMD
//
// DO NOT use rayon here
// since very slow when n is small
// though nested rayon does not raise error
//
// should implement sth in Wgt, Model to make this simple
// split for cov?
// TODO: integrate to genetics::score::add_score()
pub fn add_score(scores: &mut [f64], wgt: &WgtBoost, genot: &Genot, covs: Option<&Covs>) {
    //log::debug!("wgt {:?}", wgt);
    match wgt.wgt().kind() {
        WgtKind::Snv(_, _, mi) => {
            //log::debug!("mi {}", mi.unwrap());
            // This should never panic. => ?
            let genot_mi = genot.to_genot_snv(mi.unwrap());
            match wgt.wgt().model().coef() {
                Coef::Binary((const_ti, alpha_ti)) => {
                    let threshold = wgt.wgt().model().threshold().unwrap();
                    let score_add_high = const_ti + alpha_ti;
                    let score_add_low = const_ti - alpha_ti;

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

                    // use par_bridge(): very slow
                    //scores
                    //    .iter_mut()
                    //    .zip(genot_mi.iter())
                    //    .par_bridge()
                    //    .for_each(|(score, val)| {
                    //        let score_add = if val == 3 {
                    //            // FIXME: is this correct???
                    //            // should be `const_ti`??
                    //            // or NEVER HAPPEN?
                    //            panic!("This should never happen??.");
                    //            0.0
                    //        } else {
                    //            match (val as f64) > threshold {
                    //                true => score_add_high,
                    //                false => score_add_low,
                    //            }
                    //        };
                    //        *score += score_add;
                    //    });

                    //for (score, val) in scores.iter_mut().zip(genot_mi.iter()) {
                    //    let score_add = if val == 3 {
                    //        0.0
                    //    } else {
                    //        match (val as f64) > threshold {
                    //            true => score_add_high,
                    //            false => score_add_low,
                    //        }
                    //    };
                    //    *score += score_add;
                    //}
                }
                Coef::LinearConst((const_ti, alpha_ti)) => {
                    let s2 = const_ti + 2.0 * alpha_ti;
                    let s1 = const_ti + alpha_ti;
                    let s0 = const_ti;

                    scores
                        .iter_mut()
                        .zip(genot_mi.iter())
                        .for_each(|(score, val)| {
                            let score_add = match val {
                                2 => s2,
                                1 => s1,
                                0 => s0,
                                //3 => sm,
                                _ => panic!(""),
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
                                _ => panic!(""),
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
                                _ => panic!(""),
                            };
                            *score += score_add;
                        });

                    //scores
                    //    .iter_mut()
                    //    .zip(genot_mi.iter())
                    //    .par_bridge()
                    //    .for_each(|(score, val)| {
                    //        let score_add = match val {
                    //            2 => s2,
                    //            1 => s1,
                    //            0 => s0,
                    //            3 => sm,
                    //            _ => panic!(""),
                    //        };
                    //        *score += score_add;
                    //    });

                    //for (score, val) in scores.iter_mut().zip(genot_mi.iter()) {
                    //    let score_add = match val {
                    //        2 => s2,
                    //        1 => s1,
                    //        0 => s0,
                    //        3 => sm,
                    //        _ => panic!(""),
                    //    };
                    //    *score += score_add;
                    //}
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
                        //scores.par_iter_mut().for_each(|score| {
                        //    *score += alpha;
                        //});

                        //for score in scores.iter_mut() {
                        //    *score += alpha;
                        //}
                    }
                    Coef::Binary((const_t, alpha)) => {
                        let score_add = alpha + const_t;
                        scores.iter_mut().for_each(|score| {
                            *score += score_add;
                        });

                        //scores.par_iter_mut().for_each(|score| {
                        //    *score += score_add;
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

/// a item refer to a wgt, which is different from iter
fn calculate_write_score(
    dout: &Path,
    item_last_indexs_write: &[usize], // this should indicate the last index, not number of items
    wgts: &WgtBoosts,
    dataset: &Dataset,
    //samples_id: &[String],
    //samples_id: &[(String, String)],
    nocov: bool,
    is_nsnv: bool,                   // for fout name
    item_ns_fname: Option<&[usize]>, // when fname = items_write
    integrate: bool,
) {
    let item_ns_fname = match item_ns_fname {
        Some(v) => v,
        None => item_last_indexs_write,
    };

    let n = dataset.genot().n();
    let genot = dataset.genot();
    //let phe = dataset.samples().phe_unwrap();
    let covs = dataset.samples().covs();
    //let covs = dataset.samples().covs().unwrap();
    let samples_id = dataset.samples().names();

    let mut score_paras: Vec<Vec<f64>> = vec![];
    let mut paras: Vec<String> = vec![];

    let mut scores = vec![0.0f64; n];

    let mut items_write_pair = item_last_indexs_write.iter().zip(item_ns_fname.iter());
    //let mut iters_write = iters_write.clone().iter();
    let (mut item_next, mut item_fname) = items_write_pair.next().unwrap();
    //let (mut iter_i, mut iter_next) = iters_write.enumerate().next().unwrap();
    //let mut iter_next = *iters_write.next().unwrap();

    for (item_i, wgt) in wgts.wgts().iter().enumerate() {
        //log::debug!("wgt {:?}", wgt);

        if !(nocov && wgt.wgt().is_cov()) {
            add_score(&mut scores, wgt, genot, covs);
        }
        //log::debug!("iter_next,use{},{}", iter_next, iter_use);
        if item_i == *item_next {
            //write
            //let fout_iteration =
            //super::io::fname_score_createdir(dout, *item_fname, nocov, is_nsnv, integrate);
            log::debug!("Save iteration: {}", item_fname);
            paras.push(item_fname.to_string());

            score_paras.push(scores.clone());
            //super::io::write_scores(&fout_iteration, &scores, samples_id);

            // raise error...
            //(iter_next, iter_fname) = match iters_write_pair.next() {
            let v = match items_write_pair.next() {
                Some(v) => v,
                None => break,
            };
            // raise error...
            //(iter_next, iter_fname) = v;
            item_next = v.0;
            item_fname = v.1;
        }
    }

    // write
    let fout_concat = super::io::fname_score_concat_createdir(dout, nocov, is_nsnv, integrate);
    log::debug!("Write score: {:?}", fout_concat);
    super::io::write_scores_concat(
        &fout_concat,
        &score_paras,
        &paras,
        samples_id,
        is_nsnv,
        integrate,
    );
}

// /// a item refer to a wgt, which is different from iter
// fn calculate_write_score(
//     dout: &Path,
//     item_last_indexs_write: &[usize], // this should indicate the last index, not number of items
//     wgts: &WgtBoosts,
//     dataset: &Dataset,
//     //samples_id: &[String],
//     //samples_id: &[(String, String)],
//     nocov: bool,
//     is_nsnv: bool,                   // for fout name
//     item_ns_fname: Option<&[usize]>, // when fname = items_write
//     integrate: bool,
// ) {
//     let item_ns_fname = match item_ns_fname {
//         Some(v) => v,
//         None => item_last_indexs_write,
//     };

//     let n = dataset.genot().n();
//     let genot = dataset.genot();
//     //let phe = dataset.samples().phe_unwrap();
//     let covs = dataset.samples().covs();
//     //let covs = dataset.samples().covs().unwrap();
//     let samples_id = dataset.samples().names();

//     let mut scores = vec![0.0f64; n];

//     let mut items_write_pair = item_last_indexs_write.iter().zip(item_ns_fname.iter());
//     //let mut iters_write = iters_write.clone().iter();
//     let (mut item_next, mut item_fname) = items_write_pair.next().unwrap();
//     //let (mut iter_i, mut iter_next) = iters_write.enumerate().next().unwrap();
//     //let mut iter_next = *iters_write.next().unwrap();

//     for (item_i, wgt) in wgts.wgts().iter().enumerate() {
//         //log::debug!("wgt {:?}", wgt);

//         if !(nocov && wgt.wgt().is_cov()) {
//             add_score(&mut scores, wgt, genot, covs);
//             //add_score(&mut scores, wgt, genot, Some(covs));
//         }
//         //log::debug!("iter_next,use{},{}", iter_next, iter_use);
//         if item_i == *item_next {
//             //log::debug!("Write iter {}", iter_next);
//             //write
//             let fout_iteration =
//                 super::io::fname_score_createdir(dout, *item_fname, nocov, is_nsnv, integrate);
//             log::debug!("Write iteration {}: {:?}", item_fname, fout_iteration);

//             super::io::write_scores(&fout_iteration, &scores, samples_id);
//             //super::io::write_scores(&fout_iteration, &scores, phe, samples_id);

//             //log::debug!("Done write iteration {}: {:?}", item_fname, fout_iteration);

//             // raise error...
//             //(iter_next, iter_fname) = match iters_write_pair.next() {
//             let v = match items_write_pair.next() {
//                 Some(v) => v,
//                 None => break,
//             };
//             // raise error...
//             //(iter_next, iter_fname) = v;
//             item_next = v.0;
//             item_fname = v.1;
//         }
//     }
// }

/// Iteration and item number are different.
/// return corresponding items index to given iterations.
/// assume iter_until_item is monotonically increased
/// for nsnvs, duplicated snvs are not newly counted.
/// ex. when iters_write=[1,3,5] and iter_until_item=[1,2,2,2,3,3,4]
/// then return [0,5]
/// since write=1 -> item index=0,  write=3 -> item index=5 and write=5 is out of the bound
fn create_item_last_indexs_write(iters_write: &[usize], iter_until_item: &[usize]) -> Vec<usize> {
    let iter_max = *iter_until_item.iter().last().unwrap();
    // find corresponding iter from backward
    // exclude nsnv > max
    let mut nsnvs_write_rev = iters_write.iter().filter(|&v| *v <= iter_max).rev();

    let mut nsnv_next = *nsnvs_write_rev.next().unwrap();
    let mut item_last_indexs_write_back: Vec<usize> = Vec::new();
    for (item_i, nsnv_until_iter) in iter_until_item.iter().enumerate().rev() {
        //log::debug!("item_i, nsnv_until_iter {},{}", item_i, nsnv_until_iter);
        //log::debug!("nsnv_next {}", nsnv_next);
        if *nsnv_until_iter == nsnv_next {
            item_last_indexs_write_back.push(item_i);
            nsnv_next = match nsnvs_write_rev.next() {
                Some(v) => *v,
                None => break,
            };
        }
    }

    let item_last_indexs_write = item_last_indexs_write_back
        .iter()
        .rev()
        .map(|v| *v)
        .collect::<Vec<usize>>();

    item_last_indexs_write
}

pub fn calculate_write_score_iterations(
    dout: &Path,
    iterations_write: &[usize],
    wgts: &WgtBoosts,
    dataset: &Dataset,
    //samples_id: &[String],
    //samples_id: &[(String, String)],
    nocov: bool,
) {
    // iteration index -> number of iterations
    // so iteration+1
    let iter_until_item = wgts
        .wgts()
        .iter()
        .map(|v| v.iteration() + 1)
        .collect::<Vec<usize>>();
    assert!(vec::is_sorted(&iter_until_item));

    //log::debug!("iter_until_item {:?}", iter_until_item);
    //log::debug!("iterations_write {:?}", iterations_write);

    let item_last_indexs_write = create_item_last_indexs_write(iterations_write, &iter_until_item);

    log::debug!(
        "iters to write corresponding to nsnvs are {:?}",
        item_last_indexs_write
    );

    if item_last_indexs_write.len() == 0 {
        log::debug!("No iterations to be written.");
        return;
    }

    calculate_write_score(
        dout,
        &item_last_indexs_write,
        wgts,
        dataset,
        //samples_id,
        nocov,
        false,
        Some(iterations_write),
        false,
    );
}

pub fn calculate_write_score_nsnvs(
    dout: &Path,
    nsnvs_write: &[usize], // or snv_ns_write
    wgts: &WgtBoosts,
    dataset: &Dataset,
    //samples_id: &[(String, String)],
    //samples_id: &[String],
    nocov: bool,
) {
    // first, count when (at which iteration) to write score
    let mut nsnvs_until_item: Vec<usize> = Vec::with_capacity(wgts.wgts().len());
    // count non-duplicated used snvs
    let mut snv_used: HashSet<String> = HashSet::with_capacity(*nsnvs_write.last().unwrap());
    let mut count_unique_snv = 0;
    for wgt in wgts.wgts().iter() {
        if wgt.wgt().is_snv() {
            // use rs or sida?? -> sida should be fine since assured to be unique
            let snv_name = wgt.wgt().kind().snv_index().sida();
            if !snv_used.contains(snv_name) {
                count_unique_snv += 1;
                snv_used.insert(snv_name.to_owned());
            }
        }
        nsnvs_until_item.push(count_unique_snv);
    }

    //log::debug!("nsnvs_until_item {:?}", nsnvs_until_item);

    assert!(vec::is_sorted(&nsnvs_until_item));

    let item_last_indexs_write = create_item_last_indexs_write(nsnvs_write, &nsnvs_until_item);

    log::debug!(
        "iters to write corresponding to nsnvs are {:?}",
        item_last_indexs_write
    );

    if item_last_indexs_write.len() == 0 {
        log::debug!("No iterations to be written.");
        return;
    }

    calculate_write_score(
        dout,
        &item_last_indexs_write,
        wgts,
        dataset,
        //samples_id,
        nocov,
        true,
        Some(nsnvs_write),
        false,
    );
}

pub fn calculate_write_score_para_best(
    dout: &Path,
    wgts: &WgtBoosts,
    dataset: &Dataset,
    //samples_id: &[String],
    nocov: bool,
) {
    let iterations_write = [wgts.wgts().len()];

    // iteration index -> number of iterations
    // so iteration+1
    let iter_until_item = wgts
        .wgts()
        .iter()
        .map(|v| v.iteration() + 1)
        .collect::<Vec<usize>>();
    assert!(vec::is_sorted(&iter_until_item));

    //log::debug!("iter_until_item {:?}", iter_until_item);
    //log::debug!("iterations_write {:?}", iterations_write);

    let item_last_indexs_write = create_item_last_indexs_write(&iterations_write, &iter_until_item);

    log::debug!(
        "iters to write corresponding to nsnvs are {:?}",
        item_last_indexs_write
    );

    if item_last_indexs_write.len() == 0 {
        log::debug!("No iterations to be written.");
        return;
    }

    calculate_write_score(
        dout,
        &item_last_indexs_write,
        wgts,
        dataset,
        //samples_id,
        nocov,
        false,
        Some(&iterations_write),
        true,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_item_last_indexs_write() {
        let iters_write = [1, 2, 100];
        let iter_until_item = [0, 1, 2, 2];
        let item_last_indexs_write = create_item_last_indexs_write(&iters_write, &iter_until_item);
        assert_eq!(item_last_indexs_write, vec![1, 3]);
    }
}
