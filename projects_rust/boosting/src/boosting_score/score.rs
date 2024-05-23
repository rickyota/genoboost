use std::collections::HashSet;
use std::path::Path;

use crate::dout_file::DoutScoreParaFile;
use crate::wgt_boosts::WgtBoosts;
use genetics::genot::prelude::*;
use genetics::score as gscore;
use genetics::SampleScore;
use genetics::{vec, Dataset};
use genetics::{CovsTrait, WgtKind};

/// a item refer to a wgt, which is different from iter
fn calculate_write_score(
    dout: &DoutScoreParaFile,
    item_last_indexs_write: &[usize], // this should indicate the last index, not number of items
    wgts: &WgtBoosts,
    dataset: &Dataset,
    nocov: bool,
    is_nsnv: bool,                   // for fout name
    item_ns_fname: Option<&[usize]>, // when fname = items_write
    integrate: bool,
    allow_nonexist_snv: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
) {
    let item_ns_fname = match item_ns_fname {
        Some(v) => v,
        None => item_last_indexs_write,
    };

    let n = dataset.genot().n();
    let genot = dataset.genot();
    let covs = dataset.samples().covs();
    let samples_id = dataset.samples().names();

    let mut score_paras: Vec<SampleScore> = vec![];
    let mut paras: Vec<String> = vec![];

    let mut scores = SampleScore::new(n);

    let mut items_write_pair = item_last_indexs_write.iter().zip(item_ns_fname.iter());
    let (mut item_next, mut item_fname) = items_write_pair.next().unwrap();

    for (item_i, wgt) in wgts.wgts().iter().enumerate() {
        //log::debug!("wgt {:?}", wgt);

        if !(nocov && wgt.wgt().is_cov()) {
            // missing_to_mode=true in boosting
            gscore::add_score(
                &mut scores,
                wgt,
                Some(genot),
                covs,
                allow_nonexist_snv,
                missing_to_mode, //true,
                missing_to_mean, //false,
            );

            scores.check_no_nan();
        }
        if item_i == *item_next {
            log::debug!("Save iteration: {}", item_fname);
            paras.push(item_fname.to_string());

            score_paras.push(scores.clone_align());

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
    let fout_concat = dout.fname_score_concat_createdir(nocov, is_nsnv, integrate);
    log::debug!("Write score: {:?}", fout_concat);
    write_scores_concat(
        &fout_concat,
        &score_paras,
        &paras,
        samples_id,
        is_nsnv,
        integrate,
    );
}

pub fn write_scores_concat(
    fout: &Path,
    score_paras: &[SampleScore],
    paras: &[String],
    samples_id: &[String],
    is_nsnv: bool,
    integrate: bool,
) {
    if integrate {
        assert_eq!(score_paras.len(), 1);
        gscore::write_scores(fout, &score_paras[0], samples_id)
    } else {
        if is_nsnv {
            gscore::write_scores_paras(fout, score_paras, "n", paras, samples_id)
        } else {
            gscore::write_scores_paras(fout, score_paras, "iter", paras, samples_id)
        }
    }
}

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
    dout: &DoutScoreParaFile,
    iterations_write: &[usize],
    wgts: &WgtBoosts,
    dataset: &Dataset,
    nocov: bool,
    allow_nonexist_snv: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
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
        nocov,
        false,
        Some(iterations_write),
        false,
        allow_nonexist_snv,
        missing_to_mode, //true,
        missing_to_mean, //false,
    );
}

pub fn calculate_write_score_nsnvs(
    dout: &DoutScoreParaFile,
    nsnvs_write: &[usize], // or snv_ns_write
    wgts: &WgtBoosts,
    dataset: &Dataset,
    nocov: bool,
    allow_nonexist_snv: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
) {
    // first, count when (at which iteration) to write score
    let mut nsnvs_until_item: Vec<usize> = Vec::with_capacity(wgts.wgts().len());
    // count non-duplicated used snvs
    let mut snv_used: HashSet<String> = HashSet::with_capacity(*nsnvs_write.last().unwrap());
    let mut count_unique_snv = 0;
    for wgt in wgts.wgts().iter() {
        match wgt.wgt().kind() {
            WgtKind::Snv(snv_wgt) => {
                let snv_id = snv_wgt.snv_id();
                // here, snv name can be anything but it should be unique.
                // use rs or sida?? -> sida should be fine since assured to be unique
                let snv_name = snv_id.sida();
                if !snv_used.contains(snv_name) {
                    count_unique_snv += 1;
                    snv_used.insert(snv_name.to_owned());
                }
            }
            WgtKind::SnvInteraction(snv_inter_wgt) => {
                let (snv_id_1, snv_id_2) = snv_inter_wgt.snv_ids();
                let snv_name_1 = snv_id_1.sida();
                if !snv_used.contains(snv_name_1) {
                    count_unique_snv += 1;
                    snv_used.insert(snv_name_1.to_owned());
                }

                let snv_name_2 = snv_id_2.sida();
                if !snv_used.contains(snv_name_2) {
                    count_unique_snv += 1;
                    snv_used.insert(snv_name_2.to_owned());
                }
            }
            WgtKind::Cov(..) => {}
        }
        //if wgt.wgt().is_snv() {
        //    // use rs or sida?? -> sida should be fine since assured to be unique
        //    let snv_name = wgt.wgt().kind().snv_index().sida();
        //    if !snv_used.contains(snv_name) {
        //        count_unique_snv += 1;
        //        snv_used.insert(snv_name.to_owned());
        //    }
        //}

        // TOFIX: it snv interation has two newly appeared snvs, some monitoring snv number could be skipped.
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
        nocov,
        true,
        Some(nsnvs_write),
        false,
        allow_nonexist_snv,
        missing_to_mode, //true,
        missing_to_mean, //false,
    );
}

pub fn calculate_write_score_para_best(
    dout: &DoutScoreParaFile,
    wgts: &WgtBoosts,
    dataset: &Dataset,
    nocov: bool,
    allow_nonexist_snv: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
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
        nocov,
        false,
        Some(&iterations_write),
        true,
        allow_nonexist_snv,
        missing_to_mode, //true,
        missing_to_mean, //false,
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
