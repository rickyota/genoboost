//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
//! Loading genotype requires the long time. Extracting genotype takes negligible time.
//!
// split covs from samples

//mod correlation; // not implemented cmatrix yet
mod coefficient;
// mod compute_pred;
mod epsilon;
mod iteration;
pub mod loss; // pub for bench
mod regression_cov;
pub mod sample_weight;
mod table;
mod update_wgt;

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::iter::FromIterator;
use std::path::Path;
use std::time::Instant;

use self::sample_weight::SampleWeight;
use super::WgtBoost;
use crate::boosting_param::{
    AccMetric, BatchInteractionWay, BatchWay, BoostParam, BoostParamCommonTrait, BoostType, CovWay,
    InteractionWay, IterationNumber, LossFunc,
};
use crate::score;
use crate::wgt_boost;
use crate::wgt_boosts::WgtBoosts;
use crate::DoutParaFile;
use genetics::genot::prelude::*;
use genetics::samples::prelude::*;
use genetics::{pgs, Chrom};
use genetics::{vec, LdCriteria};
use genetics::{CovsTrait, Dataset, Snvs};
use genetics::{Model, SampleScore};
use genetics::{Wgt, WgtKind, WgtTrait};
use loss::LossStruct;
pub use table::ContingencyTable;

fn boosting_iter_snv(
    iteration: usize,
    boost_param: &BoostParam,
    loss: &mut LossStruct,
    dataset: &Dataset,
    scores: &mut SampleScore,
    sample_weight: &mut SampleWeight,
    extract_snvs: Option<&HashSet<usize>>,
    dataset_val: Option<&Dataset>,
    scores_val: &mut SampleScore,
) -> WgtBoost {
    let genot = dataset.genot();
    let phe = dataset.samples().phe_unwrap();

    sample_weight.print_stat(phe);

    loss::calculate_loss_gt(
        loss,
        dataset,
        sample_weight,
        boost_param,
        extract_snvs,
        None,
    );

    let mut wgt = loss::search_min_loss_gt(loss, iteration, dataset, boost_param, extract_snvs);

    // set coef
    update_wgt::update_coef(&mut wgt, boost_param, phe, genot, sample_weight, dataset);

    score::add_score_training(scores, &wgt, Some(genot), None);

    if let Some(dataset_val) = dataset_val {
        let genot_val = dataset_val.genot();
        score::add_score_training(scores_val, &wgt, Some(genot_val), None)
    }

    sample_weight.renew_sample_weight(scores, phe);

    // if phe.n() < 10 {
    //     log::debug!("score {:?}", scores);
    //     log::debug!("ws {:?}", ws);
    //     log::debug!("ps {:?}", ps);
    // }

    wgt
}

fn boosting_iter_snv_interaction(
    iteration: usize,
    boost_param: &BoostParam,
    loss: &mut LossStruct,
    dataset: &Dataset,
    scores: &mut SampleScore,
    sample_weight: &mut SampleWeight,
    // for additive
    extract_snvs: Option<&HashSet<usize>>,
    // for interaction
    extract_interaction: &Vec<(usize, usize)>,
    dataset_val: Option<&Dataset>,
    scores_val: &mut SampleScore,
) -> WgtBoost {
    let genot = dataset.genot();
    let phe = dataset.samples().phe_unwrap();

    sample_weight.print_stat(phe);

    loss::calculate_loss_gt_interaction(
        loss,
        dataset,
        sample_weight,
        boost_param,
        extract_snvs,
        extract_interaction,
    );

    let mut wgt = loss::search_min_loss_gt_interaction(
        loss,
        iteration,
        dataset,
        boost_param,
        extract_snvs,
        extract_interaction,
    );

    // set coef
    update_wgt::update_coef(&mut wgt, boost_param, phe, genot, sample_weight, dataset);

    score::add_score_training(scores, &wgt, Some(genot), None);

    if let Some(dataset_val) = dataset_val {
        let genot_val = dataset_val.genot();
        score::add_score_training(scores_val, &wgt, Some(genot_val), None)
    }

    sample_weight.renew_sample_weight(scores, phe);

    // if phe.n() < 10 {
    //     log::debug!("score {:?}", scores);
    //     log::debug!("ws {:?}", ws);
    //     log::debug!("ps {:?}", ps);
    // }

    wgt
}

/*
// TODO: extract common with boosting_iter_snv()
fn boosting_iter_snv_second(
    ti: usize,
    boost_param: BoostParam,
    loss: &mut LossStruct,
    genot: &Genot,
    pred: &mut [u8],
    scores: &mut [f64],
    ws: &mut [f64],
    ps: &mut [f64],
    phe: &Phe,
    snvs: &Snvs,
    wgt_first: &WgtBoost,
) -> WgtBoost {
    // in is_dom_rec

    let mut wgt_second = loss::search_min_loss_gt_second(wgt_first, loss, ti, snvs, boost_param);

    compute_pred::compute_pred(pred, &wgt_second, genot);

    // create table
    let (table_sum, is_eps) =
        table::calculate_table_eps(&pred, ps, phe, boost_param.eps(), boost_param.boost_type());
    log::debug!("table {:?}", table_sum);
    wgt_second.set_contingency_table(table_sum, is_eps);

    sample_weight::renew_score(scores, pred, &wgt_second);
    sample_weight::renew_ws(
        ws,
        scores,
        phe,
        boost_param.loss_func(),
        boost_param.sample_weight_clip(),
    );
    sample_weight::renew_ps(ps, ws);

    //io_wgt_boost::write_ws(fout, &ws_ti, ti);

    wgt_second
}
*/

pub fn boosting_covs(
    wgts: &mut WgtBoosts,
    dataset: &Dataset,
    scores: &mut SampleScore,
    sample_weight: &mut SampleWeight,
    iteration_start: usize,
    dataset_val: Option<&Dataset>,
    scores_val: &mut SampleScore,
) -> usize {
    // no weighted only
    assert_eq!(iteration_start, 0);

    // FIXME: for now, .is_type_logit(), ps can be constructed for the first iteration, so only able for CovWay::First
    let wgtcovs_logreg =
        regression_cov::logistic_regression_covs(dataset.samples(), iteration_start);
    //let wgtcovs_logreg = regression_cov::logistic_regression_covs_sampleweight(
    //    dataset.samples(),
    //    sample_weight.ps().unwrap(),
    //    iteration_start,
    //);

    log::debug!("wgtcovs_logreg: {:?}", &wgtcovs_logreg);

    let covs = dataset.samples().covs();

    let p = wgtcovs_logreg.len();
    for pi in 0..p {
        let wgtcov = wgtcovs_logreg[pi].clone();
        score::add_score_training(scores, &wgtcov, None, covs);
        wgts.add_wgt(wgtcov);
    }

    if let Some(dataset_val) = dataset_val {
        let covs = dataset_val.samples().covs();

        let p = wgtcovs_logreg.len();
        for pi in 0..p {
            let wgtcov = wgtcovs_logreg[pi].clone();
            score::add_score_training(scores_val, &wgtcov, None, covs);
        }
    }

    log::debug!("after cov");
    let phe = dataset.samples().phe_unwrap();
    sample_weight.renew_sample_weight(scores, phe);

    p
}

pub fn boosting_logit_const(
    dataset: &Dataset,
    scores: &mut SampleScore,
    sample_weight: &mut SampleWeight,
    iteration: usize,
    dataset_val: Option<&Dataset>,
    scores_val: &mut SampleScore,
) -> WgtBoost {
    let wzs = sample_weight.wzs().unwrap();
    let wls = sample_weight.wls().unwrap();
    let (coef, _) = coefficient::calculate_coef_logit_const(&wzs, &wls);
    let model = Model::new_coef(coef);
    let wgt = Wgt::new_const(model);
    let wgt_boost = WgtBoost::new_wgt_iteration(wgt, iteration);

    log::debug!("after const");
    log::debug!("wgt {:?}", wgt_boost);

    score::add_score_training(scores, &wgt_boost, None, None);

    if let Some(_) = dataset_val {
        score::add_score_training(scores_val, &wgt_boost, None, None);
    }

    let phe = dataset.samples().phe_unwrap();

    sample_weight.renew_sample_weight(&scores, phe);

    wgt_boost
}

const ITERATION_NUMBER_SNV_LIMIT: usize = 500_000;

fn run_next_iteration(ti: usize, iteration: Option<IterationNumber>, wgts: &WgtBoosts) -> bool {
    match iteration {
        Some(iteration) => {
            // even when is_monitor, if iteration is specified, stop on iteration
            match iteration {
                IterationNumber::Iteration(x) => ti < x,
                IterationNumber::Snv(x) => {
                    // Snv() might not stop, so determine ti limit
                    if ti > ITERATION_NUMBER_SNV_LIMIT {
                        log::info!(
                            "Iteration number exceeded the limit, so stop the iteration: {}",
                            ti
                        );
                        return false;
                    }
                    // use '<=' since the last snv might be selected in a row. ex. rs1, rs1, rs2, rs2. if '< 2',  rs2 will be selected only once
                    wgts.count_unique_snv() <= x
                }
            }
        }
        None => {
            // None only when is_monitor
            if ti > ITERATION_NUMBER_SNV_LIMIT {
                log::info!(
                    "Iteration number exceeded the limit, so stop the iteration: {}",
                    ti
                );
                false
            } else {
                true
            }
        }
    }
}

fn renew_batch(
    bi: usize,
    batch_way: Option<BatchWay>,
    loss_batch_stop: Option<f64>,
    wgts: &WgtBoosts,
    learning_rate: f64,
) -> bool {
    if bi == 0 {
        // for the first iteration
        return true;
    }

    if batch_way.unwrap().use_comp_loss() {
        if bi > 0 {
            let loss_iter = wgts.last_wgt().unwrap().loss();

            if let Some(loss_iter) = loss_iter {
                if loss_iter > loss_batch_stop.unwrap() {
                    log::debug!(
                        "last loss is larger than loss_batch_stop {} > {}",
                        loss_iter,
                        loss_batch_stop.unwrap()
                    );
                    return true;
                }
            }
            // otherwise, last wgt was cov
        }
    }

    let batch_max_iter = batch_way.unwrap().batch_max_iter(learning_rate);

    if bi >= batch_max_iter {
        log::debug!("bi reached batch_max_iter");
    }

    bi >= batch_max_iter
}

fn run_cov_iteration(wgts: &WgtBoosts, cov_way: Option<CovWay>) -> bool {
    match cov_way {
        None => false,
        Some(cov_way) => match cov_way {
            CovWay::First => wgts.wgts().len() == 0,
            CovWay::Every(_) => {
                unimplemented!("Cannot run weighted regression.");
                /*                 if wgts.wgts_n() == 0 {
                    return true;
                }
                // if last n items are all snv, run cov
                let is_last_n_snv = wgts.is_last_snv(x);
                match is_last_n_snv {
                    None => false,
                    Some(x) => x,
                } */
                // nightly
                //is_last_n_snv.is_some_and(|x| x);
            }
        },
    }
}

/*
// TODO: legacy
fn renew_wgt_skip(
    wgt_skip: &mut HashSet<usize>,
    wgt_last: Option<&WgtBoost>,
    wgt: &WgtBoost,
) -> bool {
    // add wgt to wgt_skip if wgt and wgt_last are the same variant and effect size (or loss) are exactly the same.

    if wgt_last.is_none() {
        log::debug!("First SNV.");
        return false;
    }

    if !(wgt.wgt().kind().snv_index() == wgt_last.unwrap().wgt().kind().snv_index()) {
        log::debug!("SNV is different from the last SNV.");
        return false;
    }

    if let Coef::Score4((s0, s1, s2, sm)) = wgt.wgt().model().coef() {
        if let Coef::Score4((s0l, s1l, s2l, sml)) = wgt_last.unwrap().wgt().model().coef() {
            if eq_float(s0, s0l) && eq_float(s1, s1l) && eq_float(s2, s2l) && eq_float(sm, sml) {
                log::debug!("***SNV and coef are the same as the last SNV, so this SNV is added to skip list.***");
                wgt_skip.insert(wgt.wgt().kind().index_snv().unwrap());
                //wgt_skip.push(wgt.wgt().kind().index_snv().unwrap());
                return true;
            }
        }
        log::debug!("SNV is the same as the last SNV but coef are different.");
    }

    return false;
}
*/

// able to delete after moving
// 1. --resume
// 2. covs in while()
// to boosting_batch()
pub fn boosting<W: std::io::Write>(
    writer: &mut BufWriter<W>,
    boost_param: &BoostParam,
    dataset: &Dataset,
    //writer_loss: &mut BufWriter<W>,
    //fout: Option<&str>, // should Some() when is_resume=T or is_write_loss=T
    //use_adjloss: bool,
    is_resume: bool,
    is_write_loss: bool,
    dout: Option<&DoutParaFile>, // for is_resume or is_write_loss
                                 //dout: Option<&Path>, // for is_resume or is_write_loss
                                 //dloss: Option<&Path>, //is_write_loss: bool, // TODO: make different dir in boosting
) -> Option<(usize, f64)> {
    let start_time = Instant::now();

    let phe = dataset.samples().phe_unwrap();
    let snvs = dataset.snvs();

    let n = dataset.samples().samples_n();
    let mut scores: SampleScore = SampleScore::new(n);

    let mut sample_weight = SampleWeight::new(
        n,
        boost_param.boost_type(),
        boost_param.loss_func(),
        boost_param.sample_weight_clip(),
        boost_param.sample_weight_wls_clip(),
    );
    sample_weight.renew_sample_weight(&scores, phe);

    // make pred u8 -> good idea!!
    // but slow when non-simd
    //let mut pred: Vec<u8> = vec![0u8; n];

    let m = snvs.snvs_n();
    let mut loss = LossStruct::new(boost_param.boost_type(), m);
    let mut wgts = WgtBoosts::new(boost_param.boost_type());

    let t_start;
    //let is_exist_wgt = wgt_boost::io::is_exist_wgt(fout);
    //if is_resume && wgt_boost::io::is_exist_wgt_nonzero(dout.unwrap()) {
    if is_resume && dout.unwrap().is_exist_wgt_nonzero() {
        //if is_exist_wgt & is_resume {
        // assume cov is done when resume
        log::info!("resume");

        wgts = WgtBoosts::new_from_file_dir(dout.unwrap(), boost_param.boost_type());
        // align alleles
        set_wgt_m(wgts.wgts_mut(), snvs);

        for wgt in wgts.wgts().iter() {
            score::add_score_training(
                &mut scores,
                wgt,
                Some(dataset.genot()),
                Some(dataset.samples().covs().unwrap()),
            );
        }

        sample_weight.renew_sample_weight(&scores, phe);

        t_start = wgts.wgts().len();
    } else {
        wgt_boost::io::write_cols(writer, boost_param.boost_type());
        t_start = 0;
    }

    if boost_param.is_dom_rec() {
        unimplemented!()
    }

    log::info!("Start Boosting.");
    let mut ti = t_start;

    //let iteration = boost_param.iteration();

    // index of wgt
    // legacy
    //let wgt_skip: HashSet<usize> = HashSet::new();
    //let mut wgt_skip: HashSet<usize> = HashSet::new();
    //let mut wgt_last: Option<WgtBoost> = None;
    while run_next_iteration(ti, boost_param.iteration(), &wgts) {
        if ti % 100 == 0 {
            log::info!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
        } else {
            log::debug!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
        }

        log::debug!("Unique counts of SNVs: {}", wgts.count_unique_snv());

        if run_cov_iteration(&wgts, boost_param.cov_way()) {
            log::debug!("Run cov");

            let p = boosting_covs(
                &mut wgts,
                dataset,
                &mut scores,
                &mut sample_weight,
                ti,
                None,
                &mut SampleScore::new_empty(),
            );

            wgts.write_wgt_writer(writer);

            ti += p;
        } else {
            log::debug!("Run snv");
            let wgt = boosting_iter_snv(
                ti,
                boost_param,
                &mut loss,
                dataset,
                &mut scores,
                &mut sample_weight,
                None,
                None,
                &mut SampleScore::new_empty(),
            );

            // legacy
            //let skip_this_snv = renew_wgt_skip(&mut wgt_skip, wgt_last.as_ref(), &wgt);
            //// is this ok?
            //if skip_this_snv {
            //    continue;
            //}

            wgts.add_wgt(wgt.clone());

            //log::debug!("wgts {:?}", wgts);
            wgts.write_wgt_writer(writer);

            //wgt_last = Some(wgt);

            ti += 1;
        }

        if is_write_loss {
            if (ti < 200) || (ti % 100 == 0) {
                let mut writer_loss = dout.unwrap().bufwriter_floss(ti);
                //let mut writer_loss = wgt_boost::io::bufwriter_floss(dout.unwrap(), ti);
                loss.write_writer(&mut writer_loss, snvs);
            }
        }

        log::debug!("\n");
    }

    None
}

//fn create_wgt_skip_snvs(
//    wgt_skip_batch: &HashSet<usize>,
//    wgt_skip: &HashSet<usize>,
//) -> HashSet<usize> {
//    let mut wgt_skip_snvs = wgt_skip_batch.clone();
//    wgt_skip_snvs.extend(wgt_skip);
//    //let wgt_skip_snvs: HashSet<usize> = wgt_skip_batch.union(wgt_skip).collect();
//    wgt_skip_snvs
//}

fn create_extract_snvs_batch(
    loss: &mut LossStruct,
    dataset: &Dataset,
    scores: &SampleScore,
    sample_weight: &mut SampleWeight,
    boost_param: &BoostParam,
    extract_snvs: Option<&HashSet<usize>>,
) -> (HashSet<usize>, Option<f64>) {
    let mut loss_decrease_by;

    let mut snv_first_ldfix: Option<usize> = None;

    let loss_here = if let Some(BatchWay::LdFix(..)) = boost_param.batch_way() {
        // select the best snv, and calculate loss again.
        // make snv set of the top n snvs which decreases the loss
        // between the first and the second time.

        let mut scores_local = scores.clone_align();
        let mut sample_weight_local = sample_weight.clone_align();

        let phe = dataset.samples().phe_unwrap();
        sample_weight_local.renew_sample_weight(&scores_local, phe);

        assert_eq!(sample_weight.wls_pad(), sample_weight_local.wls_pad());
        assert_eq!(sample_weight.wzs_pad(), sample_weight_local.wzs_pad());

        loss::calculate_loss_gt(
            loss,
            dataset,
            &mut sample_weight_local,
            boost_param,
            extract_snvs,
            None,
        );
        let loss_first = loss.clone();

        // first iteration
        let mut wgt = loss::search_min_loss_gt(
            loss,
            0, // tmp
            dataset,
            //&sample_weight_local,
            boost_param,
            extract_snvs,
        );
        let genot = dataset.genot();

        snv_first_ldfix = wgt.wgt().kind().index_snv();

        update_wgt::update_coef(
            &mut wgt,
            boost_param,
            &phe,
            genot,
            &sample_weight_local,
            dataset,
        );
        score::add_score_training(&mut scores_local, &wgt, Some(dataset.genot()), None);
        let phe = dataset.samples().phe_unwrap();
        sample_weight_local.renew_sample_weight(&scores_local, phe);

        // second iteration
        loss::calculate_loss_gt(
            loss,
            dataset,
            &mut sample_weight_local,
            boost_param,
            extract_snvs,
            None,
        );

        // the larger, the more disequilibrium
        loss_decrease_by = loss_first.decrease_by(&loss);
        // TODO: implement search_max()
        // flip sign to use search_min()
        loss_decrease_by = loss_decrease_by.flip_sign();

        &mut loss_decrease_by
        // TODO: check if the top wgt is seelcted in snv set
    } else {
        loss::calculate_loss_gt(
            loss,
            dataset,
            sample_weight,
            boost_param,
            extract_snvs,
            None,
        );
        loss
    };

    let m_top = boost_param.batch_way().unwrap().batch_size();
    log::debug!("batch m {}", m_top);

    let (use_snvs, loss_next) = loss_here.search_topprop_n(m_top, extract_snvs);
    log::debug!("loss batch next {}", loss_next);

    // for print
    let (loss_min, _, _) = loss_here.search_min(extract_snvs);
    log::debug!("loss min {}", loss_min);

    let snv_use_batch: HashSet<usize> = HashSet::from_iter(
        use_snvs
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(i, _)| i),
    );

    log::debug!("snv_use_batch {:?}", snv_use_batch);

    if let Some(BatchWay::LdFix(..)) = boost_param.batch_way() {
        if !snv_use_batch.contains(&snv_first_ldfix.unwrap()) {
            panic!("snv_first_ldfix sould be in snv_use_batch in the most case.");
        }
    };

    if boost_param.batch_way().unwrap().use_comp_loss() {
        (snv_use_batch, Some(loss_next))
    } else {
        (snv_use_batch, None)
    }
}

fn load_finitial_interaction(finitial_snvs: &Path) -> Vec<usize> {
    let mut initial_snvs = vec![];

    //.has_headers(true)
    let mut reader = csv::Reader::from_path(finitial_snvs).unwrap();
    for result in reader.records() {
        let record = result.unwrap();
        let snv = record[0].parse::<usize>().unwrap();
        initial_snvs.push(snv);
    }

    initial_snvs
}

fn find_initial_snvs_for_interaction(
    loss: &LossStruct,
    snvs: &Snvs,
    boost_param: &BoostParam,
    finitial_snvs: Option<&Path>,
) -> Vec<usize> {
    match boost_param.batch_interaction_way().unwrap() {
        BatchInteractionWay::InitialLdWidthMaxFixFilter(..) => {
            find_initial_snvs_for_interaction_maxfixfilter(loss, snvs, boost_param)
        }
        BatchInteractionWay::InitialLdWidthRandomFixFilter(..) => {
            load_finitial_interaction(finitial_snvs.unwrap())
        }
        BatchInteractionWay::InitialAllFixFilter(..) => {
            load_finitial_interaction(finitial_snvs.unwrap())
        }
        _ => unimplemented!(),
    }
}

fn find_initial_snvs_for_interaction_maxfixfilter(
    loss: &LossStruct,
    snvs: &Snvs,
    boost_param: &BoostParam,
) -> Vec<usize> {
    let mut extract_snvs_chroms: Vec<usize> = vec![];

    for chrom in Chrom::variants() {
        let use_snvs_chrom = snvs.use_snvs_chrom(&chrom);
        let indexs_chrom = snvs
            .snv_ids()
            .iter()
            .enumerate()
            .filter(|(_, x)| x.chrom() == &chrom)
            .map(|(i, _)| i)
            .collect::<Vec<usize>>();
        let positions_chrom = snvs
            .positions()
            .iter()
            .zip(use_snvs_chrom)
            .filter(|(_, b)| *b)
            .map(|(x, _)| *x)
            .collect::<Vec<usize>>();
        //let positions = snvs.positions(&chrom);

        // TODO for r2
        let ld_width = boost_param.batch_interaction_way().unwrap().ld_width();

        // ld_index of each snv
        let ld_indexs_chrom = positions_chrom
            .iter()
            .map(|x| x / ld_width)
            .collect::<Vec<usize>>();

        // assume positions are sorted
        let mut ld_indexs_uniq = ld_indexs_chrom.clone();
        ld_indexs_uniq.dedup();

        // for snvs with the same ld_index, select the one with the largest loss
        let extract_snvs = ld_indexs_uniq
            .iter()
            .map(|ld_index| {
                // snvs which has same ld
                let mut snvs_same_ld = indexs_chrom
                    .iter()
                    .zip(ld_indexs_chrom.iter())
                    .filter(|(_, ld)| *ld == ld_index)
                    .map(|(index, _)| *index)
                    .collect::<Vec<usize>>();

                // sort by loss
                snvs_same_ld.sort_by(|a, b| {
                    loss.access_loss(*a)
                        .partial_cmp(&loss.access_loss(*b))
                        .unwrap()
                });

                // select the minimum loss
                let snv = snvs_same_ld[0];
                snv
            })
            .collect::<Vec<usize>>();

        extract_snvs_chroms.extend(extract_snvs);
    }

    // TODO: assert uniq and ascending
    assert!(vec::is_sorted::<usize>(&extract_snvs_chroms));

    extract_snvs_chroms
}

// 100M interactions
const INITIAL_FILE_MAX: usize = 100_000_000;
// too large when writing down to file
// 1G interactions
//const INITIAL_FILE_MAX: usize = 1_000_000_000;

// TODO: extrct common with create_extract_snvs_interaction_batch()
fn create_initial_interaction(
    loss: &mut LossStruct,
    dataset: &Dataset,
    sample_weight: &mut SampleWeight,
    boost_param: &BoostParam,
    //extract_snvs_interaction: Vec<(usize, usize)>,
    dout: Option<&DoutParaFile>,
    finitial_snvs: Option<&Path>,
) -> Vec<(usize, usize)> {
    // first, calculate loss of single snvs
    // This uses cov adjustment

    let mut alphas_single = vec![f64::NAN; dataset.snvs().snvs_n()];

    loss::calculate_loss_gt(
        loss,
        dataset,
        sample_weight,
        boost_param,
        None,
        Some(&mut alphas_single),
    );

    let mut writer_loss = dout.unwrap().bufwriter_floss_initial_single();
    loss.write_initial_single_writer(&mut writer_loss, &alphas_single);

    let m = dataset.snvs().snvs_n();

    // next, create extract_snvs_interaction for all snvs
    let extract_snvs =
        find_initial_snvs_for_interaction(loss, dataset.snvs(), boost_param, finitial_snvs);

    log::debug!("extract_snvs {:?}", extract_snvs.len());
    match boost_param.batch_interaction_way().unwrap() {
        BatchInteractionWay::InitialAllFixFilter(..) => {
            log::debug!("initial interaction will be {:?}", extract_snvs.len() * m);
        }
        _ => {
            log::debug!(
                "initial interaction will be {:?}",
                extract_snvs.len() * extract_snvs.len()
            );
        }
    };

    if !boost_param.batch_interaction_way().unwrap().is_filter_all() {
        panic!(
            "batch_interaction_way should be filter_all due to too large extract_snvs_interaction."
        );
    }

    let filter_size_all = boost_param
        .batch_interaction_way()
        .unwrap()
        .filter_size_all();

    let snvs_num_each_file = match boost_param.batch_interaction_way().unwrap() {
        BatchInteractionWay::InitialAllFixFilter(..) => INITIAL_FILE_MAX / m,
        _ => INITIAL_FILE_MAX / extract_snvs.len(),
    };
    log::debug!("snvs_num_each_file {:?}", snvs_num_each_file);

    // for each file
    let mut loss_here = loss.clone();
    let mut alphas = vec![];

    // corresponding to loss
    let mut extract_snvs_interaction = vec![];
    for (file_i, extract_snvs_i) in extract_snvs.chunks(snvs_num_each_file).enumerate() {
        log::debug!("initial file_i {:?}", file_i);
        log::debug!("loss single {:?}", loss.inner_single().len());
        log::debug!("loss {:?}", loss.inner_interaction().len());
        log::debug!("loss_here single {:?}", loss.inner_single().len());
        log::debug!("loss_here {:?}", loss.inner_interaction().len());

        // TODO: load file if exists

        let extract_snvs_interaction_file_i = match boost_param.batch_interaction_way().unwrap() {
            BatchInteractionWay::InitialAllFixFilter(..) => {
                // extracted for file_i x all
                extract_snvs_i
                    .iter()
                    .flat_map(|i| (0..m).map(|j| (*i, j)))
                    .collect::<Vec<(usize, usize)>>()
            }
            _ => {
                // extracted for file_i x all extracted
                extract_snvs_i
                    .iter()
                    .flat_map(|i| extract_snvs.iter().map(|j| (*i, *j)))
                    .collect::<Vec<(usize, usize)>>()
            }
        };

        loss_here.resize_interaction(extract_snvs_interaction_file_i.len());
        alphas.resize(extract_snvs_interaction_file_i.len(), f64::NAN);

        let floss_initial = dout.unwrap().floss_initial(file_i, finitial_snvs);
        let mut floss_initial_gz = floss_initial.clone();
        // TODO: add .gz only?
        floss_initial_gz.set_extension("loss.gz");
        //log::debug!("floss_initial {:?}", floss_initial);
        //log::debug!("floss_initial_gz {:?}", floss_initial_gz);

        if floss_initial.exists() || floss_initial_gz.exists() {
            // if .loss or .loss.gz exists
            log::debug!("skip file_i {:?}", file_i);

            // BUG: loss.serach_topprop should be done
            continue;
        } else {
            //log::debug!("do not skip file_i {:?}", file_i);

            loss::calc::calc_loss_logit_interaction_cont_table(
                &mut loss_here,
                dataset,
                boost_param,
                &extract_snvs_interaction_file_i,
                &mut alphas,
            );
        }

        // write to file
        let mut writer_loss = dout.unwrap().bufwriter_floss_initial(file_i, finitial_snvs);
        loss_here.write_initial_writer(&mut writer_loss, &extract_snvs_interaction_file_i, &alphas);

        // merge loss_here to loss
        loss.extend_interaction(&loss_here);
        extract_snvs_interaction.extend(extract_snvs_interaction_file_i);

        // select top n snvs with the largest loss
        let (use_snvs, _) = loss.search_topprop_interaction_n(filter_size_all);
        loss.extract_interaction(&use_snvs);
        extract_snvs_interaction = extract_snvs_interaction
            .iter()
            .zip(use_snvs.iter())
            .filter(|(_, b)| **b)
            .map(|(i, _)| *i)
            .collect::<Vec<(usize, usize)>>();

        //extract_snvs_interaction.extend(extract_snvs_interaction_filter);
        // do not save alphas
    }

    let (use_snvs, _) = loss.search_topprop_interaction_n(filter_size_all);
    let extract_snvs_interaction_filter = extract_snvs_interaction
        .iter()
        .zip(use_snvs.iter())
        .filter(|(_, b)| **b)
        .map(|(i, _)| *i)
        .collect::<Vec<(usize, usize)>>();
    //} else {
    //extract_snvs_interaction
    //};

    extract_snvs_interaction_filter
}

// batch_interaction_way is required
fn create_extract_snvs_interaction_batch(
    loss: &mut LossStruct,
    dataset: &Dataset,
    _scores: &SampleScore,
    sample_weight: &SampleWeight,
    boost_param: &BoostParam,
    extract_snvs_interaction: Vec<(usize, usize)>,
) -> (Vec<(usize, usize)>, Option<f64>, Vec<(usize, usize)>) {
    //let mut loss_decrease_by;

    // TODO: check if extract_snvs_interaction is unique

    let mut _snv_first_ldfix: Option<usize> = None;

    let loss_here = if let Some(BatchWay::LdFix(..)) = boost_param.batch_way() {
        unimplemented!("ny");

        // select the best snv, and calculate loss again.
        // make snv set of the top n snvs which decreases the loss
        // between the first and the second time.

        //let mut scores_local = scores.clone_align();
        ////let mut scores_local = scores.to_vec();
        //let mut sample_weight_local = sample_weight.clone_align();
        ////let mut sample_weight_local = sample_weight.clone();

        //let phe = dataset.samples().phe_unwrap();
        //sample_weight_local.renew_sample_weight(&scores_local, phe);

        //assert_eq!(sample_weight.wls_pad(), sample_weight_local.wls_pad());
        //assert_eq!(sample_weight.wzs_pad(), sample_weight_local.wzs_pad());

        //loss::calculate_loss_gt(
        //    loss,
        //    dataset,
        //    &sample_weight_local,
        //    boost_param,
        //    extract_snvs_interaction,
        //);
        //let loss_first = loss.clone();

        //// first iteration
        //let mut wgt = loss::search_min_loss_gt(
        //    loss,
        //    0, // tmp
        //    dataset,
        //    //&sample_weight_local,
        //    boost_param,
        //    extract_snvs_interaction,
        //);
        ////let mut pred = vec![0u8; dataset.samples().samples_n()];
        //let genot = dataset.genot();
        ////let phe = dataset.samples().phe_unwrap();
        ////compute_pred::compute_pred(&mut pred, &wgt, genot);

        //snv_first_ldfix = wgt.wgt().kind().index_snv();

        //// TODO: move
        ////let mut sample_weight_local = sample_weight.clone();
        ////let mut sample_weight_local = create_clone_sample_weight(sample_weight_here);
        //update_wgt::update_coef(
        //    &mut wgt,
        //    boost_param,
        //    //&mut pred,
        //    &phe,
        //    genot,
        //    &sample_weight_local,
        //    dataset,
        //);
        //score::add_score_training(&mut scores_local, &wgt, Some(dataset.genot()), None);
        ////sample_weight::renew_score(&mut scores_local, &pred, &wgt);
        //let phe = dataset.samples().phe_unwrap();
        //sample_weight_local.renew_sample_weight(&scores_local, phe);

        //// second iteration
        //loss::calculate_loss_gt(
        //    loss,
        //    dataset,
        //    &sample_weight_local,
        //    boost_param,
        //    extract_snvs_interaction,
        //);

        //// the larger, the more disequilibrium
        //loss_decrease_by = loss_first.decrease_by(&loss);
        //// TODO: implement search_max()
        //// flip sign to use search_min()
        //loss_decrease_by = loss_decrease_by.flip_sign();

        //&mut loss_decrease_by
        //// TODO: check if the top wgt is seelcted in snv set
    } else {
        // calculating interactions only is fine
        loss::calc::calc_loss_logit_interaction(
            loss,
            dataset,
            sample_weight,
            boost_param,
            &extract_snvs_interaction,
        );

        loss
    };
    log::debug!("done calculate loss.");

    log::debug!("batch way {:?}", boost_param.batch_interaction_way());
    let m_top = boost_param.batch_interaction_way().unwrap().batch_size();
    log::debug!("batch m_top {}", m_top);

    let (use_snvs, loss_next) = loss_here.search_topprop_interaction_n(m_top);
    log::debug!("loss batch next {}", loss_next);

    // for print
    let (loss_min, _, _) = loss_here.search_min_interaction(&extract_snvs_interaction);
    log::debug!("loss min {}", loss_min);

    let snvs_use_batch = extract_snvs_interaction
        .iter()
        .zip(use_snvs.iter())
        .filter(|(_, b)| **b)
        .map(|(i, _)| *i)
        .collect::<Vec<(usize, usize)>>();

    log::debug!("snvs_use_batch {:?}", snvs_use_batch.len());

    // filter snvs_interaction_all
    let extract_snvs_interaction_filter =
        if boost_param.batch_interaction_way().unwrap().is_filter_all() {
            let filter_size_all = boost_param
                .batch_interaction_way()
                .unwrap()
                .filter_size_all();

            let (use_snvs, _) = loss_here.search_topprop_interaction_n(filter_size_all);
            extract_snvs_interaction
                .iter()
                .zip(use_snvs.iter())
                .filter(|(_, b)| **b)
                .map(|(i, _)| *i)
                .collect::<Vec<(usize, usize)>>()
        } else {
            extract_snvs_interaction
        };

    // old?
    //if let Some(BatchWay::LdFix(..)) = boost_param.batch_way() {
    //    //if !wgt_use_batch.contains(&snv_first_ldfix.unwrap()) {
    //    //    panic!("snv_first_ldfix sould be in wgt_use_batch in the most case.");
    //    //}
    //};

    return (
        snvs_use_batch,
        Some(loss_next),
        extract_snvs_interaction_filter,
    );
}

fn filter_interaction_batch(
    snvs_interaction_batch: Vec<(usize, usize)>,
    loss: &LossStruct,
    boost_param: &BoostParam,
) -> Vec<(usize, usize)> {
    let filter_size_batch = boost_param
        .batch_interaction_way()
        .unwrap()
        .filter_size_batch();

    let (use_snvs, _) = loss.search_topprop_interaction_n(filter_size_batch);
    let snvs_interaction_batch_extract = snvs_interaction_batch
        .iter()
        .zip(use_snvs.iter())
        .filter(|(_, b)| **b)
        .map(|(i, _)| *i)
        .collect::<Vec<(usize, usize)>>();

    snvs_interaction_batch_extract
}

///// allow input file and stdio for test
///// return (nsnvs, acc)
//pub fn boosting_batch_old<W: std::io::Write>(
//    writer: &mut BufWriter<W>,
//    boost_param: BoostParam,
//    dataset: &Dataset,
//    dataset_val: Option<&Dataset>,
//    use_adjloss: bool,
//    use_const_for_loss: bool,
//    is_resume: bool,
//    is_write_loss: bool,
//    dout: Option<&Path>, // for is_resume or is_write_loss
//    //dloss: Option<&Path>, //is_write_loss: bool, // TODO: make different dir in boosting
//    is_monitor: bool,
//    //nsnvs_monitor: Option<Vec<usize>>,
//    nsnvs_monitor: Option<&[usize]>,
//) -> Option<(usize, f64)> {
//    if !boost_param.cov_way().unwrap().is_first() {
//        panic!("Should indicate CovWay::First in boosting_batch()")
//    }
//    if use_const_for_loss {
//        panic!("Deprecated: use_const_for_loss");
//    }
//
//    let start_time = Instant::now();
//
//    let phe = dataset.samples().phe_unwrap();
//    let genot = dataset.genot();
//    let snvs = dataset.snvs();
//
//    let n = dataset.samples().samples_n();
//    let mut scores: Vec<f64> = vec![0.0; n];
//    let n_val = if let Some(dataset_val) = dataset_val {
//        dataset_val.samples().samples_n()
//    } else {
//        0
//    };
//    let mut scores_val: Vec<f64> = vec![0.0; n_val];
//
//    let mut sample_weight = SampleWeight::new(
//        n,
//        boost_param.boost_type(),
//        boost_param.loss_func(),
//        boost_param.sample_weight_clip(),
//        boost_param.sample_weight_wls_clip(),
//    );
//    sample_weight.renew_sample_weight(&scores, phe);
//
//    // make pred u8 -> good idea!!
//    // but slow when non-simd
//    let mut pred: Vec<u8> = vec![0u8; n];
//    let mut pred_val: Vec<u8> = vec![0u8; n_val];
//
//    let m = dataset.snvs().snvs_n();
//
//    //let nsnvs_monitor = match nsnvs_monitor {
//    //    Some(x) => x,
//    //    None => {
//    //        //let nsnvs_monitor = [
//    //        //    20usize, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 300,
//    //        //    400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
//    //        //    10000,
//    //        //];
//    //        // if exceed 20k, run every 1k
//    //        let nsnvs_monitor = [
//    //            5usize, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 300,
//    //            400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
//    //            10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000,
//    //        ];
//    //        nsnvs_monitor.to_vec()
//    //    }
//    //};
//    //let nsnvs_monitor = [
//    //    5usize, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 300, 400,
//    //    500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
//    //];
//    let mut acc_monitor: BTreeMap<usize, f64> = BTreeMap::new();
//    //let mut acc_monitor: HashMap<usize, f64> = HashMap::new();
//    //let mut acc_monitor = vec![f64::NAN; nsnvs_monitor.len()];
//    // save sample score of previous iteration
//    // for duplicated snvs; to get sample score with 5 snvs, only way is wait until 6-th snv and calculate score of previous snvs
//    let mut scores_val_monitor: Vec<f64> = vec![0.0; n_val];
//
//    let mut loss = LossStruct::new(boost_param.boost_type(), m);
//    let mut wgts = WgtBoosts::new(boost_param.boost_type());
//
//    let t_start;
//    if is_resume & wgt_boost::io::is_exist_wgt_nonzero(dout.unwrap()) {
//        // assume cov is done when resume
//        log::info!("resume");
//
//        panic!("Do not use resume since reesults will be changed.")
//
//        /*
//        wgts = WgtBoosts::new_from_file_dir(dout.unwrap(), boost_param.boost_type());
//        set_wgt_m(wgts.wgts_mut(), snvs);
//        for wgt in wgts.wgts().iter() {
//            score::add_score(
//                &mut scores,
//                wgt,
//                dataset.genot(),
//                Some(dataset.samples().covs().unwrap()),
//            );
//        }
//        sample_weight.renew_sample_weight(&scores, phe);
//        t_start = wgts.wgts().len();
//         */
//    } else {
//        wgt_boost::io::write_cols(writer, boost_param.boost_type());
//        t_start = 0;
//    }
//
//    if boost_param.is_dom_rec() {
//        unimplemented!()
//    }
//
//    log::info!("Start Boosting.");
//    let mut ti = t_start;
//
//    let iteration = boost_param.iteration();
//
//    // first cov
//    {
//        let p = boosting_covs(
//            &mut wgts,
//            dataset,
//            &mut scores,
//            &mut sample_weight,
//            ti,
//            dataset_val,
//            &mut scores_val,
//        );
//        wgts.write_wgt_writer(writer);
//        ti += p;
//    }
//
//    let scores_val_cov = scores_val.clone();
//
//    /*     // for const
//       if use_const_for_loss {
//           let wgt = boosting_logit_const(
//               dataset,
//               &mut scores,
//               &mut sample_weight,
//               ti,
//               dataset_val,
//               &mut scores_val,
//           );
//           wgts.add_wgt(wgt.clone());
//           wgts.write_wgt_writer(writer);
//
//           ti += 1;
//       }
//    */
//    // TODO: cleaner for iteration and integrate boosting()
//    // while b {
//    //       if b {
//    //                   create_wgt_skip_batch
//    //        }
//    //       if b {
//    //                    iteration
//    //        }
//    // }
//
//    // legacy
//    // index of wgt
//    //let mut wgt_skip_snvs: HashSet<usize> = HashSet::new();
//    // for illegal snvs
//    let wgt_skip: HashSet<usize> = HashSet::new();
//    //let mut wgt_skip: HashSet<usize> = HashSet::new();
//    //let mut wgt_last: Option<WgtBoost> = None;
//    // while true {
//    //while run_next_iteration_monitor(ti, boost_param.iteration(), &wgts, is_monitor) {
//    while run_next_iteration(ti, boost_param.iteration(), &wgts) {
//        log::debug!("New batch");
//        let (wgt_skip_batch, loss_batch_stop) = create_wgt_skip_batch(
//            &mut loss,
//            genot,
//            &mut sample_weight,
//            phe,
//            boost_param,
//            &wgt_skip,
//            use_adjloss,
//        );
//
//        log::debug!("Batch wgt size {}", m - wgt_skip_batch.len());
//
//        // iteration in batch
//        let mut bi = 0usize;
//
//        while run_next_iteration_in_batch(
//            bi,
//            boost_param.batch_way(),
//            loss_batch_stop,
//            &wgts,
//            boost_param.learning_rate(),
//        ) & run_next_iteration(ti, iteration, &wgts)
//        //) & run_next_iteration_monitor(ti, iteration, &wgts, is_monitor)
//        {
//            if ti % 100 == 0 {
//                log::info!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
//            } else {
//                log::debug!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
//            }
//
//            log::debug!("Batch iteration: {}", bi);
//            log::debug!("Unique counts of SNVs: {}", wgts.count_unique_snv());
//
//            //let wgt_skip_iter = create_wgt_skip_snvs(&wgt_skip_batch, &wgt_skip);
//            //log::debug!(
//            //    "Batch wgt_batch size, illegal size {} {}",
//            //    m - wgt_skip_batch.len(),
//            //    wgt_skip.len()
//            //);
//            //log::debug!("Batch wgt_iter size {}", m - wgt_skip_iter.len());
//
//            // here before w, z are renewed
//            //let loss_max_theory = compute_loss_max_theory(&zs, &wls, phe.n());
//            //log::debug!("loss max {}", loss_max_theory);
//
//            // score before adding the next wgt
//            scores_val_monitor.clone_from_slice(&scores_val);
//
//            let wgt = boosting_iter_snv(
//                ti,
//                boost_param,
//                &mut loss,
//                genot,
//                &mut pred,
//                &mut scores,
//                &mut sample_weight,
//                phe,
//                snvs,
//                &wgt_skip_batch,
//                //&wgt_skip_iter,
//                use_adjloss,
//                dataset_val,
//                &mut scores_val,
//                &mut pred_val,
//            );
//
//            // legacy
//            //let skip_this_snv = renew_wgt_skip(&mut wgt_skip, wgt_last.as_ref(), &wgt);
//            //// is this ok?
//            //if skip_this_snv {
//            //    continue;
//            //}
//
//            wgts.add_wgt(wgt.clone());
//
//            //log::debug!("wgts {:?}", wgts);
//            wgts.write_wgt_writer(writer);
//
//            //wgt_last = Some(wgt);
//
//            if is_write_loss {
//                if (ti < 200) || (ti % 100 == 0) {
//                    log::debug!("Write loss file.");
//
//                    if use_adjloss {
//                        let mut writer_loss =
//                            wgt_boost::io::bufwriter_floss_adjmax(dout.unwrap(), ti);
//                        loss.write_writer(&mut writer_loss, snvs);
//                    } else {
//                        let mut writer_loss = wgt_boost::io::bufwriter_floss(dout.unwrap(), ti);
//                        loss.write_writer(&mut writer_loss, snvs);
//
//                        fn compute_loss_max_theory(zs: &[f64], wls: &[f64]) -> f64 {
//                            zs.iter().zip(wls.iter()).map(|(z, w)| w * z * z).sum()
//                        }
//
//                        let mut loss_diff = loss.clone();
//                        let loss_diff_mut = loss_diff.inner_mut();
//                        let zs = sample_weight.zs().unwrap();
//                        let wls = sample_weight.wls().unwrap();
//                        let loss_max_theory = compute_loss_max_theory(&zs, &wls);
//                        loss_diff_mut
//                            .iter_mut()
//                            .for_each(|x| *x = 0.5 * (*x - loss_max_theory));
//
//                        let mut writer_loss =
//                            wgt_boost::io::bufwriter_floss_adjmax(dout.unwrap(), ti);
//                        loss_diff.write_writer(&mut writer_loss, snvs);
//                    }
//                }
//            }
//
//            bi += 1;
//            ti += 1;
//
//            /*             // for const
//            // below of ti+=1
//            if use_const_for_loss {
//                let wgt = boosting_logit_const(
//                    dataset,
//                    &mut scores,
//                    &mut sample_weight,
//                    ti,
//                    dataset_val,
//                    &mut scores_val,
//                );
//                wgts.add_wgt(wgt.clone());
//                wgts.write_wgt_writer(writer);
//
//                ti += 1;
//            }
//            */
//            //if let Some(dataset_val) = dataset_val {
//            if is_monitor {
//                let dataset_val = dataset_val.unwrap();
//                let nsnv = wgts.count_unique_snv();
//                // -1 for ?? -> should be legacy bug when using ti
//                //let nsnv_monitor = nsnv - 1;
//                if is_update_acc(nsnv, nsnvs_monitor) {
//                    //if nsnvs_monitor.contains(&nsnv_monitor) {
//                    log::debug!("monitor stat: nsnvs {}", wgts.count_unique_snv());
//                    update_acc_monitor(
//                        &mut acc_monitor,
//                        nsnv,
//                        dataset_val,
//                        &scores_val_monitor,
//                        //&scores_val,
//                        &scores_val_cov,
//                    );
//
//                    // write acc
//                    let mut writer_acc = wgt_boost::io::bufwriter_acc(dout.unwrap());
//                    write_acc(&acc_monitor, &mut writer_acc);
//
//                    let is_exit = monitor_acc(&acc_monitor, nsnvs_monitor, nsnv);
//                    //let is_exit = monitor_acc(&mut acc_monitor, nsnvs_monitor, nsnv_monitor);
//                    if is_exit {
//                        log::info!("Exit since monitoring accuracy exceeded the threshold.");
//                        // TODO: create wgt with maximum acc
//                        let (nsnvs_max, acc_max) = compute_acc_max(&mut acc_monitor);
//                        //compute_acc_max(&mut acc_monitor, &nsnvs_monitor);
//                        return Some((nsnvs_max, acc_max));
//                    }
//                }
//            };
//
//            log::debug!("\n");
//        }
//    }
//
//    return None;
//}

/// allow input file and stdio for test
/// return (nsnvs, acc)
pub fn boosting_batch<W: std::io::Write>(
    writer: &mut BufWriter<W>,
    boost_param: &BoostParam,
    dataset: &Dataset,
    dataset_val: Option<&Dataset>,
    is_resume: bool,
    is_write_loss: bool,
    dout: Option<&DoutParaFile>, // for is_resume or is_write_loss
    //dout: Option<&Path>, // for is_resume or is_write_loss
    //dloss: Option<&Path>, //is_write_loss: bool, // TODO: make different dir in boosting
    //is_monitor: bool,
    // now in boost_param
    //nsnvs_monitor: Option<Vec<usize>>,
    //nsnvs_monitor: Option<&[usize]>,
    fscore_start: Option<&Path>,
) -> Option<(usize, f64)> {
    if !boost_param.cov_way().unwrap().is_first() {
        panic!("Should indicate CovWay::First in boosting_batch()")
    }

    let start_time = Instant::now();

    let n = dataset.samples().samples_n();
    let n_val = if let Some(dataset_val) = dataset_val {
        dataset_val.samples().samples_n()
    } else {
        0
    };
    let (mut scores, mut scores_val) = if let Some(fscore_start) = fscore_start {
        SampleScore::new_from(
            fscore_start,
            dataset.samples().names(),
            dataset_val.map(|x| x.samples().names()),
        )
    } else {
        let scores = SampleScore::new(n);
        let scores_val = SampleScore::new(n_val);
        (scores, scores_val)
    };
    assert_eq!(scores_val.n(), n_val);

    let phe = dataset.samples().phe_unwrap();
    let snvs = dataset.snvs();

    //let mut scores = SampleScore::new(n);
    //let n_val = if let Some(dataset_val) = dataset_val {
    //    dataset_val.samples().samples_n()
    //} else {
    //    0
    //};
    //let mut scores_val = SampleScore::new(n_val);

    let mut sample_weight = SampleWeight::new(
        n,
        boost_param.boost_type(),
        boost_param.loss_func(),
        boost_param.sample_weight_clip(),
        boost_param.sample_weight_wls_clip(),
    );
    sample_weight.renew_sample_weight(&scores, phe);

    let mut acc_monitor: BTreeMap<usize, f64> = BTreeMap::new();

    // save sample score of previous iteration
    // for duplicated snvs; to get sample score with 5 snvs, only way is wait until 6-th snv and calculate score of previous snvs
    let mut scores_val_monitor;

    let m = dataset.snvs().snvs_n();
    let mut loss = LossStruct::new(boost_param.boost_type(), m);
    let mut wgts = WgtBoosts::new(boost_param.boost_type());

    let t_start;
    //if is_resume && wgt_boost::io::is_exist_wgt_nonzero(dout.unwrap()) {
    if is_resume && dout.unwrap().is_exist_wgt_nonzero() {
        // assume cov is done when resume
        log::info!("resume");

        panic!("Do not use resume since results will be changed depending on batch.")

        /*
        wgts = WgtBoosts::new_from_file_dir(dout.unwrap(), boost_param.boost_type());
        set_wgt_m(wgts.wgts_mut(), snvs);
        for wgt in wgts.wgts().iter() {
            score::add_score(
                &mut scores,
                wgt,
                dataset.genot(),
                Some(dataset.samples().covs().unwrap()),
            );
        }
        sample_weight.renew_sample_weight(&scores, phe);
        t_start = wgts.wgts().len();
         */
    } else {
        wgt_boost::io::write_cols(writer, boost_param.boost_type());
        t_start = 0;
    }

    if boost_param.is_dom_rec() {
        unimplemented!()
    }

    log::info!("Start Boosting.");
    let mut ti = t_start;

    // first cov: not in while loop to save scores_val_cov.
    {
        let p = boosting_covs(
            &mut wgts,
            dataset,
            &mut scores,
            &mut sample_weight,
            ti,
            dataset_val,
            &mut scores_val,
        );
        wgts.write_wgt_writer(writer);
        ti += p;
    }

    let scores_val_cov = scores_val.clone_align();

    // for monitor acc
    let mut nsnv_prev = 0usize;

    // iteration in batch
    let mut bi = 0usize;

    let mut loss_batch_stop: Option<f64> = None;
    let mut snvs_batch = HashSet::new();

    while run_next_iteration(ti, boost_param.iteration(), &wgts) {
        if ti % 100 == 0 {
            log::info!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
        } else {
            log::debug!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
        }

        if renew_batch(
            bi,
            boost_param.batch_way(),
            loss_batch_stop,
            &wgts,
            boost_param.learning_rate(),
        ) {
            log::debug!("New batch");
            (snvs_batch, loss_batch_stop) = create_extract_snvs_batch(
                &mut loss,
                dataset,
                &scores,
                &mut sample_weight,
                boost_param,
                None,
            );
            bi = 0;

            log::debug!("Batch wgt size {}", snvs_batch.len());
        }

        log::debug!("Batch iteration: {}", bi);
        log::debug!("Unique counts of SNVs: {}", wgts.count_unique_snv());

        // score before adding the next wgt
        scores_val_monitor = scores_val.clone_align();

        let wgt = boosting_iter_snv(
            ti,
            boost_param,
            &mut loss,
            dataset,
            &mut scores,
            &mut sample_weight,
            Some(&snvs_batch),
            dataset_val,
            &mut scores_val,
        );

        wgts.add_wgt(wgt.clone());

        //log::debug!("wgts {:?}", wgts);
        wgts.write_wgt_writer(writer);

        if is_write_loss {
            if (ti < 200) || (ti % 100 == 0) {
                log::debug!("Write loss file.");

                let mut writer_loss = dout.unwrap().bufwriter_floss_adjmax(ti);
                //let mut writer_loss = wgt_boost::io::bufwriter_floss_adjmax(dout.unwrap(), ti);
                loss.write_writer(&mut writer_loss, snvs);
            }
        }

        bi += 1;
        ti += 1;

        //if let Some(dataset_val) = dataset_val {
        //if is_monitor {
        if boost_param.is_monitor() {
            let dataset_val = dataset_val.unwrap();
            let nsnv = wgts.count_unique_snv();

            // To deal with duplicated snvs, monitor accuracy with the previous snv
            // only monitor acc when nsnv is changed
            // ex. To monitor acc of #snv=5, wait until nsnv=6 and calculate acc of previous iteration
            // This is when nsnv(=6) is changed from nsnv_prev(=5) and nsnv_prev(=5) is in nsnvs_monitor
            //
            // score_val_monitor is score of the previous iteration

            //if (nsnv != nsnv_prev) && is_update_acc(nsnv, nsnvs_monitor) {
            //log::debug!("nsnv, nsnv_prev: {}, {}", nsnv, nsnv_prev);
            if (nsnv != nsnv_prev) && is_update_acc(nsnv_prev, boost_param.monitor_nsnvs()) {
                //log::debug!("monitor stat: nsnvs {}", wgts.count_unique_snv());
                log::debug!("monitor stat: nsnv_prev {}", nsnv_prev);
                update_acc_monitor(
                    &mut acc_monitor,
                    //nsnv,
                    nsnv_prev,
                    dataset_val,
                    &scores_val_monitor,
                    &scores_val_cov,
                    boost_param.acc_metric(),
                );

                // write acc
                let mut writer_acc = dout.unwrap().bufwriter_facc();
                //let mut writer_acc = wgt_boost::io::bufwriter_acc(dout.unwrap());
                write_acc(&acc_monitor, &mut writer_acc);

                let is_exit = monitor_acc(&acc_monitor, boost_param.monitor_nsnvs(), nsnv_prev);
                //let is_exit = monitor_acc(&mut acc_monitor, nsnvs_monitor, nsnv_monitor);
                if is_exit {
                    log::info!("Exit since monitoring accuracy exceeded the threshold.");
                    // TODO: create wgt with maximum acc
                    let (nsnvs_max, acc_max) = compute_acc_max(&mut acc_monitor);
                    //compute_acc_max(&mut acc_monitor, &nsnvs_monitor);
                    return Some((nsnvs_max, acc_max));
                }
            }

            nsnv_prev = nsnv;
        };

        log::debug!("\n");
    }

    return None;
}

//fn create_one_by_selected(wgts: &WgtBoosts) -> Vec<(usize, usize)> {
//    let last_wgt = wgts.last_wgt();
//    if last_wgt.is_none() {
//        // first iteration
//        return vec![];
//    }
//    let last_wgt = last_wgt.unwrap().wgt();
//    //if !last_wgt.is_snv() {
//    if !last_wgt.is_snv_single() {
//        // last wgt is cov or interaction
//        return vec![];
//    }
//    let last_snv_index = last_wgt.kind().index_snv().unwrap();
//
//    // extract single snvs
//    let selected_snv_indexs = wgts
//        .wgts()
//        .iter()
//        .filter(|w| w.wgt().is_snv_single())
//        .map(|w| w.kind().index_snv().unwrap())
//        .collect::<Vec<usize>>();
//    // sort?
//
//    // create interaction
//    // remove self by self
//    let snvs_interaction = selected_snv_indexs
//        .iter()
//        .filter(|i| **i != last_snv_index)
//        .map(|i| (*i, last_snv_index))
//        .collect::<Vec<(usize, usize)>>();
//    snvs_interaction
//}

fn last_interaction_selected_by_selected(
    wgts: &WgtBoosts,
    //snvs_interaction: Vec<(usize, usize)>,
) -> Vec<(usize, usize)> {
    // prepare last snv
    let last_wgt = wgts.last_wgt();
    if last_wgt.is_none() {
        return vec![];
    }
    let last_wgt = last_wgt.unwrap().wgt();
    //if !last_wgt.is_snv() {
    if !last_wgt.is_snv_single() {
        return vec![];
    }
    let last_snv_index = last_wgt.kind().index_snv().unwrap();

    // extract single snvs
    let selected_snv_indexs = wgts
        .wgts()
        .iter()
        .filter(|w| w.wgt().is_snv_single())
        .map(|w| w.kind().index_snv().unwrap())
        .collect::<Vec<usize>>();

    // create interaction
    // remove self by self
    let last_snvs_interaction = selected_snv_indexs
        .iter()
        .filter(|i| **i != last_snv_index)
        .map(|i| (*i, last_snv_index))
        .collect::<Vec<(usize, usize)>>();

    //// add
    //let mut snvs_interaction_both = snvs_interaction;
    //snvs_interaction_both.extend(snvs_interaction_last);

    //snvs_interaction_both

    last_snvs_interaction
}

//fn create_selected_by_selected(wgts: &WgtBoosts) -> Vec<(usize, usize)> {
//    // extract single snvs
//    let selected_snv_indexs = wgts
//        .wgts()
//        .iter()
//        .filter(|w| w.wgt().is_snv_single())
//        .map(|w| w.kind().index_snv().unwrap())
//        .collect::<Vec<usize>>();
//
//    // create interaction
//    // remove self by self
//    let snvs_interaction = selected_snv_indexs
//        .iter()
//        .flat_map(|i| {
//            selected_snv_indexs
//                .iter()
//                .filter(|j| *i != **j)
//                .map(|j| (*i, *j))
//        })
//        .collect::<Vec<(usize, usize)>>();
//
//    snvs_interaction
//}

fn filter_snvs_interaction(
    snvs_interaction: &[(usize, usize)],
    //boost_param: &BoostParam,
    ld_criteria: LdCriteria,
    //ld_radius: usize,
    snvs: &Snvs,
    genot: &Genot,
) -> Vec<(usize, usize)> {
    let snv_ids = snvs.snv_ids();
    //let ld_radius = boost_param.ld_radius().unwrap();
    let snvs_interaction = snvs_interaction
        .iter()
        .filter(|(i, j)| {
            !snv_ids[*i].is_in_ld_criteria(
                &snv_ids[*j],
                ld_criteria,
                &genot.to_genot_snv(*i),
                &genot.to_genot_snv(*j),
            )
        })
        .copied()
        .collect::<Vec<(usize, usize)>>();

    snvs_interaction
}

/// Assume this fn is called at every iteration.
/// add interaction term of the last wgt by others.
/// interaction excludes self by self
fn add_last_snvs_interaction(
    wgts: &WgtBoosts,
    boost_param: &BoostParam,
    snvs: &Snvs,
    genot: &Genot,
    snvs_interaction: Vec<(usize, usize)>,
) -> Vec<(usize, usize)> {
    if let Some(interaction_way) = boost_param.interaction_way() {
        let last_snvs_interaction = match interaction_way {
            InteractionWay::SelectedBySelected => last_interaction_selected_by_selected(wgts),
            InteractionWay::SelectedBySelectedLd => {
                let last_snvs_interaction = last_interaction_selected_by_selected(wgts);

                let last_snvs_interaction = filter_snvs_interaction(
                    &last_snvs_interaction,
                    boost_param.ld_criteria().unwrap(),
                    //boost_param.ld_radius().unwrap(),
                    snvs,
                    genot,
                );
                last_snvs_interaction
            }
        };

        // add
        let last_snvs_interaction_uniq =
            extract_snvs_interaction_not_in_vec(last_snvs_interaction, &snvs_interaction);

        let mut snvs_interaction = snvs_interaction;
        snvs_interaction.extend(last_snvs_interaction_uniq);
        //snvs_interaction.extend(last_snvs_interaction);
        snvs_interaction
    } else {
        panic!("Should indicate interaction_way in boosting_interaction()")
    }
}

fn extract_snvs_interaction_not_in_vec(
    snvs_interaction_add: Vec<(usize, usize)>, // ~ 10k
    snvs_interaction: &[(usize, usize)],       // could be large; ~ 1m
) -> Vec<(usize, usize)> {
    //let snvs_interaction_add_reverse = snvs_interaction_add
    //    .iter()
    //    .map(|(i, j)| (*j, *i))
    //    .collect::<Vec<(usize, usize)>>();

    let snvs_interaction_common = snvs_interaction
        .iter()
        .filter(|(i, j)| {
            snvs_interaction_add.contains(&(*i, *j)) || snvs_interaction_add.contains(&(*j, *i))
        })
        .copied()
        .collect::<Vec<(usize, usize)>>();

    // intearction order follows snvs_interaction_add
    let snvs_interaction_uniq = snvs_interaction_add
        .iter()
        .filter(|(i, j)| {
            !snvs_interaction_common.contains(&(*i, *j))
                && !snvs_interaction_common.contains(&(*j, *i))
        })
        .copied()
        .collect::<Vec<(usize, usize)>>();

    //wrong
    //let snvs_interaction_uniq = snvs_interaction_common
    //    .iter()
    //    .map(|(i, j)| {
    //        if snvs_interaction_add.contains(&(*i, *j)) {
    //            (*i, *j)
    //        } else {
    //            (*j, *i)
    //        }
    //    })
    //    .collect::<Vec<(usize, usize)>>();

    snvs_interaction_uniq
}

// include duplicated interaction
//fn add_last_snvs_interaction(
//    wgts: &WgtBoosts,
//    boost_param: &BoostParam,
//    snvs: &Snvs,
//    genot: &Genot,
//    snvs_interaction: Vec<(usize, usize)>,
//) -> Vec<(usize, usize)> {
//    if let Some(interaction_way) = boost_param.interaction_way() {
//        let last_snvs_interaction = match interaction_way {
//            InteractionWay::SelectedBySelected => last_interaction_selected_by_selected(wgts),
//            InteractionWay::SelectedBySelectedLd => {
//                let last_snvs_interaction = last_interaction_selected_by_selected(wgts);
//
//                let last_snvs_interaction = filter_snvs_interaction(
//                    &last_snvs_interaction,
//                    boost_param.ld_criteria().unwrap(),
//                    //boost_param.ld_radius().unwrap(),
//                    snvs,
//                    genot,
//                );
//                last_snvs_interaction
//            }
//        };
//
//        // add
//        let mut snvs_interaction = snvs_interaction;
//        snvs_interaction.extend(last_snvs_interaction);
//        snvs_interaction
//    } else {
//        panic!("Should indicate interaction_way in boosting_interaction()")
//    }
//}
//
//// interaction excludes self by self
//fn create_snvs_interaction(
//    wgts: &WgtBoosts,
//    boost_param: &BoostParam,
//    snvs: &Snvs,
//    snvs_interaction: Vec<(usize, usize)>,
//) -> Vec<(usize, usize)> {
//    if let Some(interaction_way) = boost_param.interaction_way() {
//        match interaction_way {
//            //InteractionWay::OneBySelected => create_one_by_selected(wgts),
//            //InteractionWay::OneBySelectedLd => {
//            //    let snvs_interaction = create_one_by_selected(wgts);
//            //    // filter by ld
//            //    let snvs_interaction = filter_snvs_interaction(
//            //        &snvs_interaction,
//            //        boost_param.ld_radius().unwrap(),
//            //        snvs,
//            //    );
//            //    snvs_interaction
//            //}
//            InteractionWay::SelectedBySelected => create_selected_by_selected(wgts),
//            InteractionWay::SelectedBySelectedLd => {
//                // TODO: store current snvs_interaction and add new pair from newly selected snv
//                let snvs_interaction = create_selected_by_selected(wgts);
//
//                let snvs_interaction = filter_snvs_interaction(
//                    &snvs_interaction,
//                    boost_param.ld_radius().unwrap(),
//                    snvs,
//                );
//                snvs_interaction
//            }
//        }
//    } else {
//        panic!("Should indicate interaction_way in boosting_interaction()")
//    }
//}

/// allow input file and stdio for test
/// return (nsnvs, acc)
pub fn boosting_batch_interaction<W: std::io::Write>(
    writer: &mut BufWriter<W>,
    boost_param: &BoostParam,
    dataset: &Dataset,
    dataset_val: Option<&Dataset>,
    is_resume: bool,
    is_write_loss: bool,
    is_initial_only: bool,
    dout: Option<&DoutParaFile>, // for is_resume or is_write_loss, or saving initial loss
    finitial_snvs: Option<&Path>, //dout: Option<&Path>, // for is_resume or is_write_loss
                                 //dloss: Option<&Path>, //is_write_loss: bool, // TODO: make different dir in boosting
                                 //is_monitor: bool,
                                 // now in boost_param
                                 //nsnvs_monitor: Option<Vec<usize>>,
                                 //nsnvs_monitor: Option<&[usize]>,
) -> Option<(usize, f64)> {
    if !boost_param.cov_way().unwrap().is_first() {
        panic!("Should indicate CovWay::First in boosting_batch()")
    }
    //if use_const_for_loss {
    //    panic!("Deprecated: use_const_for_loss");
    //}

    let start_time = Instant::now();

    let phe = dataset.samples().phe_unwrap();
    let snvs = dataset.snvs();

    let n = dataset.samples().samples_n();
    let mut scores = SampleScore::new(n);

    let n_val = if let Some(dataset_val) = dataset_val {
        dataset_val.samples().samples_n()
    } else {
        0
    };
    let mut scores_val = SampleScore::new(n_val);

    let mut sample_weight = SampleWeight::new(
        n,
        boost_param.boost_type(),
        boost_param.loss_func(),
        boost_param.sample_weight_clip(),
        boost_param.sample_weight_wls_clip(),
    );
    sample_weight.renew_sample_weight(&scores, phe);

    let mut acc_monitor: BTreeMap<usize, f64> = BTreeMap::new();

    // save sample score of previous iteration
    // for duplicated snvs; to get sample score with 5 snvs, only way is wait until 6-th snv and calculate score of previous snvs
    let mut scores_val_monitor;

    let mut wgts = WgtBoosts::new(boost_param.boost_type());

    let t_start;
    if is_resume && dout.unwrap().is_exist_wgt_nonzero() {
        // assume cov is done when resume
        log::info!("resume");

        panic!("Do not use resume since results will be changed depending on batch.")
    } else {
        wgt_boost::io::write_cols(writer, boost_param.boost_type());
        t_start = 0;
    }

    if boost_param.is_dom_rec() {
        unimplemented!()
    }

    log::info!("Start Boosting.");
    let mut ti = t_start;

    // first cov: not in while loop to save scores_val_cov.
    {
        let p = boosting_covs(
            &mut wgts,
            dataset,
            &mut scores,
            &mut sample_weight,
            ti,
            dataset_val,
            &mut scores_val,
        );
        wgts.write_wgt_writer(writer);
        ti += p;
    }

    let scores_val_cov = scores_val.clone_align();

    // for monitor acc
    let mut nsnv_prev = 0usize;

    // iteration in batch
    let mut bi = 0usize;

    let mut loss_batch_stop: Option<f64> = None;
    let mut snvs_batch = HashSet::new();

    let mut snvs_interaction_batch = vec![];
    //let mut snvs_interaction_all = vec![];

    let m = dataset.snvs().snvs_n();
    // m_inter~batch
    let mut loss = LossStruct::new_interaction(boost_param.boost_type(), m, 0);
    // for batch, interactions only
    // m_inter could large. start with m
    // TOFIX: more m_inter for initial
    let mut loss_all = LossStruct::new_interaction_capacity(boost_param.boost_type(), m, 0, m);

    let mut snvs_interaction_all = if boost_param.batch_interaction_way().unwrap().is_initial() {
        create_initial_interaction(
            &mut loss_all,
            dataset,
            &mut sample_weight,
            boost_param,
            dout,
            finitial_snvs,
        )
    } else {
        vec![]
    };
    loss_all.resize_interaction(snvs_interaction_all.len());

    if is_initial_only {
        // TEMPORARY STOP HERE to get .loss
        log::info!("Stop here");
        //panic!("Stop here.");
        std::process::exit(0);
    }

    while run_next_iteration(ti, boost_param.iteration(), &wgts) {
        if ti % 100 == 0 {
            log::info!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
        } else {
            log::debug!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
        }

        if renew_batch(
            bi,
            boost_param.batch_way(),
            loss_batch_stop,
            &wgts,
            boost_param.learning_rate(),
        ) {
            log::debug!("New batch");
            (snvs_batch, loss_batch_stop) = create_extract_snvs_batch(
                // you can use loss_all
                &mut loss,
                dataset,
                &scores,
                &mut sample_weight,
                boost_param,
                None,
            );
            log::debug!("Batch wgt size {}", snvs_batch.len());

            log::debug!("New batch for interaction");
            log::debug!("snvs_interaction_all: {:?}", snvs_interaction_all.len());
            (snvs_interaction_batch, _, snvs_interaction_all) =
                create_extract_snvs_interaction_batch(
                    // loss_all for snvs_interaction_all
                    &mut loss_all,
                    dataset,
                    &scores,
                    &mut sample_weight,
                    boost_param,
                    snvs_interaction_all,
                );
            loss_all.resize_interaction(snvs_interaction_all.len());
            log::debug!("Batch wgt size {}", snvs_batch.len());
            log::debug!(
                "snvs_interaction_all filtered: {:?}",
                snvs_interaction_all.len()
            );

            bi = 0;
        }

        log::debug!("Batch iteration: {}", bi);
        log::debug!("Unique counts of SNVs: {}", wgts.count_unique_snv());

        // score before adding the next wgt
        scores_val_monitor = scores_val.clone_align();

        snvs_interaction_all = add_last_snvs_interaction(
            &wgts,
            boost_param,
            snvs,
            dataset.genot(),
            snvs_interaction_all,
        );
        loss_all.resize_interaction(snvs_interaction_all.len());
        //let snvs_interaction = create_snvs_interaction(&wgts, boost_param, snvs);
        if snvs_interaction_all.len() < 30 {
            log::debug!("snvs_interaction_all: {:?}", snvs_interaction_all);
        } else {
            log::debug!("snvs_interaction_all: {:?}", snvs_interaction_all.len());
        }

        // filter snvs_interaction_batch
        // cannot be after add_last_snvs since loss is nan for newly added interaction
        snvs_interaction_batch =
            filter_interaction_batch(snvs_interaction_batch, &loss, boost_param);
        // I will forget to do this. -> there is assert, so no problem
        loss.resize_interaction(snvs_interaction_batch.len());
        log::debug!(
            "snvs_interaction_batch filtered: {:?}",
            snvs_interaction_batch.len()
        );

        // snv interaction in the same batch will be candidate for the next interaction
        snvs_interaction_batch = add_last_snvs_interaction(
            &wgts,
            boost_param,
            snvs,
            dataset.genot(),
            snvs_interaction_batch,
        );
        loss.resize_interaction(snvs_interaction_batch.len());
        if snvs_interaction_batch.len() < 30 {
            log::debug!("snvs_interaction_batch: {:?}", snvs_interaction_batch);
        } else {
            log::debug!("snvs_interaction_batch: {:?}", snvs_interaction_batch.len());
        }

        let wgt = boosting_iter_snv_interaction(
            ti,
            boost_param,
            &mut loss,
            dataset,
            &mut scores,
            &mut sample_weight,
            Some(&snvs_batch),
            &snvs_interaction_batch,
            //&snvs_interaction_all,
            dataset_val,
            &mut scores_val,
        );

        wgts.add_wgt(wgt.clone());

        //log::debug!("wgts {:?}", wgts);
        wgts.write_wgt_writer(writer);

        if is_write_loss {
            if (ti < 200) || (ti % 100 == 0) {
                log::debug!("Write loss file.");
                let mut writer_loss = dout.unwrap().bufwriter_floss_adjmax(ti);
                loss.write_writer(&mut writer_loss, snvs);
            }
        }

        // count batch only for single snv
        if wgts.last_wgt().unwrap().wgt().is_snv_single() {
            bi += 1;
        }
        //bi += 1;
        ti += 1;

        if boost_param.is_monitor() {
            let dataset_val = dataset_val.unwrap();
            let nsnv = wgts.count_unique_snv();

            // To deal with duplicated snvs, monitor accuracy with the previous snv
            // only monitor acc when nsnv is changed
            // ex. To monitor acc of #snv=5, wait until nsnv=6 and calculate acc of previous iteration
            // This is when nsnv(=6) is changed from nsnv_prev(=5) and nsnv_prev(=5) is in nsnvs_monitor
            //
            // score_val_monitor is score of the previous iteration

            //if (nsnv != nsnv_prev) && is_update_acc(nsnv, nsnvs_monitor) {
            //log::debug!("nsnv, nsnv_prev: {}, {}", nsnv, nsnv_prev);
            if (nsnv != nsnv_prev) && is_update_acc(nsnv_prev, boost_param.monitor_nsnvs()) {
                log::debug!("monitor stat: nsnv_prev {}", nsnv_prev);
                update_acc_monitor(
                    &mut acc_monitor,
                    //nsnv,
                    nsnv_prev,
                    dataset_val,
                    &scores_val_monitor,
                    &scores_val_cov,
                    boost_param.acc_metric(),
                );

                // write acc
                let mut writer_acc = dout.unwrap().bufwriter_facc();
                write_acc(&acc_monitor, &mut writer_acc);

                let is_exit = monitor_acc(&acc_monitor, boost_param.monitor_nsnvs(), nsnv_prev);
                if is_exit {
                    log::info!("Exit since monitoring accuracy exceeded the threshold.");
                    let (nsnvs_max, acc_max) = compute_acc_max(&mut acc_monitor);
                    return Some((nsnvs_max, acc_max));
                }
            }

            nsnv_prev = nsnv;
        };

        log::debug!("\n");
    }

    return None;
}

fn write_acc<W: std::io::Write>(acc_monitor: &BTreeMap<usize, f64>, writer: &mut BufWriter<W>) {
    let str_header = "nsnv\tacc".to_owned();
    let strings = acc_monitor
        .iter()
        .map(|(nsnv, acc)| format!("{}\t{}", nsnv, acc))
        .collect::<Vec<String>>();

    let str = str_header + "\n" + &strings.join("\n");

    //log::debug!("str wgt {}", &str);
    writer.write(str.as_bytes()).unwrap();
    // capture error
    writer.flush().unwrap();
}

fn is_update_acc(nsnv: usize, nsnvs_monitor: Option<&[usize]>) -> bool {
    match nsnvs_monitor {
        Some(nsnvs_monitor) => nsnvs_monitor.contains(&nsnv),
        None => {
            if nsnv < 1000 {
                let nsnvs_monitor_first = [
                    5usize, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200,
                    300, 400, 500, 600, 700, 800, 900,
                ];
                nsnvs_monitor_first.contains(&nsnv)
            } else {
                // every 200
                nsnv % 200 == 0
                // every 1k
                //nsnv % 1000 == 0
            }
        }
    }
}

fn update_acc_monitor(
    acc_monitor: &mut BTreeMap<usize, f64>,
    //acc_monitor: &mut [f64],
    //nsnvs_monitor: &[usize],
    nsnv: usize,
    dataset_val: &Dataset,
    scores_val: &SampleScore,
    //scores_val: &[f64],
    scores_val_cov: &SampleScore,
    //scores_val_cov: &[f64],
    acc_metric: Option<AccMetric>,
) {
    acc_metric.expect("acc_metric should be indicated.");

    //let n = dataset_val.samples().samples_n();
    let acc = match acc_metric.unwrap() {
        AccMetric::CovAdjustedPseudoR2 => pgs::nagelkerke_r2(
            &dataset_val.samples().phe_unwrap().inner_bool(),
            //&dataset_val.samples().phe_unwrap().inner_i32(),
            //&dataset_val.samples().phe_unwrap().inner_f64(),
            scores_val.scores(),
            //&scores_val[..n],
            scores_val_cov.scores(),
            //&scores_val_cov[..n],
        ),
        AccMetric::AUC => unimplemented!(),
    };

    if acc.is_nan() {
        panic!("Accuracy is nan.");
    }

    //log::debug!("nsnvs_monitor {:?}", nsnvs_monitor);
    log::debug!("nsnv {:?}", nsnv);
    //let nsnv_index = nsnvs_monitor.iter().position(|x| *x == nsnv).unwrap();
    //log::debug!("nsnv_index {:?}", nsnv_index);
    log::debug!("acc {:?}", acc);
    // overwrite  if nsnv is duplicated
    acc_monitor.insert(nsnv, acc);
    //acc_monitor[nsnv_index] = acc;
    log::debug!("acc_monitor {:?}", acc_monitor);
}

//fn compute_acc_max(acc_monitor: &mut [f64], nsnvs_monitor: &[usize]) -> (usize, f64) {
//    let index_max = acc_monitor
//        .iter()
//        .enumerate()
//        .filter(|(_, a)| !a.is_nan())
//        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
//        .map(|(index, _)| index)
//        .unwrap();
//
//    let nsnv_max = nsnvs_monitor[index_max];
//    let acc_max = acc_monitor[index_max];
//    (nsnv_max, acc_max)
//}

fn compute_acc_max(acc_monitor: &BTreeMap<usize, f64>) -> (usize, f64) {
    acc_monitor
        .iter()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(nsnv, acc)| (*nsnv, *acc))
        .unwrap()

    //.filter(|(_, a)| !a.is_nan())
}

//fn monitor_acc(acc_monitor: &mut [f64], nsnvs_monitor: &[usize], nsnv: usize) -> bool {
fn monitor_acc(
    acc_monitor: &BTreeMap<usize, f64>,
    nsnvs_monitor: Option<&[usize]>,
    nsnv: usize,
) -> bool {
    log::debug!("nsnvs_monitor {:?}", nsnvs_monitor);
    log::debug!("acc_monitor {:?}", acc_monitor);

    let (nsnv_max, _) = compute_acc_max(acc_monitor);
    //let (nsnv_max, _) = compute_acc_max(acc_monitor, nsnvs_monitor);

    /*     let index_max = acc_monitor
        .iter()
        .enumerate()
        .filter(|(_, a)| !a.is_nan())
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(index, _)| index)
        .unwrap();

    let nsnv_max = nsnvs_monitor[index_max]; */
    //log::debug!("nsnvs_max {}", nsnv_max);
    //log::debug!("nsnvs_current {}", nsnv);
    log::debug!("nsnvs_current: {} vs max: {}", nsnv, nsnv_max);

    if let Some(nsnvs_monitor) = nsnvs_monitor {
        if nsnv >= *nsnvs_monitor.last().unwrap() {
            // confirm that last nsnv is in acc_monitor
            assert!(acc_monitor.contains_key(nsnvs_monitor.last().unwrap()));
            log::debug!("To exit since exceeded last of # snvs in monitoring nsnvs.");
            return true;
        }
    }
    //if !acc_monitor.last().unwrap().is_nan() {
    //    log::debug!("Exceeded last of # snvs in monitoring nsnvs.");
    //    return true;
    //}

    // if nsnv <= 500, do not stop
    // if 500 < nsnv <= 2k, stop if nsnv >= 2*nsnv_max
    // if 2k < nsnv, stop if nsnv >= nsnv_max + 1k
    if nsnv <= 500 {
        return false;
    } else if nsnv <= 2000 {
        if nsnv < 2 * nsnv_max {
            return false;
        } else {
            log::debug!("To exit since 500<nsnv<=2000 and nsnv>=2*max.");
            return true;
        }
    } else {
        if nsnv < nsnv_max + 1000 {
            return false;
        } else {
            log::debug!("To exit since nsnv>2000 and nsnv>=max+1000.");
            return true;
        }
    }

    // old: for monitor every 1k
    // if nsnv <= 500, do not stop
    // if 500 < nsnv <= 4k, stop if nsnv >= 4*nsnv_max
    // if 4k < nsnv, stop if nsnv >= nsnv_max + 3k
    //if nsnv <= 500 {
    //    return false;
    //} else if nsnv <= 4000 {
    //    if nsnv < 4 * nsnv_max {
    //        return false;
    //    } else {
    //        log::debug!("To exit since 500<nsnv<=4000 and nsnv>=4*max.");
    //        return true;
    //    }
    //} else {
    //    if nsnv < nsnv_max + 3000 {
    //        return false;
    //    } else {
    //        log::debug!("To exit since nsnv>4000 and nsnv>=max+3000.");
    //        return true;
    //    }
    //}

    //if nsnvs > nsnvs_max+2000 {
    //log::debug!("nsnvs_current: {} vs max: {}", nsnv as f64, nsnv_max as f64);
    //if (nsnv as f64) > (nsnv_max as f64) * 1.5f64 {
    //    return true;
    //}

    //false
}

fn set_wgt_m<W: WgtTrait>(wgts: &mut [W], snvs: &Snvs) {
    //let sida_to_m = wgt_sida_to_dataset_m(wgts, snvs);

    // here use &str not String
    let mut sida_to_m = HashMap::with_capacity(snvs.snvs_n());
    for (mi, snvid) in snvs.snv_ids().iter().enumerate() {
        sida_to_m.insert(snvid.sida(), mi);
    }

    for wgt in wgts.iter_mut() {
        //if wgt.kind().is_snv() {
        match wgt.kind() {
            //WgtKind::Snv(snv_id, ..) => {
            WgtKind::Snv(snv_wgt) => {
                let snv_id = snv_wgt.snv_id();
                let sida = snv_id.sida();
                //let sida = wgt.kind().snv_index().sida();
                //log::debug!("sida {}", sida);
                //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
                if let Some(mi) = sida_to_m.get(sida) {
                    wgt.set_snv_index_check(Some(*mi));
                } else {
                    panic!("Unknown SNV in wgt.");
                }
            }
            //WgtKind::SnvInteraction(snv_id_1, _, snv_id_2, _) => {
            WgtKind::SnvInteraction(snv_inter_wgt) => {
                let (snv_id_1, snv_id_2) = snv_inter_wgt.snv_ids();
                let sida_1 = snv_id_1.sida();
                let sida_2 = snv_id_2.sida();
                if let (Some(mi_1), Some(mi_2)) = (sida_to_m.get(sida_1), sida_to_m.get(sida_2)) {
                    wgt.set_snv_index_interaction_check(Some(*mi_1), Some(*mi_2));
                    //wgt.set_snv_index(Some(*mi));
                } else {
                    panic!("Unknown SNV in wgt.");
                }
            }
            WgtKind::Cov(_) => {}
        }

        // do nth for cov
        // not found
        //wgt_to_m.insert(wgti, None);
    }
}

/*
pub fn boosting_pruning(
    fout: &str,
    iter: usize,
    boost_type: &BoostType,
    predictions: &GenotBi<Vec<u8>>,
    //predictions: &[B8],
    ys: &[B8],
    m: usize,
    n: usize,
    snvs: &[Snv],
    //covs: &[Var],
    wgtcovs: &[WgtBoost],
    //wgtcovs: &[Wgt],
    is_dom_rec: bool,
    clump_r2: f64,
) {
    // TODO: this is not optimal. Is there any way to avoid create Vec?
    // use &[&SnvIndex] to input ?
    // -> solved
    //let snv_indexs: Vec<SnvIndex> = snvs.iter().map(|v| v.snv_index().clone()).collect();

    let mut score_ti: Vec<f64> = Vec::with_capacity(n);
    vec::push(&mut score_ti, 0.0, n);
    let mut ws_ti: Vec<f64> = Vec::with_capacity(n);
    vec::push(&mut ws_ti, 1.0 / (n as f64), n);

    let mut ps_ti: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
    vec::push(&mut ps_ti, 0.0, n + 32);

    let len_n = genot_bi::len_n(n);
    let mut pred: Vec<B8> = Vec::with_capacity(len_n);
    vec::push(&mut pred, 0x00, len_n);

    let mut losss_gt: Vec<f64> = Vec::with_capacity(2 * m);
    vec::push(&mut losss_gt, f64::NAN, 2 * m);

    let mut is_pruned: Vec<bool> = Vec::with_capacity(m); // or m+32
    vec::push(&mut is_pruned, false, m);

    let mut r2s: Vec<f64> = Vec::with_capacity(m);
    vec::push(&mut r2s, 0.0, m);

    //calculate_mean()
    // calculate_var
    let genotype_means = correlation::calculate_means_genotype(predictions, m, n);
    let genotype_vars = correlation::calculate_vars_genotype(&genotype_means, predictions, m, n);

    sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);

    // init fwgt
    io_wgt_boost::init_fwgt(fout);
    io_wgt_boost::init_fws(fout);

    // load wgtcovs
    let p = wgtcovs.len();

    for ti in 0..p {
        let pi = ti;
        sample_weight::renew_score_var(&mut score_ti, &wgtcovs[pi]);
        wgtcovs[pi].write_wgt(fout);
    }

    log::debug!("after cov");
    sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
    //sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
    sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
    //log::debug!("score {:?}", score_ti);
    //log::debug!("ws {:?}", ws_ti);
    //log::debug!("ps {:?}", ps_ti);

    // for constonly
    let ti = p;
    let wgt = loss::create_loss_const(ti, &ps_ti, &mut pred, &ys, n);
    wgt.write_wgt(fout);
    compute::compute_pred(&mut pred, &wgt, predictions, n);

    sample_weight::renew_score(&mut score_ti, &pred, &wgt);
    sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
    //sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
    sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
    //log::debug!("score {:?}", score_ti);
    //log::debug!("ws {:?}", ws_ti);
    //log::debug!("ps {:?}", ps_ti);

    //print max,min of ws
    log::debug!(
        "ps max,min {},{}",
        ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.max(v)),
        ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.min(v)),
    );
    //let max = iter.fold(0.0 / 0.0, f64::max); //somehow error

    // calculate loss beforehand just once
    let mut losss_gt: Vec<f64> = Vec::with_capacity(2 * m);
    vec::push(&mut losss_gt, f64::NAN, 2 * m);
    loss::calculate_loss_gt(&mut losss_gt, predictions, &ps_ti, ys, n);

    for ti in p + 1..iter {
        log::debug!("\nIteration: {}", ti);

        //log::debug!("is_pruned {:?}", is_pruned);

        //let wgt_gt = loss::search_min_loss_gt_pruned(
        let mut wgt = loss::search_min_loss_gt_pruned(
            &losss_gt,
            ti,
            predictions,
            &ps_ti,
            ys,
            n,
            snvs,
            &is_pruned,
        );

        let abcd_sum = wgt.abcd_sum().unwrap();
        log::debug!("abcd {:?}", abcd_sum);

        let coef_ti = coefficient::calculate_coefficients(&wgt, boost_type);
        log::debug!("alpha, const {:?}", coef_ti);
        //log::debug!("const, alpha: {:.4}, {:.4}", coef_ti.0, coef_ti.1);
        wgt.wgt_mut().model_mut().set_coef_binary(coef_ti);
        wgt.write_wgt(fout);

        compute::compute_pred(&mut pred, &wgt, predictions, n);

        sample_weight::renew_score(&mut score_ti, &pred, &wgt);
        sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
        //sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
        sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
        //log::debug!("score {:?}", score_ti);
        //log::debug!("ws {:?}", ws_ti);
        //log::debug!("ps {:?}", ps_ti);

        log::debug!(
            "ps max,min {},{}",
            ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.max(v)),
            ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.min(v)),
        );

        io_wgt_boost::write_ws(fout, &ws_ti, ti);

        if is_dom_rec {
            log::debug!("Choose the other inheritance model.");

            let mut wgt_second =
                loss::search_wgt_gt_second(&wgt, &losss_gt, ti, predictions, &ps_ti, ys, n, snvs);

            let abcd_sum = wgt_second.abcd_sum();
            log::debug!("abcd {:?}", abcd_sum);

            let coef_ti = coefficient::calculate_coefficients(&wgt_second, boost_type);
            log::debug!("alpha, const {:?}", coef_ti);
            //log::debug!("const, alpha: {:.4}, {:.4}", coef_ti.0, coef_ti.1);
            wgt_second.wgt_mut().model_mut().set_coef_binary(coef_ti);
            wgt_second.write_wgt(fout);

            compute::compute_pred(&mut pred, &wgt_second, predictions, n);

            sample_weight::renew_score(&mut score_ti, &pred, &wgt_second);
            sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
            sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
            //log::debug!("score {:?}", score_ti);
            //log::debug!("ws {:?}", ws_ti);
            //log::debug!("ps {:?}", ps_ti);

            log::debug!(
                "ps max,min {},{}",
                ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.max(v)),
                ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.min(v)),
            );

            io_wgt_boost::write_ws(fout, &ws_ti, ti);
        }

        // calculate r2 of chosen SNVs and all SNVs
        correlation::calculate_r2s(
            &mut r2s,
            predictions,
            &genotype_means,
            &genotype_vars,
            &is_pruned,
            &wgt,
            &snvs,
            //&snv_indexs,
            n,
        );

        // update pruned
        correlation::update_is_pruned(
            &mut is_pruned,
            &r2s,
            clump_r2,
            wgt.wgt().kind().index_snv().unwrap(),
        );

        if correlation::check_all_pruned(&is_pruned) {
            log::debug!("All SNVs are pruned at {}-th iteration.", ti);
            break;
        }
    }
}
 */

/*
pub fn boosting_ss(
    fout: &str,
    iter: usize,
    predictions: &GenotTwin,
    //predictions: &[B8],
    m: usize,
    n: usize,
    sum_stats: &[SumStat],
    //covs: &[Var],
    wgtcovs: &[WgtBoost],
    //wgtcovs: &[Wgt],
    is_dom_rec: bool,
    clump_r2: f64,
) {
    let mut is_pruned: Vec<bool> = Vec::with_capacity(m); // or m+32
    vec::push(&mut is_pruned, false, m);

    let mut r2s: Vec<f64> = Vec::with_capacity(m);
    vec::push(&mut r2s, 0.0, m);

    //calculate_mean()
    // calculate_var
    let genotype_means = correlation::calculate_means_genotype(predictions, m, n);
    let genotype_vars = correlation::calculate_vars_genotype(&genotype_means, predictions, m, n);

    // init fwgt
    io_wgt_boost::init_fwgt(fout);
    io_wgt_boost::init_fws(fout);

    // load wgtcovs
    let p = wgtcovs.len();

    for ti in 0..p {
        let pi = ti;
        wgtcovs[pi].write_wgt(fout);
    }

    log::debug!("after cov");

    //log::debug!("snvs {:?}", snvs);

    //let snv_indexs: Vec<SnvIndex> = sum_stats.iter().map(|v| v.snv_index().clone()).collect();

    for ti in p + 1..iter {
        log::debug!("\nIteration: {}", ti);

        //log::debug!("is_pruned {:?}", is_pruned);

        //let wgt_gt = loss::search_min_loss_gt_pruned(
        let wgt = loss::search_min_loss_gt_pruned_loss_ss(sum_stats, ti, &is_pruned);

        wgt.write_wgt(fout);

        if is_dom_rec {
            log::debug!("Choose the other inheritance model.");

            let wgt_second = loss::search_min_loss_gt_pruned_loss_ss_second(&wgt, sum_stats, ti);

            wgt_second.write_wgt(fout);
        }

        // calculate r2 of chosen SNVs and all SNVs
        correlation::calculate_r2s(
            &mut r2s,
            predictions,
            &genotype_means,
            &genotype_vars,
            &is_pruned,
            &wgt,
            &sum_stats,
            //&snv_indexs,
            n,
        );

        // update pruned
        correlation::update_is_pruned(
            &mut is_pruned,
            &r2s,
            clump_r2,
            wgt.wgt().kind().index_snv().unwrap(),
        );

        if correlation::check_all_pruned(&is_pruned) {
            log::debug!("All SNVs are pruned at {}-th iteration.", ti);
            break;
        }
    }
}
 */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_snvs_interaction_not_in_vec() {
        let snvs_interaciton_add: Vec<(usize, usize)> = vec![(1, 2), (10, 11), (6, 5), (13, 12)];

        let snvs_interaction_all = vec![(1, 2), (3, 4), (5, 6), (7, 8)];

        let snvs_interaction_uniq_exp = vec![(10, 11), (13, 12)];

        let snvs_interaction_uniq =
            extract_snvs_interaction_not_in_vec(snvs_interaciton_add, &snvs_interaction_all);

        assert_eq!(snvs_interaction_uniq, snvs_interaction_uniq_exp);
    }

    #[test]
    fn test_monitor_acc() {
        // nsnv=nsnvs_monitor last
        let acc_monitor: BTreeMap<usize, f64> =
            BTreeMap::from([(100, 0.5), (500, 0.2), (1000, 0.1)]);
        let nsnvs_monitor: Option<&[usize]> = Some(&[100, 500, 1000]);
        let nsnv = 1000;

        assert_eq!(monitor_acc(&acc_monitor, nsnvs_monitor, nsnv), true);
    }

    #[test]
    fn test_monitor_acc_2() {
        // nsnv<=500
        let acc_monitor: BTreeMap<usize, f64> = BTreeMap::from([(100, 0.5)]);
        let nsnv = 100;

        assert_eq!(monitor_acc(&acc_monitor, None, nsnv), false);
    }

    #[test]
    fn test_monitor_acc_3() {
        // nsnv<=2k and nsnv<2*nsnv_max
        // max nsnvs=1000
        let acc_monitor: BTreeMap<usize, f64> = BTreeMap::from([(1500, 0.5), (2000, 0.1)]);
        let nsnv = 2000;

        assert_eq!(monitor_acc(&acc_monitor, None, nsnv), false);
    }

    #[test]
    fn test_monitor_acc_4() {
        // nsnv<=2k and nsnv>=2*nsnv_max
        // max nsnvs=100
        let acc_monitor: BTreeMap<usize, f64> =
            BTreeMap::from([(100, 0.5), (500, 0.2), (1000, 0.1), (4000, 0.05)]);
        let nsnv = 2000;

        assert_eq!(monitor_acc(&acc_monitor, None, nsnv), true);
    }

    #[test]
    fn test_monitor_acc_5() {
        // nsnv>2k and nsnv<nsnv_max+1k
        // max nsnvs=3400
        let acc_monitor: BTreeMap<usize, f64> =
            BTreeMap::from([(1000, 0.1), (3400, 0.5), (4000, 0.05)]);
        let nsnv = 4000;

        assert_eq!(monitor_acc(&acc_monitor, None, nsnv), false);
    }

    #[test]
    fn test_monitor_acc_6() {
        // nsnv>2k and nsnv>=nsnv_max+1k
        // max nsnvs=3000
        let acc_monitor: BTreeMap<usize, f64> =
            BTreeMap::from([(1000, 0.1), (3000, 0.5), (4000, 0.05)]);
        let nsnv = 4000;

        assert_eq!(monitor_acc(&acc_monitor, None, nsnv), true);
    }
}
