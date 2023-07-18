//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
//! Loading genotype requires the long time. Extracting genotype takes negligible time.
//!
// split covs from samples

//mod correlation; // not implemented cmatrix yet
mod coefficient;
mod compute_pred;
mod epsilon;
mod iteration;
pub mod loss; // pub for bench
mod regression_cov;
pub mod sample_weight;
mod table;

use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufWriter;
use std::iter::FromIterator;
use std::path::Path;
use std::time::Instant;

use self::sample_weight::SampleWeight;

use super::WgtBoost;
use crate::boosting_param::{BatchWay, BoostParam, BoostType, CovWay, IterationNumber, LossFunc};
use crate::boosting_score::score;
use crate::wgt_boost;
use crate::wgt_boosts::WgtBoosts;
use genetics::genot::prelude::*;
use genetics::pgs;
use genetics::samples::prelude::*;
use genetics::wgt::{Coef, WgtKind};
use genetics::Model;
use genetics::Wgt;
use genetics::WgtTrait;
use genetics::{samples::CovsTrait, Covs, Dataset, Genot, Snvs};
use loss::LossStruct;
pub use table::ContingencyTable;

fn boosting_iter_cov(score_ti: &mut [f64], wgtcov: &WgtBoost, covs: &Covs) {
    let cov_name = wgtcov.wgt().kind().cov().name();
    let cov_vals_v;
    let cov_vals: &[f64];
    // TODO: cleaner; put these in Covs or here?
    if cov_name == "const" {
        cov_vals_v = vec![1.0; covs.vals().unwrap()[0].len()];
        cov_vals = &cov_vals_v;
    } else {
        cov_vals = covs.vals_id(cov_name);
    }
    //match wgt.wgt().model().coef() {

    sample_weight::renew_score_var(score_ti, wgtcov.wgt().model().coef(), cov_vals);
    //wgtcov.clone()
}

fn boosting_iter_snv(
    ti: usize,
    boost_param: BoostParam,
    loss: &mut LossStruct,
    genot: &Genot,
    pred: &mut [u8],
    scores: &mut [f64],
    sample_weight: &mut SampleWeight,
    phe: &Phe,
    snvs: &Snvs,
    skip_snv: &HashSet<usize>,
    use_adjloss: bool,
    dataset_val: Option<&Dataset>,
    scores_val: &mut [f64],
    pred_val: &mut [u8],
) -> WgtBoost {
    sample_weight.print_stat(phe);

    let mut wgt = loss::search_min_loss_gt(
        loss,
        ti,
        genot,
        &sample_weight,
        phe,
        snvs,
        boost_param,
        skip_snv,
        use_adjloss,
    );

    // why do I need pred??
    // now already split table and loss so should be not necessary
    // -> for renew_score() now
    compute_pred::compute_pred(pred, &wgt, genot);

    // set coef
    if boost_param.boost_type().is_type_ada() {
        //let ps_pad = sample_weight.ps_pad().unwrap();
        /*
        // create table
        let (table_sum, is_eps) =
            table::calculate_table_eps(&pred, ps, phe, boost_param.eps(), boost_param.boost_type());
        log::debug!("table {:?}", table_sum);
        wgt.set_contingency_table(table_sum, is_eps);

        // TODO: just end without error when loss==1.0

        let coef_ti = coefficient::calculate_coefficients(
            wgt.contingency_table().unwrap(),
            boost_param.boost_type(),
            boost_param.learning_rate(),
            boost_param.eps(),
        );
        log::debug!("coef {:?}", coef_ti);
        */

        // not set cont table
        //wgt.set_contingency_table(table_sum, is_eps);
        let (coef_ti, is_eps) = coefficient::calculate_coef_root_ada(
            &pred,
            sample_weight,
            //ps_pad,
            phe,
            boost_param.learning_rate(),
            boost_param.eps(),
            boost_param.boost_type(),
        );
        log::debug!("coef {:?}", coef_ti);
        wgt.set_coef(coef_ti);
        wgt.set_is_eps(is_eps);
        log::debug!("wgt {:?}", wgt);
    } else if boost_param.boost_type().is_type_logit() {
        let mi: usize = match wgt.wgt().kind() {
            WgtKind::Snv(_, _, mi) => mi.unwrap(),
            _ => panic!(),
        };

        //let coef_ti: Coef;
        // TODO: cleaner
        let (coef_ti, is_eps, is_eff_eps) = coefficient::calculate_coef_root_logit(
            &genot.to_genot_snv(mi),
            sample_weight,
            //&pred,
            //wzs_pad,
            //wls_pad,
            phe,
            boost_param.learning_rate(),
            boost_param.eps(),
            boost_param.eff_eps(),
            boost_param.boost_type(),
        );

        log::debug!("coef {:?}", coef_ti);
        wgt.set_coef(coef_ti);
        wgt.set_is_eps(is_eps);
        wgt.set_is_eff_eps(is_eff_eps);

        log::debug!("wgt {:?}", wgt);
    };

    /*
    // TMP
    // TODO: --verbose --verbose
    if count_ps_snv2(pred) < 30 {
        //log::debug!("snv2");
        print_snv2_ps(ps, phe, pred, "ps");
        //log::debug!("prob");
        print_snv2_ps(probs, phe, pred, "prob");
        //log::debug!("zs");
        print_snv2_ps(zs, phe, pred, "zs");
        //log::debug!("wzs");
        print_snv2_ps(wzs, phe, pred, "wzs");
        //log::debug!("wls");
        print_snv2_ps(wls, phe, pred, "wls");
    }
     */

    sample_weight::renew_score(scores, pred, &wgt);

    if let Some(dataset_val) = dataset_val {
        let genot_val = dataset_val.genot();
        compute_pred::compute_pred(pred_val, &wgt, genot_val);
        sample_weight::renew_score(scores_val, pred_val, &wgt);
    }

    sample_weight.renew_sample_weight(scores, phe);

    /*
    let n = phe.n();
    if n < 10 {
        log::debug!("score {:?}", scores);
        log::debug!("ws {:?}", ws);
        log::debug!("ps {:?}", ps);
    }
    */

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

    // TODO: create argument to output ws file --output ws
    //io_wgt_boost::write_ws(fout, &ws_ti, ti);

    wgt_second
}
*/

pub fn boosting_covs(
    wgts: &mut WgtBoosts,
    dataset: &Dataset,
    scores: &mut [f64],
    sample_weight: &mut SampleWeight,
    iteration_start: usize,
    dataset_val: Option<&Dataset>,
    scores_val: &mut [f64],
) -> usize {
    // no weighted only
    assert_eq!(iteration_start, 0);
    //let wgtcovs_logreg = regression_cov::logistic_regression_covs(dataset.samples(), iteration_start);
    //let wgtcovs_logreg = regression_cov::logistic_regression_covs(dataset.samples(), iteration_start);
    // FIXME: for now, .is_type_logit(), ps can be constructed for the first iteration, so only able for CovWay::First
    let wgtcovs_logreg =
        regression_cov::logistic_regression_covs(dataset.samples(), iteration_start);
    //let wgtcovs_logreg = regression_cov::logistic_regression_covs_sampleweight(
    //    dataset.samples(),
    //    sample_weight.ps().unwrap(),
    //    //ps,
    //    iteration_start,
    //);
    //let wgtcovs_logreg = regression_cov::logistic_regression_covs(dataset.samples());

    log::debug!("wgtcovs_logreg: {:?}", &wgtcovs_logreg);

    let covs = dataset.samples().covs().unwrap();

    let p = wgtcovs_logreg.len();
    // TODO: why splitting...?
    for pi in 0..p {
        let wgtcov = wgtcovs_logreg[pi].clone();
        boosting_iter_cov(scores, &wgtcov, covs);
        wgts.add_wgt(wgtcov);
    }

    if let Some(dataset_val) = dataset_val {
        let covs = dataset_val.samples().covs().unwrap();

        let p = wgtcovs_logreg.len();
        for pi in 0..p {
            let wgtcov = wgtcovs_logreg[pi].clone();
            boosting_iter_cov(scores_val, &wgtcov, covs);
        }
    }

    log::debug!("after cov");
    let phe = dataset.samples().phe();
    sample_weight.renew_sample_weight(scores, phe);

    p
}

pub fn boosting_logit_const(
    dataset: &Dataset,
    //pred: &mut [u8],
    scores: &mut [f64],
    sample_weight: &mut SampleWeight,
    iteration: usize,
    dataset_val: Option<&Dataset>,
    scores_val: &mut [f64],
) -> WgtBoost {
    let wzs = sample_weight.wzs().unwrap();
    let wls = sample_weight.wls().unwrap();
    let (coef, _) = coefficient::calculate_coef_logit_const(&wzs, &wls);
    let model = Model::new_coef(coef);
    let wgt = Wgt::construct_const(model);
    let wgt_boost = WgtBoost::construct_wgt_iteration(wgt, iteration);

    log::debug!("after const");
    log::debug!("wgt {:?}", wgt_boost);

    //wgts.add_wgt(wgt_boost.clone());

    sample_weight::renew_score(scores, &[0u8], &wgt_boost);

    if let Some(_) = dataset_val {
        sample_weight::renew_score(scores_val, &[0u8], &wgt_boost);
    }

    let phe = dataset.samples().phe();

    sample_weight.renew_sample_weight(&scores, phe);

    wgt_boost
}

const ITERATION_NUMBER_SNV_LIMIT: usize = 500_000;

fn run_next_iteration(ti: usize, iteration: IterationNumber, wgts: &WgtBoosts) -> bool {
    // TODO: Snv() might not stop, so determine ti limit
    match iteration {
        IterationNumber::Iteration(x) => ti < x,
        // use '<=' since the last snv might be selected in a row. ex. rs1, rs1, rs2, rs2. if '< 2',  rs2 will be selected only once
        IterationNumber::Snv(x) => {
            if ti > ITERATION_NUMBER_SNV_LIMIT {
                log::info!(
                    "Iteration number exceeded the limit, so stop the iteration: {}",
                    ti
                );
                return false;
            }
            wgts.count_unique_snv() <= x
        }
    }
}

fn run_next_iteration_monitor(
    ti: usize,
    iteration: IterationNumber,
    wgts: &WgtBoosts,
    is_monitor: bool,
) -> bool {
    if is_monitor {
        if ti > ITERATION_NUMBER_SNV_LIMIT {
            log::info!(
                "Iteration number exceeded the limit, so stop the iteration: {}",
                ti
            );
            return false;
        }
        return true;
    }
    // TODO: Snv() might not stop, so determine ti limit
    match iteration {
        IterationNumber::Iteration(x) => ti < x,
        // use '<=' since the last snv might be selected in a row. ex. rs1, rs1, rs2, rs2. if '< 2',  rs2 will be selected only once
        IterationNumber::Snv(x) => {
            if ti > ITERATION_NUMBER_SNV_LIMIT {
                log::info!(
                    "Iteration number exceeded the limit, so stop the iteration: {}",
                    ti
                );
                return false;
            }
            wgts.count_unique_snv() <= x
        }
    }
}

// Why BatchWay is option??
fn run_next_iteration_in_batch(
    bi: usize,
    batch_way: Option<BatchWay>,
    loss_batch_stop: f64,
    wgts: &WgtBoosts,
    learning_rate: f64,
) -> bool {
    if batch_way.unwrap().use_comp_loss() {
        if bi > 0 {
            let loss_iter = wgts.last_wgt().unwrap().loss_option();

            if let Some(loss_iter) = loss_iter {
                if loss_iter > loss_batch_stop {
                    log::debug!(
                        "last loss is larger than loss_batch_stop {} > {}",
                        loss_iter,
                        loss_batch_stop
                    );
                    return false;
                }
            }
            // otherwise, last wgt was cov
        }
    }

    let batch_max_iter = batch_way.unwrap().batch_max_iter(learning_rate);

    if bi >= batch_max_iter {
        log::debug!("bi reached batch_max_iter");
    }

    bi < batch_max_iter
}

fn run_cov_iteration(wgts: &WgtBoosts, cov_way: Option<CovWay>) -> bool {
    match cov_way {
        None => false,
        Some(cov_way) => match cov_way {
            CovWay::First => wgts.wgts().len() == 0,
            CovWay::Every(_) => {
                // TODO:implement to allow cov_way=first only
                panic!("Cannot run weighted regression.");
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

fn eq_float(f1: f64, f2: f64) -> bool {
    (f1 - f2).abs() < 1e-7
}

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

/// allow input file and stdio for test
pub fn boosting<W: std::io::Write>(
    writer: &mut BufWriter<W>,
    boost_param: BoostParam,
    dataset: &Dataset,
    //writer_loss: &mut BufWriter<W>,
    //fout: Option<&str>, // should Some() when is_resume=T or is_write_loss=T
    use_adjloss: bool,
    is_resume: bool,
    is_write_loss: bool,
    dout: Option<&Path>, // for is_resume or is_write_loss
                         //dloss: Option<&Path>, //is_write_loss: bool, // TODO: make different dir in boosting
) -> Option<(usize, f64)> {
    let start_time = Instant::now();

    let phe = dataset.samples().phe();
    let genot = dataset.genot();
    let snvs = dataset.snvs();

    let n = dataset.samples().samples_n();
    let mut scores: Vec<f64> = vec![0.0; n];

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
    let mut pred: Vec<u8> = vec![0u8; n];

    let m = dataset.snvs().snvs_n();
    let mut loss = LossStruct::new(boost_param.boost_type(), m);

    let mut wgts = WgtBoosts::new(boost_param.boost_type());

    let t_start;
    //let is_exist_wgt = wgt_boost::io::is_exist_wgt(fout);
    if is_resume & wgt_boost::io::is_exist_wgt_nonzero(dout.unwrap()) {
        //if is_exist_wgt & is_resume {
        // assume cov is done when resume
        log::info!("resume");

        wgts = WgtBoosts::new_from_file_dir(dout.unwrap(), boost_param.boost_type());
        // align alleles
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
    } else {
        wgt_boost::io::write_cols(writer, boost_param.boost_type());
        t_start = 0;
    }

    if boost_param.is_dom_rec() {
        unimplemented!()
    }

    log::info!("Start Boosting.");
    let mut ti = t_start;

    let iteration = boost_param.iteration();

    // index of wgt
    let mut wgt_skip: HashSet<usize> = HashSet::new();
    let mut wgt_last: Option<WgtBoost> = None;
    while run_next_iteration(ti, iteration, &wgts) {
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
                &mut [0.0f64; 0],
            );

            wgts.write_wgt_writer(writer);

            ti += p;
        } else {
            log::debug!("Run snv");
            let wgt = boosting_iter_snv(
                ti,
                boost_param,
                &mut loss,
                genot,
                &mut pred,
                &mut scores,
                &mut sample_weight,
                phe,
                snvs,
                &wgt_skip,
                use_adjloss,
                None,
                &mut [0.0f64; 0],
                &mut [0u8; 0],
            );

            let skip_this_snv = renew_wgt_skip(&mut wgt_skip, wgt_last.as_ref(), &wgt);
            // is this ok?
            if skip_this_snv {
                continue;
            }

            wgts.add_wgt(wgt.clone());

            //log::debug!("wgts {:?}", wgts);
            wgts.write_wgt_writer(writer);

            wgt_last = Some(wgt);

            ti += 1;
        }

        if is_write_loss {
            if (ti < 200) || (ti % 100 == 0) {
                let mut writer_loss = wgt_boost::io::bufwriter_floss(dout.unwrap(), ti);
                loss.write_writer(&mut writer_loss, snvs);
            }
        }

        log::debug!("\n");
    }

    None
}

fn create_wgt_skip_snvs(
    wgt_skip_batch: &HashSet<usize>,
    wgt_skip: &HashSet<usize>,
) -> HashSet<usize> {
    let mut wgt_skip_snvs = wgt_skip_batch.clone();
    wgt_skip_snvs.extend(wgt_skip);
    //let wgt_skip_snvs: HashSet<usize> = wgt_skip_batch.union(wgt_skip).collect();
    wgt_skip_snvs
}

fn create_wgt_skip_batch(
    loss: &mut LossStruct,
    genot: &Genot,
    sample_weight: &SampleWeight,
    phe: &Phe,
    boost_param: BoostParam,
    skip_snv: &HashSet<usize>,
    use_adjloss: bool,
) -> (HashSet<usize>, f64) {
    loss::calculate_loss_gt(
        loss,
        genot,
        &sample_weight,
        phe,
        boost_param,
        skip_snv,
        use_adjloss,
    );

    let m_top = boost_param.batch_way().unwrap().batch_size();
    log::debug!("batch m {}", m_top);

    let (use_snvs, loss_next) = loss.search_topprop_modelfree_n(m_top, Some(skip_snv));
    log::debug!("loss batch next {}", loss_next);

    // for print
    let (loss_min, _, _) = loss.search_min(skip_snv);
    log::debug!("loss min {}", loss_min);

    let wgt_use_batch: HashSet<usize> = HashSet::from_iter(
        use_snvs
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(i, _)| i),
    );

    let m = loss.inner().len();

    let wgt_skip_batch: HashSet<usize> = HashSet::from_iter(0usize..m)
        .difference(&wgt_use_batch)
        .copied()
        .collect();

    (wgt_skip_batch, loss_next)
}

/// allow input file and stdio for test
/// return (nsnvs, acc)
pub fn boosting_batch<W: std::io::Write>(
    writer: &mut BufWriter<W>,
    boost_param: BoostParam,
    dataset: &Dataset,
    dataset_val: Option<&Dataset>,
    use_adjloss: bool,
    use_const_for_loss: bool,
    is_resume: bool,
    is_write_loss: bool,
    dout: Option<&Path>, // for is_resume or is_write_loss
    //dloss: Option<&Path>, //is_write_loss: bool, // TODO: make different dir in boosting
    is_monitor: bool,
    nsnvs_monitor: Option<Vec<usize>>,
    //nsnvs_monitor: Option<&[usize]>,
) -> Option<(usize, f64)> {
    if !boost_param.cov_way().unwrap().is_first() {
        panic!("Should indicate CovWay::First in boosting_batch()")
    }
    if use_const_for_loss {
        panic!("Deprecated: use_const_for_loss");
    }

    let start_time = Instant::now();

    let phe = dataset.samples().phe();
    let genot = dataset.genot();
    let snvs = dataset.snvs();

    let n = dataset.samples().samples_n();
    let mut scores: Vec<f64> = vec![0.0; n];
    let n_val = if let Some(dataset_val) = dataset_val {
        dataset_val.samples().samples_n()
    } else {
        0
    };
    let mut scores_val: Vec<f64> = vec![0.0; n_val];

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
    let mut pred: Vec<u8> = vec![0u8; n];
    let mut pred_val: Vec<u8> = vec![0u8; n_val];

    let m = dataset.snvs().snvs_n();

    // TODO: make interval even
    let nsnvs_monitor = match nsnvs_monitor {
        Some(x) => x,
        None => {
            // TODO: make interval even
            let nsnvs_monitor = [
                5usize, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 300,
                400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                10000,
            ];
            nsnvs_monitor.to_vec()
        }
    };
    //let nsnvs_monitor = [
    //    5usize, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 300, 400,
    //    500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
    //];
    let mut acc_monitor = vec![f64::NAN; nsnvs_monitor.len()];
    // save sample score of previous iteration
    // for duplicated snvs; to get sample score with 5 snvs, only way is wait until 6-th snv and calculate score of previous snvs
    let mut scores_val_monitor: Vec<f64> = vec![0.0; n_val];

    let mut loss = LossStruct::new(boost_param.boost_type(), m);
    let mut wgts = WgtBoosts::new(boost_param.boost_type());

    let t_start;
    if is_resume & wgt_boost::io::is_exist_wgt_nonzero(dout.unwrap()) {
        // assume cov is done when resume
        log::info!("resume");

        panic!("Do not use resume since reesults will be changed.")

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

    let iteration = boost_param.iteration();

    // first cov
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

    let scores_val_cov = scores_val.clone();

    /*     // for const
       if use_const_for_loss {
           let wgt = boosting_logit_const(
               dataset,
               &mut scores,
               &mut sample_weight,
               ti,
               dataset_val,
               &mut scores_val,
           );
           wgts.add_wgt(wgt.clone());
           wgts.write_wgt_writer(writer);

           ti += 1;
       }
    */
    // TODO: cleaner for iteration and integrate boosting()
    // while b {
    //       if b {
    //                   create_wgt_skip_batch
    //        }
    //       if b {
    //                    iteration
    //        }
    // }
    // index of wgt
    //let mut wgt_skip_snvs: HashSet<usize> = HashSet::new();
    // for illegal snvs
    let mut wgt_skip: HashSet<usize> = HashSet::new();
    let mut wgt_last: Option<WgtBoost> = None;
    // while true {
    while run_next_iteration_monitor(ti, iteration, &wgts, is_monitor) {
        log::debug!("New batch");
        let (wgt_skip_batch, loss_batch_stop) = create_wgt_skip_batch(
            &mut loss,
            genot,
            &mut sample_weight,
            phe,
            boost_param,
            &wgt_skip,
            use_adjloss,
        );

        log::debug!("Batch wgt size {}", m - wgt_skip_batch.len());

        // iteration in batch
        let mut bi = 0usize;

        while run_next_iteration_in_batch(
            bi,
            boost_param.batch_way(),
            loss_batch_stop,
            &wgts,
            boost_param.learning_rate(),
        ) & run_next_iteration_monitor(ti, iteration, &wgts, is_monitor)
        {
            if ti % 100 == 0 {
                log::info!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
            } else {
                log::debug!("Iteration: {}, {} sec", ti, start_time.elapsed().as_secs());
            }

            log::debug!("Batch iteration: {}", bi);
            log::debug!("Unique counts of SNVs: {}", wgts.count_unique_snv());

            let wgt_skip_iter = create_wgt_skip_snvs(&wgt_skip_batch, &wgt_skip);
            log::debug!(
                "Batch wgt_batch size, illegal size {} {}",
                m - wgt_skip_batch.len(),
                wgt_skip.len()
            );
            log::debug!("Batch wgt_iter size {}", m - wgt_skip_iter.len());

            // here before w, z are renewed
            //let loss_max_theory = compute_loss_max_theory(&zs, &wls, phe.n());
            //log::debug!("loss max {}", loss_max_theory);

            // score before adding the next wgt
            scores_val_monitor.clone_from_slice(&scores_val);

            let wgt = boosting_iter_snv(
                ti,
                boost_param,
                &mut loss,
                genot,
                &mut pred,
                &mut scores,
                &mut sample_weight,
                phe,
                snvs,
                &wgt_skip_iter,
                use_adjloss,
                dataset_val,
                &mut scores_val,
                &mut pred_val,
            );

            let skip_this_snv = renew_wgt_skip(&mut wgt_skip, wgt_last.as_ref(), &wgt);
            // is this ok?
            if skip_this_snv {
                continue;
            }

            wgts.add_wgt(wgt.clone());

            //log::debug!("wgts {:?}", wgts);
            wgts.write_wgt_writer(writer);

            wgt_last = Some(wgt);

            if is_write_loss {
                if (ti < 200) || (ti % 100 == 0) {
                    log::debug!("Write loss file");

                    if use_adjloss {
                        let mut writer_loss =
                            wgt_boost::io::bufwriter_floss_adjmax(dout.unwrap(), ti);
                        loss.write_writer(&mut writer_loss, snvs);
                    } else {
                        let mut writer_loss = wgt_boost::io::bufwriter_floss(dout.unwrap(), ti);
                        loss.write_writer(&mut writer_loss, snvs);

                        fn compute_loss_max_theory(zs: &[f64], wls: &[f64]) -> f64 {
                            zs.iter().zip(wls.iter()).map(|(z, w)| w * z * z).sum()
                        }

                        let mut loss_diff = loss.clone();
                        let loss_diff_mut = loss_diff.inner_mut();
                        let zs = sample_weight.zs().unwrap();
                        let wls = sample_weight.wls().unwrap();
                        let loss_max_theory = compute_loss_max_theory(&zs, &wls);
                        loss_diff_mut
                            .iter_mut()
                            .for_each(|x| *x = 0.5 * (*x - loss_max_theory));

                        let mut writer_loss =
                            wgt_boost::io::bufwriter_floss_adjmax(dout.unwrap(), ti);
                        loss_diff.write_writer(&mut writer_loss, snvs);
                    }
                }
            }

            bi += 1;
            ti += 1;

            /*             // for const
                       // below of ti+=1
                       if use_const_for_loss {
                           let wgt = boosting_logit_const(
                               dataset,
                               &mut scores,
                               &mut sample_weight,
                               ti,
                               dataset_val,
                               &mut scores_val,
                           );
                           wgts.add_wgt(wgt.clone());
                           wgts.write_wgt_writer(writer);

                           ti += 1;
                       }
            */
            //if let Some(dataset_val) = dataset_val {
            if is_monitor {
                let dataset_val = dataset_val.unwrap();
                // this does not consider duplicated snvs
                // -> fixed
                let nsnv = wgts.count_unique_snv();
                let nsnv_monitor = nsnv - 1;
                if nsnvs_monitor.contains(&nsnv_monitor) {
                    log::debug!("monitor stat: nsnvs {}", wgts.count_unique_snv());
                    update_acc_monitor(
                        &mut acc_monitor,
                        &nsnvs_monitor,
                        nsnv_monitor,
                        dataset_val,
                        &scores_val_monitor,
                        //&scores_val,
                        &scores_val_cov,
                    );
                    let is_exit = monitor_acc(&mut acc_monitor, &nsnvs_monitor, nsnv_monitor);
                    if is_exit {
                        log::info!("Exit since monitoring accuracy exceeded the threshold.");
                        // TODO: create wgt with maximum acc
                        let (nsnvs_max, acc_max) =
                            compute_acc_max(&mut acc_monitor, &nsnvs_monitor);
                        return Some((nsnvs_max, acc_max));
                    }
                }
            };

            log::debug!("\n");
        }
    }

    return None;
}

fn update_acc_monitor(
    acc_monitor: &mut [f64],
    nsnvs_monitor: &[usize],
    nsnv: usize,
    dataset_val: &Dataset,
    scores_val: &[f64],
    scores_val_cov: &[f64],
) {
    //let n = dataset_val.samples().samples_n();
    let acc = pgs::nagelkerke_r2(
        &dataset_val.samples().phe().inner_i32(),
        //&dataset_val.samples().phe().inner_f64(),
        scores_val,
        //&scores_val[..n],
        scores_val_cov,
        //&scores_val_cov[..n],
    );

    let nsnv_index = nsnvs_monitor.iter().position(|x| *x == nsnv).unwrap();
    acc_monitor[nsnv_index] = acc
}

fn compute_acc_max(acc_monitor: &mut [f64], nsnvs_monitor: &[usize]) -> (usize, f64) {
    let index_max = acc_monitor
        .iter()
        .enumerate()
        .filter(|(_, a)| !a.is_nan())
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(index, _)| index)
        .unwrap();

    let nsnv_max = nsnvs_monitor[index_max];
    let acc_max = acc_monitor[index_max];
    (nsnv_max, acc_max)
}

fn monitor_acc(acc_monitor: &mut [f64], nsnvs_monitor: &[usize], nsnv: usize) -> bool {
    log::debug!("nsnvs_monitor {:?}", nsnvs_monitor);
    log::debug!("acc_monitor {:?}", acc_monitor);

    let (nsnv_max, _) = compute_acc_max(acc_monitor, nsnvs_monitor);

    /*     let index_max = acc_monitor
        .iter()
        .enumerate()
        .filter(|(_, a)| !a.is_nan())
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(index, _)| index)
        .unwrap();

    let nsnv_max = nsnvs_monitor[index_max]; */
    log::debug!("nsnvs_max {}", nsnv_max);
    log::debug!("nsnvs_current {}", nsnv);

    if !acc_monitor.last().unwrap().is_nan() {
        log::debug!("Exceeded last of # snvs in monitoring nsnvs.");
        return true;
    }

    if nsnv <= 100 {
        return false;
    }
    //if nsnv <= 1000 {
    //    return false;
    //}
    //if nsnvs > nsnvs_max+2000 {
    log::debug!("nsnvs_current {}, {}", nsnv as f64, nsnv_max as f64);
    if (nsnv as f64) > (nsnv_max as f64) * 1.5f64 {
        return true;
    }

    false
}

fn set_wgt_m<W: WgtTrait>(wgts: &mut [W], snvs: &Snvs) {
    //let sida_to_m = wgt_sida_to_dataset_m(wgts, snvs);

    // here use &str not String
    let mut sida_to_m = HashMap::with_capacity(snvs.snvs_n());
    for (mi, snvid) in snvs.snv_indexs().iter().enumerate() {
        sida_to_m.insert(snvid.sida(), mi);
    }

    for wgt in wgts.iter_mut() {
        if wgt.kind().is_snv() {
            let sida = wgt.kind().snv_index().sida();
            //log::debug!("sida {}", sida);
            //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
            if let Some(mi) = sida_to_m.get(sida) {
                wgt.set_snv_index(Some(*mi));
            } else {
                panic!("Unknown SNV in wgt.");
            }
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

/*
#[cfg(test)]
mod tests {
    use super::*;

    //use std::env;
    //use std::path::Path;

    #[test]
    fn test_boosting() {}
}
*/
