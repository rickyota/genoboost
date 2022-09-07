//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
//! Loading genotype requires the long time. Extracting genotype takes negligible time.
//!
// TODO: refactor abcd->table4, predict->genot
// TODO: change wgt_boost::Io -> io_wgt_boost
// TODO: vec::push -> resize
// split covs from samples

//mod correlation; // not implemented cmatrix yet
mod coefficient;
mod compute_pred;
mod epsilon;
mod iteration;
pub mod loss; // pub for bench
mod regression_cov;
mod sample_weight;
mod table;

use std::collections::HashMap;
use std::io::BufWriter;
use std::path::Path;
use std::time::Instant;

use super::WgtBoost;
use crate::boosting_param::{BoostParam, BoostType, LossFunc};
use crate::boosting_score::score;
use crate::wgt_boost;
use crate::wgt_boosts::WgtBoosts;
use genetics::samples::prelude::*;
use genetics::WgtTrait;
use genetics::{alloc, samples::CovsTrait, Covs, Dataset, Genot, Snvs};
use loss::LossStruct;
pub use table::ContingencyTable;

// SHOULD MOVE IN `impl WgtBoost{}`
// like https://doc.rust-lang.org/src/std/net/ip.rs.html#478
//const THRESHOLD_SNV: [f64; 2] = [0.5, 1.5];

fn boosting_iter_cov(score_ti: &mut [f64], wgtcov: &WgtBoost, covs: &Covs) -> WgtBoost {
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
    wgtcov.clone()
}

fn boosting_iter_snv(
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
) -> WgtBoost {
    println!("ps stat case/control");
    println!(
        "ps sum {:.4e}, {:.4e}",
        sample_weight::sum_ps(ps, phe, true),
        sample_weight::sum_ps(ps, phe, false)
    );
    println!(
        "ps mean {:.4e}, {:.4e}",
        sample_weight::mean_ps(ps, phe, true),
        sample_weight::mean_ps(ps, phe, false)
    );
    println!(
        "ps med {:.4e}, {:.4e}",
        sample_weight::med_ps(ps, phe, true),
        sample_weight::med_ps(ps, phe, false)
    );
    println!(
        "ps max {:.4e}, {:.4e}",
        sample_weight::max_ps(ps, phe, true),
        sample_weight::max_ps(ps, phe, false)
    );
    println!(
        "ps min {:.4e}, {:.4e}",
        sample_weight::min_ps(ps, phe, true),
        sample_weight::min_ps(ps, phe, false)
    );

    let mut wgt = loss::search_min_loss_gt(loss, ti, genot, &ps, phe, snvs, boost_param);

    compute_pred::compute_pred(pred, &wgt, genot);

    // create table
    let (table_sum, is_eps) =
        table::calculate_table_sum(&pred, ps, phe, boost_param.eps(), boost_param.boost_type());
    println!("table {:?}", table_sum);
    wgt.set_contingency_table(table_sum, is_eps);

    // TODO: just end without error when loss==1.0
    //let table = wgt.contingency_table();
    //println!("table {:?}", table);
    //let abcd_sum = wgt.abcd_sum();
    //println!("abcd {:?}", abcd_sum);

    let coef_ti = coefficient::calculate_coefficients(
        //&wgt,
        wgt.contingency_table().unwrap(),
        boost_param.boost_type(),
        boost_param.learning_rate(),
    );
    println!("coef {:?}", coef_ti);

    wgt.set_coef(coef_ti);
    println!("wgt {:?}", wgt);
    // write wgt outside
    //wgt.write_wgt(fout);

    sample_weight::renew_score(scores, pred, &wgt);
    sample_weight::renew_ws(
        ws,
        scores,
        phe,
        boost_param.loss_func(),
        boost_param.sample_weight_clip(),
    );
    sample_weight::renew_ps(ps, ws);

    let n = phe.n();
    if n < 10 {
        println!("score {:?}", scores);
        println!("ws {:?}", ws);
        println!("ps {:?}", ps);
    }

    // TODO: create argument to output ws file --output ws
    //io_wgt_boost::write_ws(fout, &ws_ti, ti);

    wgt
}

// TODO: extract common with boosting_iter_snv()
fn boosting_iter_snv_second(
    ti: usize,
    boost_param: BoostParam,
    loss: &mut LossStruct,
    predictions: &Genot,
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

    compute_pred::compute_pred(pred, &wgt_second, predictions);

    // create table
    let (table_sum, is_eps) =
        table::calculate_table_sum(&pred, ps, phe, boost_param.eps(), boost_param.boost_type());
    println!("table {:?}", table_sum);
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

pub fn boosting_covs(
    wgts: &mut WgtBoosts,
    boost_param: BoostParam,
    dataset: &Dataset,
    scores: &mut [f64],
    ws: &mut [f64],
    ps: &mut [f64],
    pred: &mut [u8],
) -> usize {
    let wgtcovs_logreg = regression_cov::logistic_regression_covs(dataset.samples());

    let covs = dataset.samples().covs().unwrap();

    let p = wgtcovs_logreg.len();
    for ti in 0..p {
        let pi = ti;
        let wgtcov = boosting_iter_cov(scores, &wgtcovs_logreg[pi], covs);
        wgts.add_wgt(wgtcov);
    }

    let phe = dataset.samples().phe();
    let predictions = dataset.genot();
    let n = dataset.samples().samples_n();

    println!("after cov");
    sample_weight::renew_ws(
        ws,
        scores,
        phe,
        boost_param.loss_func(),
        boost_param.sample_weight_clip(),
    );
    sample_weight::renew_ps(ps, ws);
    //println!("score {:?}", score_ti);
    //println!("ws {:?}", ws_ti);
    //println!("ps {:?}", ps_ti);

    // for constonly
    let ti = p;
    let wgt = loss::create_loss_const(ti, ps, pred, phe, n, boost_param);
    compute_pred::compute_pred(pred, &wgt, predictions);

    sample_weight::renew_score(scores, pred, &wgt);
    sample_weight::renew_ws(
        ws,
        scores,
        phe,
        boost_param.loss_func(),
        boost_param.sample_weight_clip(),
    );
    sample_weight::renew_ps(ps, ws);
    //println!("score {:?}", score_ti);
    //println!("ws {:?}", ws_ti);
    //println!("ps {:?}", ps_ti);

    wgts.add_wgt(wgt);
    //wgts.write_wgt_writer(writer);

    p
}

/// allow input file and stdio for test
pub fn boosting<W: std::io::Write>(
    //fout: &str,
    //file: &File,
    writer: &mut BufWriter<W>,
    iteration: usize,
    boost_param: BoostParam,
    dataset: &Dataset,
    //writer_loss: &mut BufWriter<W>,
    //fout: Option<&str>, // should Some() when is_resume=T or is_write_loss=T
    is_resume: bool,
    is_write_loss: bool,
    dout: Option<&Path>, // for is_resume or is_write_loss
                         //dloss: Option<&Path>, //is_write_loss: bool, // TODO: make different dir in boosting
) {
    let start_time = Instant::now();

    let n = dataset.samples().samples_n();
    let mut scores: Vec<f64> = vec![0.0; n];
    let mut ws: Vec<f64> = vec![1.0 / (n as f64); n];
    // len != n only for ps
    let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
    ps.resize(n + 32, 0.0f64);

    // make pred u8 -> good idea!!
    // but slow when non-simd
    let mut pred: Vec<u8> = vec![0u8; n];

    let m = dataset.snvs().snvs_n();

    let mut loss = LossStruct::new(boost_param.boost_type(), m);

    sample_weight::renew_ps(&mut ps, &ws);

    let phe = dataset.samples().phe();
    let genot = dataset.genot();
    let snvs = dataset.snvs();

    let mut wgts;

    let t_start;
    //let is_exist_wgt = wgt_boost::io::is_exist_wgt(fout);
    if is_resume & wgt_boost::io::is_exist_wgt(dout.unwrap()) {
        //if is_exist_wgt & is_resume {
        // assume cov is done when resume
        println!("resume");

        wgts = WgtBoosts::new_from_file(dout.unwrap(), boost_param.boost_type());
        set_wgt_m(wgts.wgts_mut(), snvs);

        for wgt in wgts.wgts().iter() {
            score::add_score(
                &mut scores,
                wgt,
                dataset.genot(),
                dataset.samples().covs().unwrap(),
            );
        }

        sample_weight::renew_ws(
            &mut ws,
            &scores,
            phe,
            boost_param.loss_func(),
            boost_param.sample_weight_clip(),
        );
        sample_weight::renew_ps(&mut ps, &mut ws);

        t_start = wgts.wgts().len();
    } else {
        wgts = WgtBoosts::new(boost_param.boost_type());
        let p = boosting_covs(
            &mut wgts,
            boost_param,
            dataset,
            &mut scores,
            &mut ws,
            &mut ps,
            &mut pred,
        );

        wgt_boost::io::write_cols(writer, boost_param.boost_type());
        wgts.write_wgt_writer(writer);

        t_start = p + 1;
    }

    //print max,min of ws
    println!(
        "ps max,min {},{}",
        ps[..n].iter().fold(0.0 / 0.0, |v, v1| v1.max(v)),
        ps[..n].iter().fold(0.0 / 0.0, |v, v1| v1.min(v)),
    );
    //let max = iter.fold(0.0 / 0.0, f64::max); //somehow error
    //let max = iter.fold(0.0 / 0.0, |m, v| v.max(m));

    //for ti in p + 1..iteration {
    //for ti in t_start..iteration {

    let mut ti = t_start;
    while ti < iteration {
        println!(
            "\nIteration: {}, {} sec",
            ti,
            start_time.elapsed().as_secs()
        );

        let wgt = boosting_iter_snv(
            ti,
            boost_param,
            &mut loss,
            genot,
            &mut pred,
            &mut scores,
            &mut ws,
            &mut ps,
            phe,
            snvs,
        );
        wgts.add_wgt(wgt.clone());
        //wgts.add_wgt(wgt);

        //println!("wgts {:?}", wgts);
        wgts.write_wgt_writer(writer);

        ti += 1;

        if is_write_loss {
            if (ti < 200) || (ti % 100 == 0) {
                let mut writer_loss = wgt_boost::io::bufwriter_floss(dout.unwrap(), ti);
                loss.write_writer(&mut writer_loss, snvs);
            }
        }

        if boost_param.is_dom_rec() {
            let wgt_second = boosting_iter_snv_second(
                ti,
                boost_param,
                &mut loss,
                genot,
                &mut pred,
                &mut scores,
                &mut ws,
                &mut ps,
                phe,
                snvs,
                &wgt,
            );
            wgts.add_wgt(wgt_second);

            //println!("wgts {:?}", wgts);
            wgts.write_wgt_writer(writer);

            ti += 1;
        }
    }
}

// TODO: mv
//fn wgt_sida_to_dataset_m<W: WgtTrait>(wgts: &mut [W], snvs: &Snvs) -> HashMap<String, usize> {
//}

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
            //println!("sida {}", sida);
            //println!("m_in {:?}", sida_in_to_m_in.get(sida));
            if let Some(mi) = sida_to_m.get(sida) {
                wgt.set_snv_index(Some(*mi));
            } else {
                panic!("Unknown SNV in wgt.")
            }
        }

        // do nth for cov
        // not found
        //wgt_to_m.insert(wgti, None);
    }
}

/*
pub fn boosting(
    fout: &str,
    iter: usize,
    boost_type: &BoostType,
    predictions: &[B8],
    ys: &[B8],
    m: usize,
    n: usize,
    snvs: &[Snv],
    wgtcovs: &[WgtBoost],
    //wgtcovs: &[Wgt],
    is_dom_rec: bool,
) {
    let mut score_ti: Vec<f64> = Vec::with_capacity(n);
    vec::push(&mut score_ti, 0.0, n);

    let mut ws_ti: Vec<f64> = Vec::with_capacity(n);
    vec::push(&mut ws_ti, 1.0 / (n as f64), n);

    // len != n only for ps
    let mut ps_ti: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
    vec::push(&mut ps_ti, 0.0, n + 32);

    let len_n = predict::len_n(n);
    //let len_n = mistake::get_len_n(n);

    let mut pred: Vec<B8> = Vec::with_capacity(len_n);
    vec::push(&mut pred, 0x00, len_n);

    let mut losss_gt: Vec<f64> = Vec::with_capacity(2 * m);
    vec::push(&mut losss_gt, f64::NAN, 2 * m);

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

    println!("after cov");
    //sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
    sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
    sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
    //println!("score {:?}", score_ti);
    //println!("ws {:?}", ws_ti);
    //println!("ps {:?}", ps_ti);

    // for constonly
    let ti = p;
    let wgt = loss::create_loss_const(ti, &ps_ti, &mut pred, &ys, n);
    wgt.write_wgt(fout);
    compute::compute_pred(&mut pred, &wgt, predictions, n);

    sample_weight::renew_score(&mut score_ti, &pred, &wgt);
    sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
    //sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
    sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
    //println!("score {:?}", score_ti);
    //println!("ws {:?}", ws_ti);
    //println!("ps {:?}", ps_ti);

    //print max,min of ws
    println!(
        "ps max,min {},{}",
        ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.max(v)),
        ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.min(v)),
    );
    //let max = iter.fold(0.0 / 0.0, f64::max); //somehow error
    //let max = iter.fold(0.0 / 0.0, |m, v| v.max(m));

    for ti in p + 1..iter {
        println!("\nIteration: {}", ti);
        let mut wgt = loss::search_min_loss_gt(&mut losss_gt, ti, predictions, &ps_ti, ys, n, snvs);
        /*
        let wgt_gt = loss::search_min_loss_gt(&mut losss_gt, ti, predictions, &ps_ti, ys, n, snvs);
        let mut wgt = wgt_gt;
         */

        // TODO: just end without error when loss==1.0

        let abcd_sum = wgt.abcd_sum();
        println!("abcd {:?}", abcd_sum);

        let coef_ti = coefficient::calculate_coefficients(&wgt, boost_type);
        println!("const, alpha: {:.4}, {:.4}", coef_ti.0, coef_ti.1);
        wgt.wgt_mut().model_mut().set_coef_binary(coef_ti);
        wgt.write_wgt(fout);

        compute::compute_pred(&mut pred, &wgt, predictions, n);

        sample_weight::renew_score(&mut score_ti, &pred, &wgt);
        sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
        sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
        //println!("score {:?}", score_ti);
        //println!("ws {:?}", ws_ti);
        //println!("ps {:?}", ps_ti);

        println!(
            "ps max,min {},{}",
            ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.max(v)),
            ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.min(v)),
        );

        // TODO: create argument to output ws file --output ws
        io_wgt_boost::write_ws(fout, &ws_ti, ti);

        if is_dom_rec {
            println!("Choose the other inheritance model.");

            let mut wgt_second =
                loss::search_wgt_gt_second(&wgt, &losss_gt, ti, predictions, &ps_ti, ys, n, snvs);
            //let mut wgt = loss::search_min_loss_gt(&mut losss_gt, ti, predictions, &ps_ti, ys, n, snvs);

            let abcd_sum = wgt_second.abcd_sum();
            println!("abcd {:?}", abcd_sum);

            let coef_ti = coefficient::calculate_coefficients(&wgt_second, boost_type);
            println!("const, alpha: {:.4}, {:.4}", coef_ti.0, coef_ti.1);
            wgt_second.wgt_mut().model_mut().set_coef_binary(coef_ti);
            wgt_second.write_wgt(fout);

            compute::compute_pred(&mut pred, &wgt_second, predictions, n);

            sample_weight::renew_score(&mut score_ti, &pred, &wgt_second);
            sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
            sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
            //println!("score {:?}", score_ti);
            //println!("ws {:?}", ws_ti);
            //println!("ps {:?}", ps_ti);

            println!(
                "ps max,min {},{}",
                ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.max(v)),
                ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.min(v)),
            );

            io_wgt_boost::write_ws(fout, &ws_ti, ti);
        }
    }
}
*/

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

    println!("after cov");
    sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
    //sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
    sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
    //println!("score {:?}", score_ti);
    //println!("ws {:?}", ws_ti);
    //println!("ps {:?}", ps_ti);

    // for constonly
    let ti = p;
    let wgt = loss::create_loss_const(ti, &ps_ti, &mut pred, &ys, n);
    wgt.write_wgt(fout);
    compute::compute_pred(&mut pred, &wgt, predictions, n);

    sample_weight::renew_score(&mut score_ti, &pred, &wgt);
    sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
    //sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
    sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
    //println!("score {:?}", score_ti);
    //println!("ws {:?}", ws_ti);
    //println!("ps {:?}", ps_ti);

    //print max,min of ws
    println!(
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
        println!("\nIteration: {}", ti);

        //println!("is_pruned {:?}", is_pruned);

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
        println!("abcd {:?}", abcd_sum);

        let coef_ti = coefficient::calculate_coefficients(&wgt, boost_type);
        println!("alpha, const {:?}", coef_ti);
        //println!("const, alpha: {:.4}, {:.4}", coef_ti.0, coef_ti.1);
        wgt.wgt_mut().model_mut().set_coef_binary(coef_ti);
        wgt.write_wgt(fout);

        compute::compute_pred(&mut pred, &wgt, predictions, n);

        sample_weight::renew_score(&mut score_ti, &pred, &wgt);
        sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
        //sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type, n);
        sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
        //println!("score {:?}", score_ti);
        //println!("ws {:?}", ws_ti);
        //println!("ps {:?}", ps_ti);

        println!(
            "ps max,min {},{}",
            ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.max(v)),
            ps_ti[..n].iter().fold(0.0 / 0.0, |v, v1| v1.min(v)),
        );

        io_wgt_boost::write_ws(fout, &ws_ti, ti);

        if is_dom_rec {
            println!("Choose the other inheritance model.");

            let mut wgt_second =
                loss::search_wgt_gt_second(&wgt, &losss_gt, ti, predictions, &ps_ti, ys, n, snvs);

            let abcd_sum = wgt_second.abcd_sum();
            println!("abcd {:?}", abcd_sum);

            let coef_ti = coefficient::calculate_coefficients(&wgt_second, boost_type);
            println!("alpha, const {:?}", coef_ti);
            //println!("const, alpha: {:.4}, {:.4}", coef_ti.0, coef_ti.1);
            wgt_second.wgt_mut().model_mut().set_coef_binary(coef_ti);
            wgt_second.write_wgt(fout);

            compute::compute_pred(&mut pred, &wgt_second, predictions, n);

            sample_weight::renew_score(&mut score_ti, &pred, &wgt_second);
            sample_weight::renew_ws(&mut ws_ti, &score_ti, ys, boost_type);
            sample_weight::renew_ps(&mut ps_ti, &ws_ti, n);
            //println!("score {:?}", score_ti);
            //println!("ws {:?}", ws_ti);
            //println!("ps {:?}", ps_ti);

            println!(
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
            println!("All SNVs are pruned at {}-th iteration.", ti);
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

    println!("after cov");

    //println!("snvs {:?}", snvs);

    //let snv_indexs: Vec<SnvIndex> = sum_stats.iter().map(|v| v.snv_index().clone()).collect();

    for ti in p + 1..iter {
        println!("\nIteration: {}", ti);

        //println!("is_pruned {:?}", is_pruned);

        //let wgt_gt = loss::search_min_loss_gt_pruned(
        let wgt = loss::search_min_loss_gt_pruned_loss_ss(sum_stats, ti, &is_pruned);

        wgt.write_wgt(fout);

        if is_dom_rec {
            println!("Choose the other inheritance model.");

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
            println!("All SNVs are pruned at {}-th iteration.", ti);
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
