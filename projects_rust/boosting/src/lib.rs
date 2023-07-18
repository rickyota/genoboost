//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
//

mod boosting_param;
mod boosting_score;
pub mod boosting_train; // pub for bench
mod wgt_boost;
mod wgt_boosts;

//use itertools::Itertools;
use rand::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use mysmartcore::linalg::basic::arrays::Array2;
// mysmartcore v0.3.1
// https://smartcorelib.org/user_guide/model_selection.html
use mysmartcore::linalg::basic::matrix::DenseMatrix;
use mysmartcore::model_selection::KFold;
use mysmartcore::model_selection::*;

pub use boosting_param::{
    BoostMethod, BoostParam, BoostParams, BoostParamsTypes, BoostType, EffEps, Eps, IterationNumber,
};
use boosting_train::{loss::LossStruct, ContingencyTable};
use genetics::snv;
use genetics::textfile;
use genetics::GenotFormat;
use genetics::{io_genot, vec, Dataset, Snvs};
use wgt_boost::WgtBoost;
use wgt_boosts::WgtBoosts;

use crate::boosting_train::sample_weight::SampleWeight;

#[macro_use]
extern crate assert_float_eq;

/// Run GenoBoost. Input plink file.
pub fn run_boosting(
    dout: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    cov_name: &str,
    boost_method: BoostMethod,
    boost_params: &BoostParams,
    fin_snv: Option<&Path>,
    fin_sample: Option<&Path>,
    //fin_cov: Option<&Path>,
    fin_sample_val: Option<&Path>,
    //is_dom_rec: bool,
    use_adjloss: bool,
    use_const_for_loss: bool,
    is_resume: bool,
    is_write_loss: bool,
    prune_snv: Option<f64>,
    //learning_rates: &[f64],
    is_monitor: bool,
    //nsnvs_monitor: Option<Vec<usize>>,
) {
    // check fwgt does not exist.
    if !is_resume {
        for learning_rate in boost_params.learning_rates().iter() {
            wgt_boost::io::check_file_wgt_not_exist(dout, *learning_rate);
        }
    }

    // TODO: if all wgt exceeds #SNVs, then exit here.

    let phe_buf = fin_phe.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
    let snv_buf = fin_snv.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
    let sample_buf = fin_sample.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
    let sample_val_buf =
        fin_sample_val.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());

    io_genot::check_valid_fin(fin, gfmt);

    let start_time = Instant::now();
    let (dataset, dataset_val) = load_dataset_boosting(
        dout,
        fin,
        gfmt,
        phe_buf.as_deref(),
        //fin_phe,
        phe_name,
        cov_name,
        boost_params,
        snv_buf.as_deref(),
        //fin_snv,
        sample_buf.as_deref(),
        //fin_sample,
        sample_val_buf.as_deref(),
        //fin_sample_val,
        use_adjloss,
        prune_snv,
    );
    log::info!("Created dataset: {} sec", start_time.elapsed().as_secs());

    wgt_boost::io::create_dir(&dout);

    let mut nsvs_acc_use: Vec<Option<(usize, f64)>> = vec![];
    for boost_param_lr in boost_params.clone().into_iter() {
        let learning_rate = boost_param_lr.learning_rate();
        log::debug!("learning rate: {:?}", learning_rate);
        // add --integrate
        let dout_para = wgt_boost::io::get_dname_para(dout, learning_rate);
        wgt_boost::io::create_dir(&dout_para);
        let mut writer = wgt_boost::io::bufwriter_fwgt_append(&dout_para);
        let nsnv_acc_use = if boost_param_lr.batch_way().is_none() {
            match boost_method {
                BoostMethod::Classic => {
                    log::info!("Run boosting");
                    boosting_train::boosting(
                        &mut writer,
                        boost_param_lr,
                        &dataset,
                        use_adjloss,
                        is_resume,
                        is_write_loss,
                        Some(&dout_para),
                    )
                }
                _ => unimplemented!(),
            }
        } else {
            match boost_method {
                BoostMethod::Classic => {
                    log::info!("Run boosting");
                    boosting_train::boosting_batch(
                        &mut writer,
                        boost_param_lr,
                        &dataset,
                        dataset_val.as_ref(),
                        use_adjloss,
                        use_const_for_loss,
                        is_resume,
                        is_write_loss,
                        Some(&dout_para),
                        is_monitor,
                        None,
                    )
                }
                _ => unimplemented!(),
            }
        };
        nsvs_acc_use.push(nsnv_acc_use);
    }

    if nsvs_acc_use.iter().any(|&x| x.is_none()) {
        log::info!("Do not create the most accurate weight file across learning rates since training was suspended in some learning rates.")
    } else {
        // TODO: create best .wgt across lr
        let learning_rates = boost_params.learning_rates();
        let (nsnv_acc_max, lr) = nsvs_acc_use
            .iter()
            .map(|x| x.unwrap())
            .zip(learning_rates.iter())
            .max_by(|((_, acc_a), _), ((_, acc_b), _)| acc_a.partial_cmp(acc_b).unwrap())
            .map(|((nsnv, _), lr)| (nsnv, lr))
            .expect("Cannot determine the highest accuracy.");
        log::info!(
            "The best parameter is #SNVs={:?}, learning rate={:?}.",
            nsnv_acc_max,
            lr
        );

        let fpara_best = wgt_boost::io::get_file_wgt_best_para(dout);
        let para_best_string = create_best_para_buf(nsnv_acc_max, *lr);
        write_file_string(&fpara_best, para_best_string);

        let wgt_best_string = create_best_wgt_buf(dout, nsnv_acc_max, *lr, None);
        let fwgt_best = wgt_boost::io::get_file_wgt_best(dout);
        write_file_string(&fwgt_best, wgt_best_string);
    }

    /*
    match boost_method {
        BoostMethod::Classic => {
            log::info!("Run boosting");
            boosting::boosting(
                //&file,
                &mut writer,
                iteration,
                boost_param,
                &dataset_ext,
                //&mut writer_loss,
                //fout,
                is_resume,
                is_write_loss,
                Some(dout), //dout_option,
                            //dloss,
            )
            //                 is_dom_rec,
            //boosting::boosting(fout, iteration, boost_param, &dataset)
        }

        BoostMethod::Pruning(_) => {
            unimplemented!();
        }
        //BoostMethod::Pruning(clump_r2) => {
        //    log::debug!("Run boosting pruning");
        //    boosting::boosting_pruning(
        //        fout,
        //        iteration,
        //        boost_type,
        //        &predictions,
        //        &ys,
        //        m,
        //        n,
        //        &snvs,
        //        &wgtcovs,
        //        is_dom_rec,
        //        clump_r2,
        //    )
        //}
        BoostMethod::Ss(_) => {
            unimplemented!();
        } //BoostMethod::Ss(clump_r2) => {
          //    log::debug!("Run boosting ss.");
          //    // load loss
          //    //summary_statistics::set_ss_loss(fin_loss.unwrap(), &mut snvs);
          //    let sum_stats =
          //        sum_stat::io::convert_sum_stat_from_snv_consume(fin_loss.unwrap(), snvs);

          //    boosting::boosting_ss(
          //        fout,
          //        iteration,
          //        &predictions,
          //        m,
          //        n,
          //        &sum_stats,
          //        &wgtcovs,
          //        is_dom_rec,
          //        clump_r2,
          //    )
          //}
    }
    */
}

fn extract_sample_cvi(sample_idx_cvi: &[usize], samples_in: &[String]) -> Vec<String> {
    sample_idx_cvi
        .iter()
        .map(|idx| samples_in[*idx].clone())
        .collect()
}
fn sample_string_to_buf(sample_string: Vec<String>) -> Vec<u8> {
    sample_string.join("\n").as_bytes().to_vec()
}

pub fn run_boosting_integrate_cv(
    dout: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    cov_name: &str,
    boost_method: BoostMethod,
    boost_params_types: &BoostParamsTypes,
    fin_snv: Option<&Path>,
    fin_sample: Option<&Path>,
    //fin_cov: Option<&Path>,
    fin_sample_val: Option<&Path>,
    //is_dom_rec: bool,
    use_adjloss: bool,
    use_const_for_loss: bool,
    is_resume: bool,
    is_write_loss: bool,
    prune_snv: Option<f64>,
    //learning_rates: &[f64],
    is_monitor: bool,
    cross_validation: Option<usize>,
    seed: Option<u64>,
) {
    match cross_validation {
        Some(cvn) => {
            let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
            // create cv samples
            let sample_idx_cvs: Vec<(Vec<usize>, Vec<usize>)> = if cvn == 1 {
                if let Some(seed) = seed {
                    StdRng::seed_from_u64(seed);
                }
                let mut vec: Vec<usize> = (0..n_in).collect();
                vec.shuffle(&mut rand::thread_rng());
                // TODO: mv up
                let prop_tr = 0.8;
                let n_tr = (prop_tr * (n_in as f64)) as usize;
                vec![(vec[..n_tr].to_vec(), vec[n_tr..].to_vec())]
            } else {
                let sample_matrix: DenseMatrix<f64> = DenseMatrix::zeros(n_in, 1);
                KFold::default()
                    .with_n_splits(cvn)
                    .with_shuffle(true)
                    .with_seed(seed)
                    .split(&sample_matrix)
                    .into_iter()
                    .collect()
            };

            let samples_in: Vec<String> = io_genot::load_samples_id(fin, gfmt, None);

            // run each
            let phe_buf = fin_phe.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            let snv_buf = fin_snv.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            for cvi in 0..cvn {
                //let (sample_idx_cvi, sample_idx_val_cvi) = cv.next().unwrap();
                let (sample_idx_cvi, sample_idx_val_cvi) = &sample_idx_cvs[cvi];

                // create cv-specific dir
                let dout_cv = wgt_boost::io::get_dir_cv(dout, cvi);
                let sample_string = extract_sample_cvi(&sample_idx_cvi, &samples_in);
                let sample_val_string = extract_sample_cvi(&sample_idx_val_cvi, &samples_in);

                let sample_buf = sample_string_to_buf(sample_string);
                let sample_val_buf = sample_string_to_buf(sample_val_string);

                run_boosting_integrate(
                    &dout_cv,
                    &fin,
                    gfmt,
                    phe_buf.as_deref(),
                    //fin_phe.as_deref(),
                    phe_name.as_deref(),
                    &cov_name,
                    boost_method,
                    &boost_params_types,
                    snv_buf.as_deref(),
                    //fin_snv,
                    Some(&sample_buf),
                    //fin_sample,
                    Some(&sample_val_buf),
                    //fin_sample_val,
                    //fin_snv.as_deref(),
                    //fin_sample.as_deref(),
                    //fin_cov.as_deref(),
                    //fin_sample_val.as_deref(),
                    use_adjloss,
                    use_const_for_loss,
                    is_resume,
                    is_write_loss,
                    prune_snv,
                    //&learning_rates,
                    true, //is_monitor,
                );
            }
        }
        None => {
            // TODO:
            // load once first for stdin input
            let phe_buf = fin_phe.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            let snv_buf = fin_snv.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            let sample_buf =
                fin_sample.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            let sample_val_buf =
                fin_sample_val.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            // Genoboost
            run_boosting_integrate(
                &dout,
                &fin,
                gfmt,
                phe_buf.as_deref(),
                //fin_phe.as_deref(),
                phe_name.as_deref(),
                &cov_name,
                boost_method,
                &boost_params_types,
                snv_buf.as_deref(),
                //fin_snv,
                sample_buf.as_deref(),
                //fin_sample,
                sample_val_buf.as_deref(),
                //fin_sample_val,
                //fin_snv.as_deref(),
                //fin_sample.as_deref(),
                //fin_cov.as_deref(),
                //fin_sample_val.as_deref(),
                use_adjloss,
                use_const_for_loss,
                is_resume,
                is_write_loss,
                prune_snv,
                //&learning_rates,
                is_monitor,
            );
        }
    }
}

// TODO: concat with run_boosting()
/// Integrate Genoboost. Select best model from add and non-add
pub fn run_boosting_integrate(
    dout: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    phe_buf: Option<&[u8]>,
    //fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    cov_name: &str,
    boost_method: BoostMethod,
    boost_params_types: &BoostParamsTypes,
    snv_buf: Option<&[u8]>,
    //fin_snv: Option<&Path>,
    sample_buf: Option<&[u8]>,
    //fin_sample: Option<&Path>,
    //fin_cov: Option<&Path>,
    sample_val_buf: Option<&[u8]>,
    //fin_sample_val: Option<&Path>,
    //is_dom_rec: bool,
    use_adjloss: bool,
    use_const_for_loss: bool,
    is_resume: bool,
    is_write_loss: bool,
    prune_snv: Option<f64>,
    //learning_rates: &[f64],
    is_monitor: bool,
) {
    // TODO: if all wgt exceeds #SNVs, then exit here.
    io_genot::check_valid_fin(fin, gfmt);

    assert_eq!(sample_val_buf.is_some(), is_monitor);

    let mut nsnvss_acc_use: Vec<Vec<Option<(usize, f64)>>> = vec![];
    for boost_params in boost_params_types.clone().into_iter() {
        log::info!("boost type: {:?}", boost_params.boost_type());
        // check fwgt does not exist.
        if !is_resume {
            for learning_rate in boost_params_types.learning_rates().iter() {
                wgt_boost::io::check_file_wgt_not_exist(dout, *learning_rate);
            }
        }

        let start_time = Instant::now();
        // should be inside of boost type since prune_snv differs
        let (dataset, dataset_val) = load_dataset_boosting(
            dout,
            fin,
            gfmt,
            phe_buf.as_deref(),
            //fin_phe,
            phe_name,
            cov_name,
            &boost_params,
            snv_buf.as_deref(),
            //fin_snv,
            sample_buf.as_deref(),
            //fin_sample,
            sample_val_buf.as_deref(),
            //fin_sample_val,
            use_adjloss,
            prune_snv,
        );
        log::info!("Created dataset: {} sec", start_time.elapsed().as_secs());

        wgt_boost::io::create_dir(&dout);

        let mut nsnvs_acc_use: Vec<Option<(usize, f64)>> = vec![];
        for boost_param_lr in boost_params.clone().into_iter() {
            let learning_rate = boost_param_lr.learning_rate();
            log::info!("learning rate: {:?}", learning_rate);
            let dout_para = wgt_boost::io::get_dname_para_integrate(
                dout,
                learning_rate,
                boost_param_lr.boost_type(),
            );
            wgt_boost::io::create_dir(&dout_para);
            let mut writer = wgt_boost::io::bufwriter_fwgt_append(&dout_para);
            let nsnv_acc_use = if boost_param_lr.batch_way().is_none() {
                match boost_method {
                    BoostMethod::Classic => {
                        log::info!("Run boosting");
                        boosting_train::boosting(
                            &mut writer,
                            boost_param_lr,
                            &dataset,
                            use_adjloss,
                            is_resume,
                            is_write_loss,
                            Some(&dout_para),
                        )
                    }
                    _ => unimplemented!(),
                }
            } else {
                match boost_method {
                    BoostMethod::Classic => {
                        log::info!("Run boosting batch");
                        boosting_train::boosting_batch(
                            &mut writer,
                            boost_param_lr,
                            &dataset,
                            dataset_val.as_ref(),
                            use_adjloss,
                            use_const_for_loss,
                            is_resume,
                            is_write_loss,
                            Some(&dout_para),
                            is_monitor,
                            None,
                        )
                    }
                    _ => unimplemented!(),
                }
            };
            // TODO: should write all acc not only best nsnv
            //if let Some(nsnv_acc_use) = nsnv_acc_use {
            //    let facc = wgt_boost::io::get_file_acc(&dout_para);
            //    let nsnv_acc_string = create_acc_string(nsnv_acc_use);
            //    write_file_string(&facc, nsnv_acc_string);
            //}
            nsnvs_acc_use.push(nsnv_acc_use);
        }
        nsnvss_acc_use.push(nsnvs_acc_use);
    }

    if is_monitor {
        log::debug!("nsnvss_acc_use {:?}", nsnvss_acc_use);
        best_wgt_integrate(nsnvss_acc_use, boost_params_types, dout);
    }
}

fn best_wgt_integrate(
    nsnvss_acc_use: Vec<Vec<Option<(usize, f64)>>>,
    boost_params_types: &BoostParamsTypes,
    dout: &Path,
) {
    let (mut nsnv_max, mut lr_max, mut boost_type_max): (
        Option<usize>,
        Option<f64>,
        Option<BoostType>,
    ) = (None, None, None);
    let mut acc_max: Option<f64> = None;
    for (boost_params, nsnvs_acc_use) in boost_params_types
        .clone()
        .into_iter()
        .zip(nsnvss_acc_use.iter())
    {
        for (boost_param_lr, nsnv_acc_use) in
            boost_params.clone().into_iter().zip(nsnvs_acc_use.iter())
        {
            //if nsnvs_acc_use.iter().any(|&x| x.is_none()) {
            // if any is None, do not create
            match nsnv_acc_use {
                None => {
                    log::info!("Do not create the most accurate weight file across learning rates since training was suspended in some learning rates.");
                    return;
                }
                Some(nsnv_acc_use) => {
                    let (nsnv, acc) = nsnv_acc_use;
                    if acc_max.is_none()
                        || (acc.partial_cmp(&acc_max.unwrap()).unwrap()
                            == std::cmp::Ordering::Greater)
                    {
                        acc_max = Some(*acc);
                        nsnv_max = Some(*nsnv);
                        lr_max = Some(boost_param_lr.learning_rate());
                        boost_type_max = Some(boost_param_lr.boost_type());
                    }
                }
            }
        }
    }

    if acc_max.is_none() {
        log::info!("Cannot determine the highest accuracy.");
        return;
    }
    //let acc_max = acc_max.unwrap();
    let nsnv_max = nsnv_max.unwrap();
    let lr_max = lr_max.unwrap();
    let boost_type_max = boost_type_max.unwrap();

    log::info!(
        "The best parameter is #SNVs={:?}, learning rate={:?}.",
        nsnv_max,
        lr_max
    );

    let fpara_best = wgt_boost::io::get_file_wgt_best_para(dout);
    let para_best_string = create_best_para_integrate_buf(nsnv_max, lr_max, boost_type_max);
    write_file_string(&fpara_best, para_best_string);

    let wgt_best_string = create_best_wgt_buf(dout, nsnv_max, lr_max, Some(boost_type_max));
    let fwgt_best = wgt_boost::io::get_file_wgt_best(dout);
    write_file_string(&fwgt_best, wgt_best_string);
}

fn create_best_wgt_buf(
    dout: &Path,
    nsnv_acc_max: usize,
    lr: f64,
    boost_type: Option<BoostType>,
) -> String {
    let dout_para = match boost_type {
        Some(boost_type) => wgt_boost::io::get_dname_para_integrate(dout, lr, boost_type),
        None => wgt_boost::io::get_dname_para(dout, lr),
    };
    //wgt_boost::io::create_dir(&dout_para);
    let fwgt = wgt_boost::io::get_fname_wgt(&dout_para);

    let wgts = WgtBoosts::new_from_file(&fwgt);

    // get iteration number
    let mut itern: Option<usize> = None;
    //let mut itern: usize;

    let mut set_rs = HashSet::new();
    for (wgti, wgt) in wgts.wgts().iter().enumerate() {
        if !wgt.wgt().kind().is_snv() {
            continue;
        }
        let rs = wgt.wgt().kind().snv_index().rs().to_string();
        set_rs.insert(rs);

        // first get wgt of nsnv+1
        if set_rs.len() == nsnv_acc_max + 1 {
            // iter is just before
            itern = Some(wgti - 1);
            break;
        }
    }

    if itern.is_none() {
        panic!("Error in creating best_wgt_buf.");
    }

    // +1 for header
    let linen = itern.unwrap() + 1;

    let buf_wgt = textfile::read_file_to_end(&fwgt, None).unwrap();

    let lines = BufReader::new(&buf_wgt[..]).lines().take(linen);
    let wgt_best_vec_string: Vec<String> = lines.into_iter().map(|x| x.unwrap()).collect();

    let wgt_best_string = wgt_best_vec_string.join("\n");

    wgt_best_string
}

fn write_file_string(fin: &Path, wgt_best: String) {
    let file = match File::create(&fin) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {:?}",
            &fin
        ),
    };

    let mut writer = BufWriter::new(file);
    writer.write(wgt_best.as_bytes()).unwrap();
}

fn create_best_para_buf(nsnv_acc_max: usize, lr: f64) -> String {
    /*     fn lr_to_str(lr: &Option<f64>) -> String {
        match lr {
            None => "1.0".to_string(),
            Some(x) => x.to_string(),
        }
    } */
    let wgt_best_para = vec![
        String::from("nsnv\t") + nsnv_acc_max.to_string().as_str(),
        String::from("learning_rate\t") + lr.to_string().as_str(),
    ];

    let wgt_best_para = wgt_best_para.join("\n");

    wgt_best_para
}

fn create_best_para_integrate_buf(nsnv_acc_max: usize, lr: f64, boost_type: BoostType) -> String {
    let wgt_best_para = vec![
        String::from("nsnv\t") + nsnv_acc_max.to_string().as_str(),
        String::from("learning_rate\t") + lr.to_string().as_str(),
        String::from("genetic_model\t") + boost_type.genetic_model().as_str(),
    ];

    let wgt_best_para = wgt_best_para.join("\n");

    wgt_best_para
}

fn load_dataset_boosting(
    dout: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    //fin_phe: Option<&Path>,
    phe_buf: Option<&[u8]>,
    phe_name: Option<&str>,
    cov_name: &str,
    boost_params: &BoostParams,
    snv_buf: Option<&[u8]>,
    //fin_snv: Option<&Path>,
    sample_buf: Option<&[u8]>,
    //fin_sample: Option<&Path>,
    //fin_cov: Option<&Path>,
    sample_val_buf: Option<&[u8]>,
    //fin_sample_val: Option<&Path>,
    use_adjloss: bool,
    prune_snv: Option<f64>,
) -> (Dataset, Option<Dataset>) {
    //let phe_buf = fin_phe.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
    //let snv_buf = fin_snv.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
    //let sample_buf = fin_sample.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
    //let sample_val_buf = fin_sample_val.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());

    // FIXME: check if in .fam or .phe, col=0 is not duplicated
    let boost_param = boost_params.param_lr_none();
    let use_missing = boost_param.boost_type().use_missing();
    // create dataset
    // extract snvs by loss function
    if let Some(prop_prune_snv) = prune_snv {
        log::info!("Prune SNVs by decreasing loss: {}", prop_prune_snv);
        let start = Instant::now();

        sample_val_buf.expect("Not Implemented");
        //fin_sample_val.expect("Not Implemented");

        // TODO: better
        // to get m
        let m;
        {
            //let m_in: usize = plink::compute_num_snv(fin, gfmt).unwrap();
            let snvs_in = io_genot::load_snvs(fin, gfmt);
            //(m, _) = snv::make_use_snvs(fin_snv, &snvs_in);
            (m, _) = snv::make_use_snvs_buf(snv_buf, &snvs_in);
        }

        // TODO: depends on available memory
        //let m_range = 200_000usize;
        //let m_range = 400_000usize;
        let m_range = 1_700_000usize;
        //let m_range = 2_000_000usize;
        //let m_range = 800_000usize;

        let mut losss = Vec::new();
        for snv_i_start in (0..m).step_by(m_range) {
            let mut filt_snv = vec![false; m];
            let snv_i_end = (snv_i_start + m_range).min(m);
            let m_ = snv_i_end - snv_i_start;
            filt_snv[snv_i_start..snv_i_end].fill(true);
            assert_eq!(vec::count_true(&filt_snv), m_);

            let dataset = Dataset::new(
                fin,
                gfmt,
                phe_buf.as_deref(),
                phe_name,
                cov_name,
                snv_buf.as_deref(),
                sample_buf.as_deref(),
                use_missing,
                Some(&filt_snv),
                None,
            );

            /*             let dataset = Dataset::new_old(
                           fin,
                           gfmt,
                           fin_phe,
                           phe_name,
                           cov_name,
                           fin_snv,
                           fin_sample,
                           fin_cov,
                           use_missing,
                           Some(&filt_snv),
                           None,
                       );
            */
            //TODO: mv boostmethod::classic
            let n = dataset.samples().samples_n();
            //let mut pred: Vec<u8> = vec![0u8; n];

            let mut scores: Vec<f64> = vec![0.0; n];

            let mut sample_weight = SampleWeight::new(
                n,
                boost_param.boost_type(),
                boost_param.loss_func(),
                boost_param.sample_weight_clip(),
                boost_param.sample_weight_wls_clip(),
            );
            sample_weight.renew_sample_weight(&scores, dataset.samples().phe());

            let mut wgts = WgtBoosts::new(boost_param.boost_type());
            let _ = boosting_train::boosting_covs(
                &mut wgts,
                &dataset,
                &mut scores,
                &mut sample_weight,
                0,
                None,
                &mut [0.0f64; 0],
            );

            let mut loss = LossStruct::new(boost_param.boost_type(), m_);

            boosting_train::loss::calculate_loss_gt(
                &mut loss,
                dataset.genot(),
                &sample_weight,
                dataset.samples().phe(),
                boost_param,
                &HashSet::<usize>::new(),
                use_adjloss,
            );

            losss.push(loss);
        }

        let mut v = Vec::new();
        for loss_ in losss {
            //let LossStruct::ModelFree(v_in, _) = loss_;
            let v_in = loss_.inner();
            v.extend(v_in);
            //v.extend(loss_.inner());
        }

        let loss = LossStruct::new_vec(boost_param.boost_type(), v);
        let mut writer_loss = wgt_boost::io::bufwriter_floss(dout, 0);
        let snvs = Snvs::new_plink(fin, gfmt);
        loss.write_writer(&mut writer_loss, &snvs);

        let use_snvs_loss = loss.search_topprop(prop_prune_snv);

        log::debug!("created use_snvs_loss {}", use_snvs_loss.len());

        let dataset = Dataset::new(
            fin,
            gfmt,
            phe_buf.as_deref(),
            //fin_phe,
            phe_name,
            cov_name,
            snv_buf.as_deref(),
            sample_buf.as_deref(),
            //fin_snv,
            //fin_sample,
            use_missing,
            Some(&use_snvs_loss),
            None,
        );

        log::debug!(
            "It took {} seconds to prune Dataset.",
            start.elapsed().as_secs()
        );

        //dataset_ext = dataset.extract_snvs(&use_snvs_loss);
        (dataset, None)
    } else {
        //let phe_buf = fin_phe.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
        //let snv_buf = fin_snv.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
        //let sample_buf = fin_sample.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
        let dataset = Dataset::new(
            fin,
            gfmt,
            phe_buf.as_deref(),
            //fin_phe,
            phe_name,
            cov_name,
            snv_buf.as_deref(),
            sample_buf.as_deref(),
            //fin_snv,
            //fin_sample,
            use_missing,
            None,
            None,
        );
        //dataset_ext = dataset;

        //let dataset_val = if fin_sample_val.is_some() {
        let dataset_val = if sample_val_buf.is_some() {
            // if fin for training and validation are different file,
            // need to align and sort snv
            //let sample_val_buf =
            //    fin_sample_val.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            let dataset_val = Dataset::new(
                fin,
                gfmt,
                phe_buf.as_deref(),
                //fin_phe,
                phe_name,
                cov_name,
                snv_buf.as_deref(),
                sample_val_buf.as_deref(),
                //fin_snv,
                //fin_sample_val,
                use_missing,
                None,
                Some(dataset.snvs()),
            );
            //dataset_ext_val = Some(dataset_val);
            Some(dataset_val)
        } else {
            None
        };
        (dataset, dataset_val)
    }
    //let dataset_ext_val = dataset_ext_val;
}

/*
enum WgtArgFileType {
    Dir(PathBuf, Vec<Option<f64>>),
    File(PathBuf),
}
 */

pub fn run_boosting_score_cv(
    dout_score: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    fin_phe: Option<&Path>,
    //fin_phe: &Path,
    phe_name: Option<&str>,
    cov_name: Option<&str>,
    is_every_para: bool,
    iterations_in: Option<&[usize]>,
    dout_wgt: Option<&Path>,
    fout_wgt: Option<&Path>,
    fin_sample: Option<&Path>,
    //boost_param: BoostParam,
    learning_rates: &[f64],
    use_iter: bool,
    cross_vali: Option<usize>,
) {
    match cross_vali {
        Some(cvn) => {
            if fout_wgt.is_some() {
                panic!("Use --dout_wgt for cross-validation");
            }

            let phe_buf = fin_phe.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            let sample_buf =
                fin_sample.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());

            for cvi in 0..cvn {
                let dout_wgt_cv = wgt_boost::io::get_dir_cv(dout_wgt.unwrap(), cvi);
                let dout_score_cv = boosting_score::io::get_dir_score_cv(dout_score, cvi);
                //let dout_score_cv = wgt_boost::io::get_dir_score_cv(dout_score, cvi);

                run_boosting_score(
                    &dout_score_cv,
                    &fin,
                    gfmt,
                    phe_buf.as_deref(),
                    //fin_phe.as_deref(),
                    phe_name.as_deref(),
                    cov_name,
                    is_every_para,
                    iterations_in,
                    Some(&dout_wgt_cv), // use enum?
                    None,
                    //fin_cov.as_deref(),
                    sample_buf.as_deref(),
                    //fin_sample.as_deref(),
                    //boost_param,
                    &learning_rates,
                    use_iter,
                );
            }
        }
        None => {
            let phe_buf = fin_phe.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
            let sample_buf =
                fin_sample.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());

            run_boosting_score(
                &dout_score,
                &fin,
                gfmt,
                phe_buf.as_deref(),
                //fin_phe.as_deref(),
                phe_name.as_deref(),
                cov_name,
                is_every_para,
                iterations_in,
                dout_wgt.as_deref(), // use enum?
                fout_wgt.as_deref(),
                //fin_cov.as_deref(),
                sample_buf.as_deref(),
                //fin_sample.as_deref(),
                //boost_param,
                &learning_rates,
                use_iter,
            );
        }
    }
}

/// if SNV in fin_wgt is not in fin, then
pub fn run_boosting_score(
    dout_score: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    phe_buf: Option<&[u8]>,
    //fin_phe: Option<&Path>,
    //fin_phe: &Path,
    phe_name: Option<&str>,
    cov_name: Option<&str>,
    is_every_para: bool,
    iterations_in: Option<&[usize]>,
    dout_wgt: Option<&Path>,
    fout_wgt: Option<&Path>,
    sample_buf: Option<&[u8]>,
    //fin_sample: Option<&Path>,
    //boost_param: BoostParam,
    learning_rates: &[f64],
    use_iter: bool,
) {
    io_genot::check_valid_fin(fin, gfmt);

    //let phe_buf = fin_phe.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
    //let sample_buf = fin_sample.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());

    if !is_every_para {
        let dout_wgt = dout_wgt.unwrap();

        let file_wgt_para = wgt_boost::io::get_file_wgt_best(dout_wgt);
        let exist_fwgt = textfile::exist_file(&file_wgt_para);
        if !exist_fwgt {
            log::info!("fwgt does not exist: {:?}.", &file_wgt_para);
            return;
        }

        boosting_score::run_boosting_score_para_best(
            &dout_score,
            fin,
            gfmt,
            phe_buf.as_deref(),
            //fin_phe,
            phe_name,
            cov_name,
            &file_wgt_para,
            sample_buf.as_deref(),
            //fin_sample,
            //boost_param,
        );
    } else {
        if let Some(dout_wgt) = dout_wgt {
            // TODO: load dataset before lrs like classic score
            for learning_rate in learning_rates.iter() {
                log::debug!("learning rate: {:?}", learning_rate);

                let dout_score_lr = wgt_boost::io::get_dir_score(dout_score, *learning_rate);
                log::debug!("output to: {:?}", &dout_score_lr);

                let file_wgt_para = wgt_boost::io::get_file_wgt(dout_wgt, *learning_rate);
                let exist_fwgt = textfile::exist_file(&file_wgt_para);
                if !exist_fwgt {
                    log::info!("fwgt does not exist: {:?}.", &file_wgt_para);
                    continue;
                }

                boosting_score::run_boosting_score_para(
                    &dout_score_lr,
                    fin,
                    gfmt,
                    phe_buf.as_deref(),
                    //fin_phe,
                    phe_name,
                    cov_name,
                    iterations_in.unwrap(),
                    &file_wgt_para,
                    sample_buf.as_deref(),
                    //fin_sample,
                    //boost_param,
                    use_iter,
                );
                //let dout_wgt_para = wgt_boost::io::get_dname_para(dout_wgt, learning_rate);
                /*
                // check fwgt exist.
                wgt_boost::io::check_file_wgt_exist(&file_wgt_para);
                //wgt_boost::io::check_file_wgt_exist_dir(&dout_wgt_para);

                // TODO: check gfmt of wgt

                let iterations = wgt_boost::io::valid_iterations(iterations_in, &file_wgt_para);
                //let iterations = wgt_boost::io::valid_iterations_dir(iterations_in, &dout_wgt_para);
                log::debug!("valid iters {:?}", iterations);

                // FIXME
                // assume a1 and a2 could be flipped
                // -> DO NOT flip alpha, rather flip prediction

                // input cov or not
                // FIXME: if fin_cov is None, then has_cov=false
                let has_cov = true;
                //let has_cov = covs.is_some();

                let n_in: usize = plink::compute_num_sample(fin).unwrap();
                let (_, use_samples) = plink::make_use_samples(fin_sample, fin, n_in);
                let samples_id = plink::load_samples_id(fin, &use_samples);

                let mut wgts = WgtBoosts::new_from_file(&file_wgt_para, boost_param.boost_type());
                //let mut wgts = WgtBoosts::new_from_file_dir(&dout_wgt_para, boost_param.boost_type());
                let dataset = Dataset::new_score(fin, fin_sample, fin_cov, wgts.wgts_mut());

                boosting_score::boosting_score(
                    dout_score,
                    &iterations,
                    &wgts,
                    &dataset,
                    //&genotypes,
                    //&ys_bool,
                    &samples_id,
                    //n,
                    has_cov,
                    use_iter,
                ); */
            }
        } else if let Some(file_wgt) = fout_wgt {
            boosting_score::run_boosting_score_para(
                dout_score,
                fin,
                gfmt,
                phe_buf.as_deref(),
                //fin_phe,
                phe_name,
                cov_name,
                iterations_in.unwrap(),
                &file_wgt,
                sample_buf.as_deref(),
                //fin_sample,
                //boost_param,
                use_iter,
            );
        } else {
            panic!("sth wrong.")
        }
    }
}

/*
// moved to integration test ./tests/boosting.rs
#[cfg(test)]
mod tests {
    use super::*;

    use std::env;
    //use std::path::Path;

    #[test]
    fn test_run_boosting() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy1/unittest/toy1");
        //let fout = String::from("../../tests/boosting/result/toy1/tests_toy1");
        use io_rust::text;
        use std::fs;
        use std::path::Path;
        let fwgt = fout.to_owned() + ".wgt";
        if text::exist_file(&fwgt) {
            // delete file
            fs::remove_file(fwgt).unwrap_or_else(|why| {
                println!("! {:?}", why.kind());
            });
        }

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy1/toy1");
        //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

        //use std::path::Path;
        let fin_bim: String = path.to_str().unwrap().to_owned() + "/" + &fin + ".bim"; // + ".bim";
        println!("fin_bim {:?}", fin_bim);
        let path = Path::new(&fin_bim);
        println!("cano {:?}", path.canonicalize().unwrap());

        let t = 10;
        let boost_type = BoostType::Logit;

        run_boosting(
            &fout,
            &fin,
            t,
            &boost_type,
            None,
            None,
            None,
            None,
            false,
            None,
        );

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_withregcov() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy1/unittest/toy1_regcov");
        //let fout = String::from("../../tests/boosting/result/toy1/tests_regcov_toy1");
        use io_rust::text;
        use std::fs;
        let fwgt = fout.to_owned() + ".wgt";
        if text::exist_file(&fwgt) {
            // delete file
            fs::remove_file(fwgt).unwrap_or_else(|why| {
                println!("! {:?}", why.kind());
            });
        }

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy1/toy1");
        //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

        let fin_wgt_cov = Some(String::from(
            "../../tests/boosting/data/toy1/toy1.regcov.wgtcov",
        ));

        let fin_cov = Some(String::from("../../tests/boosting/data/toy1/toy1.cov"));

        println!("fin_wgt_cov {:?}", fin_wgt_cov);

        let t = 10;
        let boost_type = BoostType::Logit;

        run_boosting(
            &fout,
            &fin,
            t,
            &boost_type,
            None,
            None,
            fin_cov.as_deref(),
            fin_wgt_cov.as_deref(),
            false,
            None,
        );

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_toy2_withregcov() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy2/unittest/toy2_regcov");
        //let fout = String::from("../../tests/boosting/result/toy2/tests_regcov_toy2");
        use io_rust::text;
        use std::fs;
        let fwgt = fout.to_owned() + ".wgt";
        if text::exist_file(&fwgt) {
            // delete file
            fs::remove_file(fwgt).unwrap_or_else(|why| {
                println!("! {:?}", why.kind());
            });
        }

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy2/toy2.chrom%");
        //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

        /*
        //use std::path::Path;
        let fin_bim: String = path.to_str().unwrap().to_owned() + "/" + &fin + ".bim"; // + ".bim";
        println!("fin_bim {:?}", fin_bim);
        let path = Path::new(&fin_bim);
        println!("cano {:?}", path.canonicalize().unwrap());
        */

        let fin_wgt_cov = Some(String::from(
            "../../tests/boosting/data/toy2/toy2.regcov.wgtcov",
        ));

        let fin_cov = Some(String::from("../../tests/boosting/data/toy2/toy2.cov"));
        let fin_snv = Some(String::from("../../tests/boosting/data/toy2/toy2.snvs"));
        let fin_sample = Some(String::from(
            "../../tests/boosting/data/toy2/toy2.cv0.samples",
        ));

        println!("fin_wgt_cov {:?}", fin_wgt_cov);

        let t = 10;
        let boost_type = BoostType::Logit;

        run_boosting(
            &fout,
            &fin,
            t,
            &boost_type,
            fin_snv.as_deref(),
            fin_sample.as_deref(),
            fin_cov.as_deref(),
            fin_wgt_cov.as_deref(),
            false,
            None,
        );

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_score() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy1/unittest/toy1_regcov");
        let fin_wgt = String::from("../../tests/boosting/data/toy1/toy1_for_score");
        //let fin_wgt = String::from("../../tests/boosting/result/toy1/unittest/toy1_for_score");

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy1/toy1");
        let fin_cov = Some(String::from("../../tests/boosting/data/toy1/toy1.cov"));

        let iters = [1, 2, 10];

        run_boosting_score(&fout, &fin, &iters, &fin_wgt, fin_cov.as_deref(), None);

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_score_withregcov() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy1/unittest/toy1");
        let fin_wgt = String::from("../../tests/boosting/data/toy1/toy1_regcov_for_score");
        //let fin_wgt = String::from("../../tests/boosting/result/toy1/unittest/toy1_for_score");

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy1/toy1");
        //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

        let iters = [1, 2, 10];

        run_boosting_score(&fout, &fin, &iters, &fin_wgt, None, None);

        println!("Done!!");
    }

    #[test]
    fn test_run_boosting_score_toy2_withregcov() {
        let path = env::current_dir().unwrap();
        println!("current dir: {}", path.display());

        let fout = String::from("../../tests/boosting/result/toy2/unittest/toy2");
        //let fin_wgt = String::from("../../tests/boosting/result/toy2/unittest/toy2_regcov_for_score");
        let fin_wgt = String::from("../../tests/boosting/data/toy2/toy2_regcov_for_score");

        // tmp: better to put the file in this crate
        let fin = String::from("../../tests/boosting/data/toy2/toy2.chrom%");
        let fin_cov = Some(String::from("../../tests/boosting/data/toy2/toy2.cov"));
        let fin_sample = Some(String::from(
            "../../tests/boosting/data/toy2/toy2.cv0.samples",
        ));

        let iters = [1, 2, 10];

        run_boosting_score(
            &fout,
            &fin,
            &iters,
            &fin_wgt,
            fin_cov.as_deref(),
            fin_sample.as_deref(),
        );

        println!("Done!!");
    }
}
*/
