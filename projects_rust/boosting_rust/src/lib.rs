//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
//

mod boosting_param;
mod boosting_score;
pub mod boosting_train; // pub for bench
mod wgt_boost;
mod wgt_boosts;

use std::collections::HashSet;
use std::{path::Path, time::Instant};

pub use boosting_param::{BoostMethod, BoostParam, BoostType, EffEps, Eps, IterationNumber};
use boosting_train::{loss::LossStruct, ContingencyTable};
use genetics::text;
use genetics::{plink, vec, Dataset, Snvs};
use wgt_boost::WgtBoost;
use wgt_boosts::WgtBoosts;

use crate::boosting_train::sample_weight::SampleWeight;

#[macro_use]
extern crate assert_float_eq;

/// Run GenoBoost. Input plink file.
pub fn run_boosting(
    dout: &Path,
    fin: &Path,
    fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    boost_method: BoostMethod,
    boost_param: BoostParam,
    fin_snv: Option<&Path>,
    fin_sample: Option<&Path>,
    fin_cov: Option<&Path>,
    fin_sample_val: Option<&Path>,
    //is_dom_rec: bool,
    use_adjloss: bool,
    use_const_for_loss: bool,
    is_resume: bool,
    is_write_loss: bool,
    prune_snv: Option<f64>,
    learning_rates: &[Option<f64>],
) {
    // check fwgt does not exist.
    if !is_resume {
        for learning_rate in learning_rates.iter() {
            wgt_boost::io::check_file_wgt_not_exist(dout, learning_rate);
        }
    }

    // TODO: if all wgt exceeds #SNVs, then exit here.
    plink::check_valid_fin(fin);

    let (dataset, dataset_val) = load_dataset_boosting(
        dout,
        fin,
        fin_phe,
        phe_name,
        boost_param,
        fin_snv,
        fin_sample,
        fin_cov,
        fin_sample_val,
        use_adjloss,
        prune_snv,
    );

    wgt_boost::io::create_dir(&dout);

    for learning_rate in learning_rates.iter() {
        log::debug!("learning rate: {:?}", learning_rate);
        let boost_param_lr = boost_param.set_learning_rate(learning_rate.clone());
        let dout_para = wgt_boost::io::get_dname_para(dout, learning_rate);
        wgt_boost::io::create_dir(&dout_para);
        let mut writer = wgt_boost::io::bufwriter_fwgt_append(&dout_para);
        if boost_param_lr.batch_way().is_none() {
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
                    )
                }
                _ => unimplemented!(),
            }
        }
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

fn load_dataset_boosting(
    dout: &Path,
    fin: &Path,
    fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    boost_param: BoostParam,
    fin_snv: Option<&Path>,
    fin_sample: Option<&Path>,
    fin_cov: Option<&Path>,
    fin_sample_val: Option<&Path>,
    use_adjloss: bool,
    prune_snv: Option<f64>,
) -> (Dataset, Option<Dataset>) {
    let use_missing = boost_param.boost_type().use_missing();
    // create dataset
    // extract snvs by loss function
    if let Some(prop_prune_snv) = prune_snv {
        log::info!("Prune SNVs by decreasing loss: {}", prop_prune_snv);
        let start = Instant::now();

        fin_sample_val.expect("Not Implemented");

        // TODO: better
        // to get m
        let m;
        {
            let m_in: usize = plink::compute_num_snv(fin).unwrap();
            let snvs_in = plink::load_snvs(fin, m_in);
            (m, _) = plink::make_use_snvs(fin_snv, &snvs_in);
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
                fin_phe,
                phe_name,
                fin_snv,
                fin_sample,
                fin_cov,
                use_missing,
                Some(&filt_snv),
                None,
            );

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
        let snvs = Snvs::new_plink(fin);
        loss.write_writer(&mut writer_loss, &snvs);

        let use_snvs_loss = loss.search_topprop(prop_prune_snv);

        log::debug!("created use_snvs_loss {}", use_snvs_loss.len());

        let dataset = Dataset::new(
            fin,
            fin_phe,
            phe_name,
            fin_snv,
            fin_sample,
            fin_cov,
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
        let dataset = Dataset::new(
            fin,
            fin_phe,
            phe_name,
            fin_snv,
            fin_sample,
            fin_cov,
            use_missing,
            None,
            None,
        );
        //dataset_ext = dataset;

        let dataset_val=if fin_sample_val.is_some() {
            // if fin for training and validation are different file,
            // need to align and sort snv
            let dataset_val = Dataset::new(
                fin,
                fin_phe,
                phe_name,
                fin_snv,
                fin_sample_val,
                fin_cov,
                use_missing,
                None,
                Some(dataset.snvs()),
            );
            //dataset_ext_val = Some(dataset_val);
            Some(dataset_val)
        }else{
            None
        };
        (dataset, dataset_val)
    }
    //let dataset_ext_val = dataset_ext_val;
}

pub fn run_boosting_score_para(
    dout_score: &Path,
    fin: &Path,
    fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    iterations_in: &[usize],
    file_wgt: &Path,
    fin_cov: Option<&Path>,
    fin_sample: Option<&Path>,
    //boost_param: BoostParam,
    use_iter: bool,
) {
    // check fwgt exist.
    wgt_boost::io::check_file_wgt_exist(&file_wgt);
    //wgt_boost::io::check_file_wgt_exist_dir(&dout_wgt_para);

    // TODO: check format of wgt

    let iterations = wgt_boost::io::valid_iterations(iterations_in, &file_wgt);
    //let iterations = wgt_boost::io::valid_iterations_dir(iterations_in, &dout_wgt_para);
    log::debug!("valid iters {:?}", iterations);

    // input cov or not
    let has_cov = fin_cov.is_some();
    //let has_cov = true;

    let n_in: usize = plink::compute_num_sample(fin).unwrap();
    let (_, use_samples) = plink::make_use_samples(fin_sample, fin, n_in);
    let samples_id = plink::load_samples_id(fin, &use_samples);

    let mut wgts = WgtBoosts::new_from_file(&file_wgt);
    //let mut wgts = WgtBoosts::new_from_file(&file_wgt, boost_param.boost_type());
    //let mut wgts = WgtBoosts::new_from_file_dir(&dout_wgt_para, boost_param.boost_type());
    let use_missing = wgts.use_missing();
    let dataset = Dataset::new_score(
        fin,
        fin_phe,
        phe_name,
        fin_sample,
        fin_cov,
        wgts.wgts_mut(),
        use_missing,
    );

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
    );
}

/*
enum WgtArgFileType {
    Dir(PathBuf, Vec<Option<f64>>),
    File(PathBuf),
}
 */

/// if SNV in fin_wgt is not in fin, then
pub fn run_boosting_score(
    dout_score: &Path,
    fin: &Path,
    fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    iterations_in: &[usize],
    dout_wgt: Option<&Path>,
    fout_wgt: Option<&Path>,
    fin_cov: Option<&Path>,
    fin_sample: Option<&Path>,
    //boost_param: BoostParam,
    learning_rates: &[Option<f64>],
    use_iter: bool,
) {
    plink::check_valid_fin(fin);

    if let Some(dout_wgt) = dout_wgt {
        // TODO: load dataset before lrs like classic score
        for learning_rate in learning_rates.iter() {
            log::debug!("learning rate: {:?}", learning_rate);

            let dout_score_lr = wgt_boost::io::get_dir_score(dout_score, learning_rate);
            log::debug!("output to: {:?}", &dout_score_lr);

            let file_wgt_para = wgt_boost::io::get_file_wgt(dout_wgt, learning_rate);
            let exist_fwgt = text::exist_file(&file_wgt_para);
            if !exist_fwgt {
                log::info!("fwgt does not exist: {:?}.", &file_wgt_para);
                continue;
            }

            run_boosting_score_para(
                &dout_score_lr,
                fin,
                fin_phe,
                phe_name,
                iterations_in,
                &file_wgt_para,
                fin_cov,
                fin_sample,
                //boost_param,
                use_iter,
            );
            //let dout_wgt_para = wgt_boost::io::get_dname_para(dout_wgt, learning_rate);
            /*
            // check fwgt exist.
            wgt_boost::io::check_file_wgt_exist(&file_wgt_para);
            //wgt_boost::io::check_file_wgt_exist_dir(&dout_wgt_para);

            // TODO: check format of wgt

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
        run_boosting_score_para(
            dout_score,
            fin,
            fin_phe,
            phe_name,
            iterations_in,
            &file_wgt,
            fin_cov,
            fin_sample,
            //boost_param,
            use_iter,
        );
    } else {
        panic!("sth wrong.")
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
