//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
//

pub mod boosting; // pub for bench
mod boosting_param;
mod boosting_score;
mod wgt_boost;
mod wgt_boosts;

use std::path::Path;

use boosting::{loss::LossStruct, ContingencyTable};
pub use boosting_param::{BoostMethod, BoostParam, BoostType, Eps};
use genetics::alloc;
use genetics::{plink, Dataset};
use wgt_boost::WgtBoost;
use wgt_boosts::WgtBoosts;

/// Run GenoBoost. Input plink file.
pub fn run_boosting(
    dout: &Path,
    //dout: &str,
    fin: &Path,
    iteration: usize,
    boost_method: BoostMethod,
    boost_param: BoostParam,
    //boost_type: &BoostType,
    fin_snv: Option<&Path>,
    fin_sample: Option<&Path>,
    fin_cov: Option<&Path>,
    //fin_wgt_cov: Option<&str>,
    //is_dom_rec: bool,
    is_resume: bool,
    is_write_loss: bool,
    //dloss: Option<&Path>,
    prune_snv: Option<f64>,
    learning_rates: &[Option<f64>],
) {
    // check fwgt does not exist.
    if !is_resume {
        for learning_rate in learning_rates.iter() {
            wgt_boost::io::check_file_wgt_not_exist(dout, learning_rate);
        }
    }

    plink::check_valid_fin(fin);

    let use_missing = boost_param.boost_type().use_missing();
    let dataset = Dataset::new(fin, fin_snv, fin_sample, fin_cov, use_missing);

    // extract snvs by loss function
    let dataset_ext;
    if let Some(prop_prune_snv) = prune_snv {
        //TODO: mv boostmethod::classic
        let n = dataset.samples().samples_n();
        let mut scores: Vec<f64> = vec![0.0; n];
        let mut ws: Vec<f64> = vec![1.0 / (n as f64); n];
        // len != n only for ps
        let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        ps.resize(n + 32, 0.0f64);
        let mut pred: Vec<u8> = vec![0u8; n];

        let mut wgts = WgtBoosts::new(boost_param.boost_type());
        let _ = boosting::boosting_covs(
            &mut wgts,
            boost_param,
            &dataset,
            &mut scores,
            &mut ws,
            &mut ps,
            &mut pred,
        );

        let mut loss = LossStruct::new(boost_param.boost_type(), dataset.snvs().snvs_n());
        boosting::loss::calculate_loss_gt(
            &mut loss,
            dataset.genot(),
            &mut ps,
            dataset.samples().phe(),
            boost_param,
        );

        let use_snvs_loss = loss.search_topprop(prop_prune_snv);

        dataset_ext = dataset.extract_snvs(&use_snvs_loss);
    } else {
        dataset_ext = dataset;
    }

    wgt_boost::io::create_dir(&dout);

    //let (file, _) = wgt_boost::io::bufwriter_file_wgt(fout);
    //let mut writer = wgt_boost::io::bufwriter_fwgt_append(dout);
    //let mut writer = wgt_boost::io::bufwriter_fwgt(fout);
    //wgt_boost::io::write_cols(&mut writer, boost_param.boost_type());

    //let mut writer_loss = wgt_boost::io::bufwriter_floss(fout);

    //let dout_option = if is_resume || is_write_loss {
    //    //Some(wgt_boost::io::get_fname_wgt_pathbuf(dout))
    //    Some(dout)
    //} else {
    //    None
    //};

    for learning_rate in learning_rates.iter() {
        println!("learning rate: {:?}", learning_rate);
        let dout_para = wgt_boost::io::get_dname_para(dout, learning_rate);
        wgt_boost::io::create_dir(&dout_para);
        let mut writer = wgt_boost::io::bufwriter_fwgt_append(&dout_para);
        match boost_method {
            BoostMethod::Classic => {
                println!("Run boosting");
                boosting::boosting(
                    &mut writer,
                    iteration,
                    boost_param,
                    &dataset_ext,
                    is_resume,
                    is_write_loss,
                    Some(&dout_para),
                )
            }
            _ => unimplemented!(),
        }
    }

    /*
    match boost_method {
        BoostMethod::Classic => {
            println!("Run boosting");
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
        //    println!("Run boosting pruning");
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
          //    println!("Run boosting ss.");
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

pub fn run_boosting_score_para(
    dout_score: &Path,
    fin: &Path,
    iterations_in: &[usize],
    file_wgt: &Path,
    fin_cov: Option<&Path>,
    fin_sample: Option<&Path>,
    boost_param: BoostParam,
    use_iter: bool,
) {
    // check fwgt exist.
    wgt_boost::io::check_file_wgt_exist(&file_wgt);
    //wgt_boost::io::check_file_wgt_exist_dir(&dout_wgt_para);

    // TODO: check format of wgt

    let iterations = wgt_boost::io::valid_iterations(iterations_in, &file_wgt);
    //let iterations = wgt_boost::io::valid_iterations_dir(iterations_in, &dout_wgt_para);
    println!("valid iters {:?}", iterations);

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

    let mut wgts = WgtBoosts::new_from_file(&file_wgt, boost_param.boost_type());
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
    iterations_in: &[usize],
    dout_wgt: Option<&Path>,
    fout_wgt: Option<&Path>,
    fin_cov: Option<&Path>,
    fin_sample: Option<&Path>,
    boost_param: BoostParam,
    learning_rates: &[Option<f64>],
    use_iter: bool,
) {
    plink::check_valid_fin(fin);

    if let Some(dout_wgt) = dout_wgt {
        for learning_rate in learning_rates.iter() {
            println!("learning rate: {:?}", learning_rate);

            let file_wgt_para = wgt_boost::io::get_file_wgt(dout_wgt, learning_rate);
            run_boosting_score_para(
                dout_score,
                fin,
                iterations_in,
                &file_wgt_para,
                fin_cov,
                fin_sample,
                boost_param,
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
            println!("valid iters {:?}", iterations);

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
            iterations_in,
            &file_wgt,
            fin_cov,
            fin_sample,
            boost_param,
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
