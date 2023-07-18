//! Application of **Genoboost**.
//! Input plink file to run Genoboost.

pub mod io;
mod run_scores;
pub mod score;

pub use crate::boosting_param::{BoostMethod, BoostParam, BoostType, EffEps, Eps, IterationNumber};
use crate::wgt_boost;
use crate::wgt_boosts::WgtBoosts;
use genetics::sample;
use genetics::GenotFormat;
use genetics::{io_genot, Dataset};
pub use run_scores::*;
use std::path::Path;

pub fn run_boosting_score_para_best(
    dout_score: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    phe_buf: Option<&[u8]>,
    //fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    cov_name: Option<&str>,
    file_wgt: &Path,
    extract_sample_buf: Option<&[u8]>,
    //fin_sample: Option<&Path>,
    //boost_param: BoostParam,
) {

    // check fwgt exist.
    wgt_boost::io::check_file_wgt_exist(&file_wgt);

    // TODO: check gfmt of wgt

    //let iterations = wgt_boost::io::valid_iterations(iterations_in, &file_wgt);
    ////let iterations = wgt_boost::io::valid_iterations_dir(iterations_in, &dout_wgt_para);
    //log::debug!("valid iters {:?}", iterations);

    // input cov or not
    let has_cov = cov_name.is_some();
    //let has_cov = fin_cov.is_some();
    //let has_cov = true;

    //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
    let (_, use_samples) = sample::make_use_samples_buf(extract_sample_buf, fin, gfmt);
    let samples_id = io_genot::load_samples_id(fin, gfmt, Some(&use_samples));

    let mut wgts = WgtBoosts::new_from_file(&file_wgt);
    //let mut wgts = WgtBoosts::new_from_file(&file_wgt, boost_param.boost_type());
    //let mut wgts = WgtBoosts::new_from_file_dir(&dout_wgt_para, boost_param.boost_type());
    let use_missing = wgts.use_missing();
    let dataset = Dataset::new_score(
        fin,
        gfmt,
        phe_buf,
        phe_name,
        cov_name,
        extract_sample_buf.as_deref(),
        wgts.wgts_mut(),
        use_missing,
    );

    boosting_score_para_best(
        dout_score,
        &wgts,
        &dataset,
        //&genotypes,
        //&ys_bool,
        &samples_id,
        //n,
        has_cov,
    ); 

    
}

pub fn run_boosting_score_para(
    dout_score: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    phe_buf: Option<&[u8]>,
    //fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    // Option<> for nocov only
    cov_name: Option<&str>,
    iterations_in: &[usize],
    file_wgt: &Path,
    extract_sample_buf: Option<&[u8]>,
    //fin_sample: Option<&Path>,
    //boost_param: BoostParam,
    use_iter: bool,
) {
    // check fwgt exist.
    wgt_boost::io::check_file_wgt_exist(&file_wgt);
    //wgt_boost::io::check_file_wgt_exist_dir(&dout_wgt_para);

    // TODO: check gfmt of wgt

    let iterations = wgt_boost::io::valid_iterations(iterations_in, &file_wgt);
    //let iterations = wgt_boost::io::valid_iterations_dir(iterations_in, &dout_wgt_para);
    log::debug!("valid iters {:?}", iterations);

    // input cov or not
    let has_cov = cov_name.is_some();
    //let has_cov = fin_cov.is_some();
    //let has_cov = true;

    //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
    let (_, use_samples) = sample::make_use_samples_buf(extract_sample_buf, fin, gfmt);
    //let (_, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
    let samples_id = io_genot::load_samples_id(fin, gfmt, Some(&use_samples));

    let mut wgts = WgtBoosts::new_from_file(&file_wgt);
    //let mut wgts = WgtBoosts::new_from_file(&file_wgt, boost_param.boost_type());
    //let mut wgts = WgtBoosts::new_from_file_dir(&dout_wgt_para, boost_param.boost_type());
    let use_missing = wgts.use_missing();
    let dataset = Dataset::new_score(
        fin,
        gfmt,
        phe_buf.as_deref(),
        //fin_phe,
        phe_name,
        cov_name,
        extract_sample_buf.as_deref(),
        //fin_sample,
        wgts.wgts_mut(),
        use_missing,
    );

    boosting_score(
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
