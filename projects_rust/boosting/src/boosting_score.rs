//! Application of **Genoboost**.
//! Input plink file to run Genoboost.

//pub mod io;
mod run_scores;
pub mod score;

use crate::dout_file::DoutScoreParaFile;
use crate::wgt_boost;
use crate::wgt_boosts::WgtBoosts;
use genetics::dataset_file::DatasetFile;
use genetics::{textfile, Dataset};
pub use run_scores::*;
use std::path::Path;

pub fn run_boosting_score_para_best(
    dout_score: &DoutScoreParaFile,
    dfile: &DatasetFile,
    file_wgt: &Path,
    fill_missing_in_dataset: bool,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
    mem: Option<usize>,
) {
    // check fwgt exist.
    textfile::check_exist_file(&file_wgt);

    // input cov or not
    let has_cov = dfile.cov_name().is_some();

    let mut wgts = WgtBoosts::new_from_file(&file_wgt);
    //let mut wgts = WgtBoosts::new_from_file(&file_wgt, boost_param.boost_type());

    let dataset = Dataset::new_datasetfile_score(
        dfile,
        wgts.wgts_mut(),
        fill_missing_in_dataset,
        allow_nonexist_snv,
        use_snv_pos,
        mem,
    );

    boosting_score_para_best(
        dout_score,
        &wgts,
        &dataset,
        has_cov,
        allow_nonexist_snv,
        missing_to_mode,
        missing_to_mean,
    );
}

pub fn run_boosting_score_para(
    dout_score: &DoutScoreParaFile,
    dfile: &DatasetFile,
    iterations_in: &[usize],
    file_wgt: &Path,
    use_iter: bool,
    fill_missing_in_dataset: bool,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
    mem: Option<usize>,
) {
    // check fwgt exist.
    textfile::check_exist_file(&file_wgt);

    let iterations = wgt_boost::io::valid_iterations(iterations_in, &file_wgt);
    log::debug!("valid iters {:?}", iterations);

    // input cov or not
    let has_cov = dfile.cov_name().is_some();

    let mut wgts = WgtBoosts::new_from_file(&file_wgt);

    // in boosting, misisng -> mode using maf info
    let dataset = Dataset::new_datasetfile_score(
        dfile,
        wgts.wgts_mut(),
        fill_missing_in_dataset,
        allow_nonexist_snv,
        use_snv_pos,
        mem,
    );

    boosting_score(
        dout_score,
        &iterations,
        &wgts,
        &dataset,
        has_cov,
        use_iter,
        allow_nonexist_snv,
        missing_to_mode,
        missing_to_mean,
    );
}
