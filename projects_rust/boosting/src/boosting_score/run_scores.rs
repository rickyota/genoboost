//use crossbeam;

use super::score;
use crate::{dout_file::DoutScoreParaFile, wgt_boosts::WgtBoosts};
use genetics::Dataset;

/// not use crossbeam
/// since using rayon inside
pub fn boosting_score(
    dout: &DoutScoreParaFile,
    iterations: &[usize],
    wgts: &WgtBoosts,
    dataset: &Dataset,
    has_cov: bool,
    use_iter: bool,
    allow_nonexist_snv: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
) {
    log::debug!("score for iterations");

    if use_iter && has_cov {
        log::debug!("Start iters with covariates.");
        score::calculate_write_score_iterations(
            dout,
            iterations,
            wgts,
            dataset,
            false,
            allow_nonexist_snv,
            missing_to_mode,
            missing_to_mean,
        );
    }

    if has_cov {
        log::debug!("Start nsnvs with covariates.");
        score::calculate_write_score_nsnvs(
            dout,
            iterations,
            wgts,
            dataset,
            false,
            allow_nonexist_snv,
            missing_to_mode,
            missing_to_mean,
        );
    }

    if use_iter {
        log::debug!("Start iters without covariates.");
        score::calculate_write_score_iterations(
            dout,
            iterations,
            wgts,
            dataset,
            true,
            allow_nonexist_snv,
            missing_to_mode,
            missing_to_mean,
        );
    }

    log::debug!("Start nsnvs without covariates.");
    score::calculate_write_score_nsnvs(
        dout,
        iterations,
        wgts,
        dataset,
        true,
        allow_nonexist_snv,
        missing_to_mode,
        missing_to_mean,
    );
}

pub fn boosting_score_para_best(
    dout: &DoutScoreParaFile,
    wgts: &WgtBoosts,
    dataset: &Dataset,
    has_cov: bool,
    allow_nonexist_snv: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
) {
    log::debug!("score for iterations");

    if has_cov {
        log::debug!("Start nsnvs with covariates.");
        score::calculate_write_score_para_best(
            dout,
            wgts,
            dataset,
            false,
            allow_nonexist_snv,
            missing_to_mode,
            missing_to_mean,
        );
    }

    log::debug!("Start nsnvs without covariates.");
    score::calculate_write_score_para_best(
        dout,
        wgts,
        dataset,
        true,
        allow_nonexist_snv,
        missing_to_mode,
        missing_to_mean,
    );
}

/*
// test for crossbeam
pub fn boosting_score(
    dout: &Path,
    iterations: &[usize],
    wgts: &WgtBoosts,
    dataset: &Dataset,
    samples_id: &[(String, String)],
    has_cov: bool,
    use_iter: bool,
) {
    log::debug!("score for iterations");


    log::debug!("Start nsnvs without covariates.");
    score::calculate_write_score_nsnvs(dout, iterations, wgts, dataset, samples_id, true);


    if has_cov {
        log::debug!("Start nsnvs with covariates.");
        score::calculate_write_score_nsnvs(
            dout, iterations, wgts, dataset, samples_id, false,
        );
    }

    if use_iter {
        log::debug!("Start iters without covariates.");
        score::calculate_write_score_iterations(
            dout, iterations, wgts, dataset, samples_id, true,
        );
    }

    if use_iter && has_cov {
        log::debug!("Start iters with covariates.");
        score::calculate_write_score_iterations(
            dout, iterations, wgts, dataset, samples_id, false,
        );
    }
}
*/

//// use crossbeam
//// Probably, use crossbeam here is better than in adding score
//// but not tried yet
//// compare speed with crossbeam or not
////pub fn boosting_score_crossbeam(
//pub fn boosting_score(
//    dout: &DoutScoreParaFile,
//    //dout: &Path,
//    iterations: &[usize],
//    wgts: &WgtBoosts,
//    dataset: &Dataset,
//    // why not samples_id in dataset.samples?
//    //samples_id: &[(String, String)],
//    //samples_id: &[String],
//    has_cov: bool,
//    use_iter: bool,
//) {
//    log::debug!("score for iterations");
//
//    //https://rust-lang-nursery.github.io/rust-cookbook/concurrency/threads.html
//    crossbeam::scope(|s| {
//        let thread0 = s.spawn(|_| {
//            if use_iter && has_cov {
//                log::debug!("Start iters with covariates.");
//                score::calculate_write_score_iterations(
//                    dout, iterations, wgts, dataset,
//                    false,
//                    //dout, iterations, wgts, dataset, samples_id, false,
//                );
//            }
//        });
//
//        let thread1 = s.spawn(|_| {
//            if has_cov {
//                log::debug!("Start nsnvs with covariates.");
//                score::calculate_write_score_nsnvs(
//                    dout, iterations, wgts, dataset,
//                    false,
//                    //dout, iterations, wgts, dataset, samples_id, false,
//                );
//            }
//        });
//
//        let thread2 = s.spawn(|_| {
//            if use_iter {
//                log::debug!("Start iters without covariates.");
//                score::calculate_write_score_iterations(dout, iterations, wgts, dataset, true);
//            }
//        });
//
//        let thread3 = s.spawn(|_| {
//            log::debug!("Start nsnvs without covariates.");
//            score::calculate_write_score_nsnvs(dout, iterations, wgts, dataset, true);
//            //score::calculate_write_score_nsnvs(dout, iterations, wgts, dataset, samples_id, true);
//        });
//
//        thread0.join().unwrap();
//        thread1.join().unwrap();
//        thread2.join().unwrap();
//        thread3.join().unwrap();
//    })
//    .unwrap();
//}
