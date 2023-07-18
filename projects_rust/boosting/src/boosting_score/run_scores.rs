use crossbeam;
use std::path::Path;

use super::score;
use crate::wgt_boosts::WgtBoosts;
use genetics::Dataset;

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

    if use_iter & has_cov {
        log::debug!("Start iters with covariates.");
        score::calculate_write_score_iterations(
            dout, iterations, wgts, dataset, samples_id, false,
        );
    }
}
*/

// use crossbeam
// Probably, use crossbeam here is better than in adding score
// but not tried yet
// TODO: compare speed with crossbeam or not
//pub fn boosting_score_crossbeam(
pub fn boosting_score(
    dout: &Path,
    iterations: &[usize],
    wgts: &WgtBoosts,
    dataset: &Dataset,
    // why not samples_id in dataset.samples?
    //samples_id: &[(String, String)],
    samples_id: &[String],
    has_cov: bool,
    use_iter: bool,
) {
    log::debug!("score for iterations");

    //https://rust-lang-nursery.github.io/rust-cookbook/concurrency/threads.html
    crossbeam::scope(|s| {
        let thread0 = s.spawn(|_| {
            if use_iter & has_cov {
                log::debug!("Start iters with covariates.");
                score::calculate_write_score_iterations(
                    dout, iterations, wgts, dataset, samples_id, false,
                );
            }
        });

        let thread1 = s.spawn(|_| {
            if has_cov {
                log::debug!("Start nsnvs with covariates.");
                score::calculate_write_score_nsnvs(
                    dout, iterations, wgts, dataset, samples_id, false,
                );
            }
        });

        let thread2 = s.spawn(|_| {
            if use_iter {
                log::debug!("Start iters without covariates.");
                score::calculate_write_score_iterations(
                    dout, iterations, wgts, dataset, samples_id, true,
                );
            }
        });

        let thread3 = s.spawn(|_| {
            log::debug!("Start nsnvs without covariates.");
            score::calculate_write_score_nsnvs(dout, iterations, wgts, dataset, samples_id, true);
        });

        thread0.join().unwrap();
        thread1.join().unwrap();
        thread2.join().unwrap();
        thread3.join().unwrap();
    })
    .unwrap();
}

pub fn boosting_score_para_best(
    dout: &Path,
    wgts: &WgtBoosts,
    dataset: &Dataset,
    // why not samples_id in dataset.samples?
    //samples_id: &[(String, String)],
    samples_id: &[String],
    has_cov: bool,
) {
    log::debug!("score for iterations");

    //https://rust-lang-nursery.github.io/rust-cookbook/concurrency/threads.html
    crossbeam::scope(|s| {
        let thread1 = s.spawn(|_| {
            if has_cov {
                log::debug!("Start nsnvs with covariates.");
                score::calculate_write_score_para_best(dout, wgts, dataset, samples_id, false);
            }
        });

        let thread3 = s.spawn(|_| {
            log::debug!("Start nsnvs without covariates.");
            score::calculate_write_score_para_best(dout, wgts, dataset, samples_id, true);
        });

        thread1.join().unwrap();
        thread3.join().unwrap();
    })
    .unwrap();
}

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
