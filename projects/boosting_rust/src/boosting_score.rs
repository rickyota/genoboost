//! Application of **Genoboost**.
//! Input plink file to run Genoboost.

mod io;
pub mod score;

use std::path::Path;

use crossbeam;
use genetics::Dataset;

use super::WgtBoost;
use crate::wgt_boosts::WgtBoosts;

// use crossbeam
// Probably, use multithreads here is better than in adding score
pub fn boosting_score(
    dout: &Path,
    iterations: &[usize],
    wgts: &WgtBoosts,
    dataset: &Dataset,
    samples_id: &[(String, String)],
    has_cov: bool,
    use_iter: bool,
) {
    println!("score for iterations");

    //https://rust-lang-nursery.github.io/rust-cookbook/concurrency/threads.html

    crossbeam::scope(|s| {
        let thread0 = s.spawn(|_| {
            if use_iter & has_cov {
                println!("Start iters with covariates.");
                score::calculate_write_score_iterations(
                    dout, iterations, wgts, dataset, samples_id, false,
                );
            }
        });

        let thread1 = s.spawn(|_| {
            if has_cov {
                println!("Start nsnvs with covariates.");
                score::calculate_write_score_nsnvs(
                    dout, iterations, wgts, dataset, samples_id, false,
                );
            }
        });

        let thread2 = s.spawn(|_| {
            if use_iter {
                println!("Start iters without covariates.");
                score::calculate_write_score_iterations(
                    dout, iterations, wgts, dataset, samples_id, true,
                );
            }
        });

        let thread3 = s.spawn(|_| {
            println!("Start nsnvs without covariates.");
            score::calculate_write_score_nsnvs(dout, iterations, wgts, dataset, samples_id, true);
        });

        thread0.join().unwrap();
        thread1.join().unwrap();
        thread2.join().unwrap();
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
