//!
//! use test bench
//! nightly
//!
//! criterion seems better
//!
//! commeted out since this raise error on non-nightly env.
//!

/*
#![feature(test)]

//use test; // error due to old file system in test
extern crate test;

use boosting_rust::{boosting::loss::calculate, BoostParam};

// 8 x 1 bit
//type B8 = u8;

mod common;

/// simd is ~5x faster
/// ~29ms
#[bench]
fn bench_calculate_loss_gt_simd(b: &mut test::Bencher) {
    let (dataset, ps, mut losss) = common::setup_test_n100k();

    /*
    let num_cpus = 1;
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_cpus)
        .build_global()
        .unwrap();
        */
    println!("num_thread set: {}", rayon::current_num_threads());

    b.iter(|| {
        calculate::calculate_loss_gt_simd(
            &mut losss,
            &dataset.genot(),
            &ps,
            &dataset.samples().ys(),
            dataset.samples().samples_n(),
            BoostParam::new_type1(),
        )
    })
}

/// ~160ms
#[bench]
fn bench_calculate_loss_gt_nosimd(b: &mut test::Bencher) {
    let (dataset, ps, mut losss) = common::setup_test_n100k();
    //let (predictions, ys, _, n, ps, mut losss) = common::setup_test_n100k();

    /*
    let num_cpus = 1;
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_cpus)
        .build_global()
        .unwrap();
         */
    println!("num_thread set: {}", rayon::current_num_threads());

    b.iter(|| {
        calculate::calculate_loss_gt_nosimd(
            &mut losss,
            &dataset.genot(),
            &ps,
            &dataset.samples().ys(),
            dataset.samples().samples_n(),
            BoostParam::new_type1(),
        )
    })
}

 */
