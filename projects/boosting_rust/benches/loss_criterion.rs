//! for criterion
//! output is in `./target/criterion/report`
//! opt-level is 3 on `cargo bench`
//!
//!
//! https://bheisler.github.io/criterion.rs/book/user_guide/benchmarking_with_inputs.html
//!
//!
// TODO: use $[criterion] to avoid group etc. https://bheisler.github.io/criterion.rs/book/user_guide/custom_test_framework.html
// also add [[bench]] in Cargo.toml

#[macro_use]
extern crate criterion;

//use criterion::black_box;
use criterion::{BenchmarkId, Criterion};

use boosting_rust::{boosting::loss::calc, BoostParam};

mod common;
//mod set_count_criterion;

// 8 x 1 bit
//type B8 = u8;

/*
fn bench_calculate_loss_gt_simd(c: &mut Criterion) {
    let (predictions, ys, _, n, ps, mut losss) = common::setup_test();

    println!("num_thread set: {}", rayon::current_num_threads());

    c.bench_function("simd", |b| {
        b.iter(|| {
            calculate::calculate_loss_gt_simd(
                black_box(&mut losss),
                black_box(&predictions),
                black_box(&ps),
                black_box(&ys),
                black_box(n),
            )
        })
        //b.iter(|| calculate::calculate_loss_gt_simd(&mut losss, &predictions, &ps, &ys, n))
    });
}
*/

/*
fn bench_calculate_loss_gt_nosimd(c: &mut Criterion) {
    let (predictions, ys, _, n, ps, mut losss) = common::setup_test();

    /*
    let num_cpus = 1;
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_cpus)
        .build_global()
        .unwrap();
         */
    println!("num_thread set: {}", rayon::current_num_threads());

    c.bench_function("no simd", |b| {
        b.iter(|| {
            calculate::calculate_loss_gt_simd(
                black_box(&mut losss),
                black_box(&predictions),
                black_box(&ps),
                black_box(&ys),
                black_box(n),
            )
        })
        //b.iter(|| calculate::calculate_loss_gt_nosimd(&mut losss, &predictions, &ps, &ys, n))
    });
}
*/

// https://bheisler.github.io/criterion.rs/book/user_guide/comparing_functions.html
// output plot could be regressed time or scatter plot. If fast, regressed time is output
// https://bheisler.github.io/criterion.rs/book/user_guide/advanced_configuration.html
fn bench_calculate_loss_gt_comp(c: &mut Criterion) {
    println!("num_thread set: {}", rayon::current_num_threads());

    // m=100k
    let (dataset, ps, mut losss) = common::setup_test_n100k();
    let m = dataset.snvs().snvs_n();

    let mut group = c.benchmark_group("calculate_loss_gt");
    // no need to call black_box in bench_with_input
    // simd
    group.bench_with_input(BenchmarkId::new("Simd", m), &m, |b, _| {
        b.iter(|| {
            calc::calculate_loss_gt_constada_simd(
                &mut losss,
                &dataset.genot(),
                &ps,
                &dataset.samples().phe(),
                BoostParam::new_type1(),
            )
        })
    });
    // no simd
    group.bench_with_input(BenchmarkId::new("No Simd", m), &m, |b, _| {
        b.iter(|| {
            calc::calculate_loss_gt_constada_nosimd(
                &mut losss,
                &dataset.genot(),
                &ps,
                &dataset.samples().phe(),
                BoostParam::new_type1(),
            )
        })
    });

    // m=10k
    let (dataset, ps, mut losss) = common::setup_test_n10k();
    let m = dataset.snvs().snvs_n();
    // simd
    group.bench_with_input(BenchmarkId::new("Simd", m), &m, |b, _| {
        b.iter(|| {
            calc::calculate_loss_gt_constada_simd(
                &mut losss,
                &dataset.genot(),
                &ps,
                &dataset.samples().phe(),
                BoostParam::new_type1(),
            )
        })
    });
    // no simd
    group.bench_with_input(BenchmarkId::new("No Simd", m), &m, |b, _| {
        b.iter(|| {
            calc::calculate_loss_gt_constada_nosimd(
                &mut losss,
                &dataset.genot(),
                &ps,
                &dataset.samples().phe(),
                BoostParam::new_type1(),
            )
        })
    });

    group.finish();

    /*
    c.bench_with_input(BenchmarkId::new("input_example", size), &size, |b, &s| {
        b.iter(|| do_something(s));
    });
    */
}

/*
criterion_group!(
    benches,
    bench_calculate_loss_gt_simd,
    bench_calculate_loss_gt_nosimd
);
*/

criterion_group!(benches, bench_calculate_loss_gt_comp);
criterion_main!(benches);
