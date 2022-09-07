//!
//! Compare loading whole file first, or only load each piece.
//!
//! load part cannot use rayon
//!
//! load whole is fastest even for prop=1%
//!
//! Compare two functions varying variables.
//! [ref](https://bheisler.github.io/criterion.rs/book/user_guide/comparing_functions.html#comparing-functions)
//!

/*
#[macro_use]
extern crate criterion;

//use criterion::black_box;
use criterion::{BenchmarkId, Criterion};

use genetics_rust::genot::BaseGenot;
use genetics_rust::plink;

use rand::Rng;

mod common;
//use super::set_count_criterion::common;

fn setup(fin: &str, prop: f64) -> (usize, usize, Vec<bool>, Vec<bool>) {
    let m_in: usize = plink::compute_num_snv(fin).unwrap();
    println!("m_in: {}", m_in);
    let n_in: usize = plink::compute_num_sample(fin).unwrap();
    println!("n_in: {}", n_in);

    let use_samples = vec![true; n_in];

    let mut rng = rand::thread_rng();

    let mut use_snvs = Vec::with_capacity(n_in);
    for _ in 0..m_in {
        let x: f64 = rng.gen();
        let b = if x < prop { true } else { false };
        use_snvs.push(b);
    }

    let m = use_snvs.iter().filter(|&&x| x == true).count();
    println!("m: {}", m_in);

    (m, n_in, use_snvs, use_samples)
}

pub fn bench_set_bed(c: &mut Criterion) {
    let props = [0.01, 0.1, 1.0];

    //let fin = String::from("../../test/data/1kg_n100000/1kg_n100000");
    //let fin = String::from("../../data/1kg_phase1_all/1kg_phase1_all");
    let fin = String::from("../../data/dataset.4/cad/ukb_imp_chr1.obs");

    //let genot = genot_twin::generate_predictions(&fin, m, n, &use_snvs, &use_samples, use_missing);

    let mut group = c.benchmark_group("load_bed");
    group.sample_size(10);

    for prop in props {
        let (m, n, use_snvs, use_samples) = setup(&fin, prop);
        //let genot = genot_twin::generate_predictions(&fin, m, n, &use_snvs, &use_samples, use_missing);

        group.bench_with_input(BenchmarkId::new("Whole", prop), &prop, |b, _| {
            b.iter(|| {
                let genot = plink::load::generate_predictions_config(
                    &fin,
                    m,
                    n,
                    &use_snvs,
                    &use_samples,
                    false,
                    false,
                );
                genot.len_n();
            })
        });
        group.bench_with_input(BenchmarkId::new("Part", prop), &prop, |b, _| {
            b.iter(|| {
                let genot = plink::load::generate_predictions_config(
                    &fin,
                    m,
                    n,
                    &use_snvs,
                    &use_samples,
                    false,
                    true,
                );
                genot.len_n();
            })
        });
    }

    group.finish();

    /*
    c.bench_with_input(BenchmarkId::new("input_example", size), &size, |b, &s| {
        b.iter(|| do_something(s));
    });
    */
}

criterion_group!(benches, bench_set_bed);
criterion_main!(benches);

 */

// avoid error by comment out
fn main() {}
