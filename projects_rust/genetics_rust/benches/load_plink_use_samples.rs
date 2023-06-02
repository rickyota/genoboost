//!
//! Compare time to load snv using `use_samples` or `None`
//!
//! Almost the same. ~1% faster
//!
//!
//!
//!

#[macro_use]
extern crate criterion;

//use criterion::black_box;
use criterion::{BenchmarkId, Criterion};

use genetics::genot::prelude::*;
use genetics::genot::BaseGenot;
use genetics::plink;

mod common;
//use super::set_count_criterion::common;

pub fn bench_set_bed(c: &mut Criterion) {
    let n = 500_000;
    let mut g = GenotSnv::new_empty(n);

    let buf_mi = vec![255u8; n];
    let use_samples = vec![true; n];

    let mut group = c.benchmark_group("use_samples");

    group.bench_with_input(BenchmarkId::new("Whole", n), &n, |b, _| {
        b.iter(|| {
            plink::load::assign_pred_from_bed(&mut g.as_genot_snv_mut_snv(), &buf_mi, None);
            g.len_n();
        })
    });
    group.bench_with_input(BenchmarkId::new("Part", n), &n, |b, _| {
        b.iter(|| {
            plink::load::assign_pred_from_bed(
                &mut g.as_genot_snv_mut_snv(),
                &buf_mi,
                Some(&use_samples),
            );
            g.len_n();
        })
    });

    group.finish();

    /*
    c.bench_with_input(BenchmarkId::new("input_example", size), &size, |b, &s| {
        b.iter(|| do_something(s));
    });
    */
}

criterion_group!(benches, bench_set_bed);
criterion_main!(benches);

// avoid error by comment out
//fn main() {}
