//!
//!
//! Compare array and bit when converting code.
//!
//! Bit is faster(~9%).
//! This bit is a little slower than set_count() since bit operation is one more.
//!
//!

/*
#[macro_use]
extern crate criterion;

//use criterion::black_box;
use criterion::{BenchmarkId, Criterion};

use genetics_rust::genot;
type B8 = u8;

mod common;
//use super::set_count_criterion::common;

/// for plink bed
/// 00(2); 1,1
/// 10(1); 1,0
/// 11(0); 0,0
/// 01(NA); 0,1
pub fn set_predictions_both_bed(pred_m: &mut [B8], bcode: u8, ni: usize, len_n: usize) {
    bit::bset(&mut pred_m[..len_n], ni, (bcode & 1) == 0);
    bit::bset(&mut pred_m[len_n..], ni, (bcode & 2) == 0);
}

// plink bed code to pred
// 00(2); 1,1
// 01(NA); 0,1
// 10(1); 1,0
// 11(0); 0,0
const BED_PRED0_AR: [bool; 4] = [true, false, true, false];
const BED_PRED1_AR: [bool; 4] = [true, true, false, false];

// ~10% slower
pub fn set_predictions_both_bed2(pred_m: &mut [B8], bcode: u8, ni: usize, len_n: usize) {
    bit::bset(&mut pred_m[..len_n], ni, BED_PRED0_AR[bcode as usize]);
    bit::bset(&mut pred_m[len_n..], ni, BED_PRED1_AR[bcode as usize]);
}

pub fn bench_set_bed(c: &mut Criterion) {
    let m = 100;
    let n = 100_000;

    let mut predictions = common::setup_test_empty(n, m);
    let count = [0u8, 1, 2, 3];
    let counts = count
        .iter()
        .cycle()
        .take(n)
        .map(|x| *x)
        .collect::<Vec<u8>>();

    let len_n = n / 8 + 5;
    //let len_n = genot_twin::len_n(n);

    let mut group = c.benchmark_group("set_bed");

    group.bench_with_input(BenchmarkId::new("Bit", m), &m, |b, _| {
        b.iter(|| {
            predictions.chunks_mut(2 * len_n).for_each(|pred| {
                (0..n).for_each(|ni| {
                    set_predictions_both_bed(pred, counts[ni], ni, len_n)
                    //genot_twin::set_predictions_both_bed(pred, counts[ni], ni, len_n)
                })
            });
        })
    });
    group.bench_with_input(BenchmarkId::new("Array", m), &m, |b, _| {
        b.iter(|| {
            predictions.chunks_mut(2 * len_n).for_each(|pred| {
                (0..n).for_each(|ni| {
                    set_predictions_both_bed2(pred, counts[ni], ni, len_n)
                    //genot_twin::set_predictions_both_bed2(pred, counts[ni], ni, len_n)
                })
            });
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
 */
fn main() {}
