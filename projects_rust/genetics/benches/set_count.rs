//!
//! Compare array and bit when converting code.
//!
//! Bit is faster(~12%).
//!

/*
#[macro_use]
extern crate criterion;

//use criterion::black_box;
use criterion::{BenchmarkId, Criterion};

type B8 = u8;

mod common;
//use super::set_count_criterion::common;

// for benchmark
// both almost the same speed
// count to pred
// 0(b00) => (0,0)
// 1(b01) => (1,0)
// 2(b10) => (1,1)
// 3(b11) => (0,1)
pub fn set_predictions_both_count(pred_m: &mut [B8], count: u8, ni: usize, len_n: usize) {
    bit::bset(&mut pred_m[..len_n], ni, (count >> 1) != (count & 1));
    bit::bset(&mut pred_m[len_n..], ni, (count & 2) == 2);
}

// count to pred
// 0 => (0,0)
// 1 => (1,0)
// 2 => (1,1)
// 3 => (0,1)
const COUNT_PRED0_AR: [bool; 4] = [false, true, true, false];
const COUNT_PRED1_AR: [bool; 4] = [false, false, true, true];

// for benchmark
// ~5% slower
pub fn set_predictions_both_count2(pred_m: &mut [B8], count: u8, ni: usize, len_n: usize) {
    bit::bset(&mut pred_m[..len_n], ni, COUNT_PRED0_AR[count as usize]);
    bit::bset(&mut pred_m[len_n..], ni, COUNT_PRED1_AR[count as usize]);
}

// for benchmark
// comparable
pub fn set_predictions_both_count3(pred_m: &mut [B8], count: u8, ni: usize, len_n: usize) {
    let val0 = if (count == 0) || (count == 3) {
        false
    } else {
        true
    };
    let val1 = if (count == 0) || (count == 1) {
        false
    } else {
        true
    };
    bit::bset(&mut pred_m[..len_n], ni, val0);
    bit::bset(&mut pred_m[len_n..], ni, val1);
}

pub fn bench_set_count(c: &mut Criterion) {
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

    let mut group = c.benchmark_group("set_count");

    group.bench_with_input(BenchmarkId::new("Bit", m), &m, |b, _| {
        b.iter(|| {
            predictions.chunks_mut(2 * len_n).for_each(|pred| {
                (0..n).for_each(|ni| {
                    set_predictions_both_count(pred, counts[ni], ni, len_n)
                    //genot_twin::set_predictions_both_count(pred, counts[ni], ni, len_n)
                })
            });
        })
    });
    group.bench_with_input(BenchmarkId::new("Array", m), &m, |b, _| {
        b.iter(|| {
            predictions.chunks_mut(2 * len_n).for_each(|pred| {
                (0..n).for_each(|ni| {
                    set_predictions_both_count2(pred, counts[ni], ni, len_n)
                    //genot_twin::set_predictions_both_count2(pred, counts[ni], ni, len_n)
                })
            });
        })
    });
    group.bench_with_input(BenchmarkId::new("If", m), &m, |b, _| {
        b.iter(|| {
            predictions.chunks_mut(2 * len_n).for_each(|pred| {
                (0..n).for_each(|ni| {
                    set_predictions_both_count3(pred, counts[ni], ni, len_n)
                    //genot_twin::set_predictions_both_count3(pred, counts[ni], ni, len_n)
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

criterion_group!(benches, bench_set_count);
criterion_main!(benches);
 */

fn main() {}
