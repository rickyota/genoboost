//!
//! Compare if and match.
//!
//! Almost the same but If is slightly faster(~2%).
//!
//

/*
#[macro_use]
extern crate criterion;

//use criterion::black_box;
use criterion::{BenchmarkId, Criterion};

use genetics_rust::genot;
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

pub fn predictions_both_to_count(pred_m: &[B8], ni: usize, len_n: usize) -> u8 {
    let pred_s0m = &pred_m[..len_n];
    let pred_s1m = &pred_m[len_n..];
    let pred0 = bit::bget(pred_s0m, ni);
    let pred1 = bit::bget(pred_s1m, ni);

    if (!pred0) & (!pred1) {
        0
    } else if pred0 & (!pred1) {
        1
    } else if pred0 & pred1 {
        2
    } else {
        3
    }
}

// ~5% slower
pub fn predictions_both_to_count2(pred_m: &[B8], ni: usize, len_n: usize) -> u8 {
    let pred_s0m = &pred_m[..len_n];
    let pred_s1m = &pred_m[len_n..];
    let pred0 = bit::bget(pred_s0m, ni);
    let pred1 = bit::bget(pred_s1m, ni);

    match (pred0, pred1) {
        (false, false) => 0,
        (true, false) => 1,
        (true, true) => 2,
        (false, true) => 3,
    }
}

pub fn bench_pred_to_count(c: &mut Criterion) {
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
    // set predicitons
    predictions.chunks_mut(2 * len_n).for_each(|pred| {
        (0..n).for_each(|ni| set_predictions_both_count(pred, counts[ni], ni, len_n))
        //(0..n).for_each(|ni| genot_twin::set_predictions_both_count(pred, counts[ni], ni, len_n))
    });

    let mut cs = vec![0; n];

    let mut group = c.benchmark_group("pred_to_count");

    group.bench_with_input(BenchmarkId::new("If", m), &m, |b, _| {
        b.iter(|| {
            predictions.chunks_mut(2 * len_n).for_each(|pred| {
                (0..n).for_each(|ni| cs[ni] = predictions_both_to_count(pred, ni, len_n))
                //.for_each(|ni| cs[ni] = genot_twin::predictions_both_to_count(pred, ni, len_n))
            });
        })
    });
    group.bench_with_input(BenchmarkId::new("Match", m), &m, |b, _| {
        b.iter(|| {
            predictions.chunks_mut(2 * len_n).for_each(|pred| {
                (0..n).for_each(|ni| cs[ni] = predictions_both_to_count2(pred, ni, len_n))
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

criterion_group!(benches, bench_pred_to_count);
criterion_main!(benches);

*/

fn main() {}
