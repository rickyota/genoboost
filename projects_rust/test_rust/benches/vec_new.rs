//!
//! Fill value to vec.
//! Compare set_len(), resize() and push()
//!
//! resize and set_len are almost ignorable.
//! <1e-12 ps / 10^11
//! -> Is it true?
//!
//!
//! push() is devestatingly slow.
//! 300 ms / 10^8
//!
//! use resize since safe
//!
//!

#[macro_use]
extern crate criterion;

//use criterion::black_box;
use criterion::{BenchmarkId, Criterion};

mod common;
//use super::set_count_criterion::common;

pub fn bench_vec_new(c: &mut Criterion) {
    let ns = [10usize.pow(8), 10usize.pow(9), 10usize.pow(11)];

    let mut group = c.benchmark_group("vec_new");

    for &n in ns.iter() {
        group.bench_with_input(BenchmarkId::new("set_len", n), &n, |b, _| {
            b.iter(|| {
                let mut vec: Vec<u8> = Vec::with_capacity(n);
                unsafe {
                    vec.set_len(n);
                }
                vec.fill(0x00);
                println!("vec {}", vec[100]);
            })
        });
        group.bench_with_input(BenchmarkId::new("resize", n), &n, |b, _| {
            b.iter(|| {
                let mut vec: Vec<u8> = Vec::with_capacity(n);
                vec.resize(n, 0x00);
                println!("vec {}", vec[100]);
            })
        });

        // push is too slow
        if n <= 10usize.pow(8) {
            group.bench_with_input(BenchmarkId::new("push", n), &n, |b, _| {
                b.iter(|| {
                    let mut vec: Vec<u8> = Vec::with_capacity(n);
                    for _ in 0..n {
                        vec.push(0x00);
                    }
                })
            });
        }
    }

    group.finish();
}

criterion_group!(benches, bench_vec_new);
criterion_main!(benches);
