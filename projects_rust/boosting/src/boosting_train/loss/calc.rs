use super::epsilon;
use super::table;
use super::LossStruct;
use crate::boosting_train::coefficient;
use crate::boosting_train::sample_weight::SampleWeight;
use crate::{BoostParam, ContingencyTable, Eps};
use genetics::genot::prelude::*;
use genetics::samples::prelude::*;
use genetics::wgt::Coef;
use rayon::prelude::*;

use std::collections::HashSet;

pub fn calculate_loss_table4(table: ContingencyTable) -> f64 {
    let (a, b, c, d) = table.four();
    //let (a, b, c, d) = table.four();
    2.0 * ((a * b).sqrt() + (c * d).sqrt())
}

pub fn calculate_loss_table5(table: ContingencyTable) -> f64 {
    let (d1, n1, d0, n0, m) = table.five();
    m + 2.0 * ((d1 * n1).sqrt() + (d0 * n0).sqrt())
}

pub fn calculate_loss_table7(table: ContingencyTable) -> f64 {
    let (d2, n2, d1, n1, d0, n0, m) = table.seven();
    m + 2.0 * ((d2 * n2).sqrt() + (d1 * n1).sqrt() + (d0 * n0).sqrt())
}

pub fn calculate_loss_table7_or_5(table: ContingencyTable) -> f64 {
    if let ContingencyTable::Seven(_) = table {
        calculate_loss_table7(table)
    } else if let ContingencyTable::Five(_) = table {
        calculate_loss_table5(table)
    } else {
        panic!("");
    }
}

// TODO
pub fn calculate_loss_ab(_table: ContingencyTable) -> f64 {
    f64::NAN
}

// cannot run on M1
// ref: https://doc.rust-lang.org/stable/core/arch/#examples
// aarch64: mac M1
//#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
// use --target x86_64-apple-darwin instead
/// ps should be aligned. alignment of predict is not necessary.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn calculate_loss_gt_constada_simd_sm(
    gsnv: &GenotSnvRef,
    ps_pad: &[f64],
    phe: &Phe,
    ef_: (f64, f64),
    epsilons: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    eps: Option<Eps>,
) -> (f64, f64) {
    // aarch64 is completely different: https://moshg.github.io/rust-std-ja/core/arch/aarch64/index.html
    //#[cfg(target_arch = "aarch64")]
    //use std::arch::aarch64::*;
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    //use std::convert::TryInto;

    let n = phe.n();

    let ys = phe.inner();
    let pred_s0m = gsnv.predict_s(0);
    let pred_s1m = gsnv.predict_s(1);

    // b: 0x80808080 = -x: value as i32
    // ~b = x - 1
    // x = ~b + 1
    //log::debug!("0x7f7f7f7f + 1: {:#010x}", 0x7f7f7f7f + 1 as i32);
    //log::debug!("0x7f7f7f7f + 1: {}", 0x7f7f7f7f + 1 as i32);
    //log::debug!("0x7f7f7f80: {}", 0x7f7f7f80 as i32);

    let bit_ext_mask: __m256i = _mm256_set_epi32(
        0x10101010,
        0x01010101,
        0x20202020,
        0x02020202,
        0x40404040,
        0x04040404,
        -0x7f7f7f80, // =0x80808080
        0x08080808,
    );
    let zerod: __m256d = _mm256_setzero_pd();

    let mut a_sum_s0_acc = _mm256_setzero_pd();
    let mut b_sum_s0_acc = _mm256_setzero_pd();
    let mut a_sum_s1_acc = _mm256_setzero_pd();
    let mut b_sum_s1_acc = _mm256_setzero_pd();

    // bi=0-3, shift=24,16,8,0
    //let shift_v = _mm256_set1_epi32(shift);
    let shifts: [__m256i; 4] = [
        _mm256_set1_epi32(24),
        _mm256_set1_epi32(16),
        _mm256_set1_epi32(8),
        _mm256_set1_epi32(0),
    ];

    // when n=33, (n/32+1)=2
    // sample index to look through is [0..64)
    // size of y: n/8+5=9: [0..64), which is included
    for ni in 0..(n / 32 + 1) {
        //log::debug!("ni {}", ni);

        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        // 1. use _mm256_set_epi8
        // o2. use from_be() and use set1 <- assume to be fast??
        let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        //log::debug!("mis 0: {:#010b}", mis_s0m[4 * ni]);
        //log::debug!("mis 1: {:#010b}", mis_s0m[4 * ni + 1]);
        //log::debug!("mis 2: {:#010b}", mis_s0m[4 * ni + 2]);
        //log::debug!("mis 3: {:#010b}", mis_s0m[4 * ni + 3]);
        //log::debug!("mis b32: {:#034b}", mis_s0_b32);
        //log::debug!("misv {:?}", misv_s0_32);

        let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let yv_32 = _mm256_set1_epi32(ys_b32 as i32);

        // Asum = y & pred
        // Bsum = !y & pred
        let flagv_a_s0_32 = _mm256_and_si256(yv_32, predv_s0_32);
        let flagv_b_s0_32 = _mm256_andnot_si256(yv_32, predv_s0_32);
        let flagv_a_s1_32 = _mm256_and_si256(yv_32, predv_s1_32);
        let flagv_b_s1_32 = _mm256_andnot_si256(yv_32, predv_s1_32);

        // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
        // abcdefgh(x4)|...
        // -> extracted (at highest position)
        // 000d0000(x4)|...
        // -> mask
        // 11111111|00000000|00000000|11111111|...
        let flagv_a_s0_32_ext = _mm256_and_si256(flagv_a_s0_32, bit_ext_mask);
        let flagv_b_s0_32_ext = _mm256_and_si256(flagv_b_s0_32, bit_ext_mask);
        let flagv_a_s1_32_ext = _mm256_and_si256(flagv_a_s1_32, bit_ext_mask);
        let flagv_b_s1_32_ext = _mm256_and_si256(flagv_b_s1_32, bit_ext_mask);

        let take_mask_a_s0_32 = _mm256_cmpeq_epi8(flagv_a_s0_32_ext, bit_ext_mask);
        let take_mask_b_s0_32 = _mm256_cmpeq_epi8(flagv_b_s0_32_ext, bit_ext_mask);
        let take_mask_a_s1_32 = _mm256_cmpeq_epi8(flagv_a_s1_32_ext, bit_ext_mask);
        let take_mask_b_s1_32 = _mm256_cmpeq_epi8(flagv_b_s1_32_ext, bit_ext_mask);

        // bi=0-3, shift=24,16,8,0
        //const SHIFTS: &'static [i32] = &[24, 16, 8, 0];

        for bi in 0usize..4 {
            // DCBADCBA...DCBA
            // -> b=1: for B
            // BA00BA00...BA00

            let shift_v = shifts[bi];
            //let shift: i32 = 8 * (3 - bi as i32);
            //let shift_v = _mm256_set1_epi32(shift);
            let take_mask_a_s0 = _mm256_sllv_epi32(take_mask_a_s0_32, shift_v);
            let take_mask_b_s0 = _mm256_sllv_epi32(take_mask_b_s0_32, shift_v);
            let take_mask_a_s1 = _mm256_sllv_epi32(take_mask_a_s1_32, shift_v);
            let take_mask_b_s1 = _mm256_sllv_epi32(take_mask_b_s1_32, shift_v);

            // TODO: now raised issue on
            // somehow, only nightly requires `shift` to be constant
            // look like slli requires `shift` to be constant...
            /*
            let shift: i32 = 8 * (3 - bi as i32);
            let take_mask_a_s0 = _mm256_slli_epi32(take_mask_a_s0_32, shift);
            let take_mask_b_s0 = _mm256_slli_epi32(take_mask_b_s0_32, shift);
            let take_mask_a_s1 = _mm256_slli_epi32(take_mask_a_s1_32, shift);
            let take_mask_b_s1 = _mm256_slli_epi32(take_mask_b_s1_32, shift);
            */

            //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
            //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
            //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
            //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);

            // ps[i].as_ptr() should be enough
            let psv_lo_ptr = ps_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let psv_hi_ptr = ps_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();

            let psv_lo: __m256d = _mm256_load_pd(psv_lo_ptr as *const _);
            let psv_hi: __m256d = _mm256_load_pd(psv_hi_ptr as *const _);

            // first for low
            let ps_masked_a_s0_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_a_s0));
            let ps_masked_b_s0_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_b_s0));
            let ps_masked_a_s1_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_a_s1));
            let ps_masked_b_s1_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_b_s1));

            //log::debug!("ps a s0 lo {:?}", ps_masked_a_s0_lo);
            //log::debug!("ps a s1 lo {:?}", ps_masked_a_s1_lo);
            //log::debug!("ps b s0 lo {:?}", ps_masked_b_s0_lo);
            //log::debug!("ps b s1 lo {:?}", ps_masked_b_s1_lo);

            // for high
            let take_mask_a_s0_hi = _mm256_slli_epi64(take_mask_a_s0, 32);
            let take_mask_b_s0_hi = _mm256_slli_epi64(take_mask_b_s0, 32);
            let take_mask_a_s1_hi = _mm256_slli_epi64(take_mask_a_s1, 32);
            let take_mask_b_s1_hi = _mm256_slli_epi64(take_mask_b_s1, 32);

            let ps_masked_a_s0_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_a_s0_hi));
            let ps_masked_b_s0_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_b_s0_hi));
            let ps_masked_a_s1_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_a_s1_hi));
            let ps_masked_b_s1_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_b_s1_hi));

            //log::debug!("a s0 hi {:?}", ps_masked_a_s0_hi);

            a_sum_s0_acc = _mm256_add_pd(a_sum_s0_acc, ps_masked_a_s0_lo);
            a_sum_s0_acc = _mm256_add_pd(a_sum_s0_acc, ps_masked_a_s0_hi);
            b_sum_s0_acc = _mm256_add_pd(b_sum_s0_acc, ps_masked_b_s0_lo);
            b_sum_s0_acc = _mm256_add_pd(b_sum_s0_acc, ps_masked_b_s0_hi);
            a_sum_s1_acc = _mm256_add_pd(a_sum_s1_acc, ps_masked_a_s1_lo);
            a_sum_s1_acc = _mm256_add_pd(a_sum_s1_acc, ps_masked_a_s1_hi);
            b_sum_s1_acc = _mm256_add_pd(b_sum_s1_acc, ps_masked_b_s1_lo);
            b_sum_s1_acc = _mm256_add_pd(b_sum_s1_acc, ps_masked_b_s1_hi);
        }
    }

    // sum 4 double horizontally to get the whole sum
    a_sum_s0_acc = _mm256_hadd_pd(a_sum_s0_acc, a_sum_s0_acc);
    b_sum_s0_acc = _mm256_hadd_pd(b_sum_s0_acc, b_sum_s0_acc);
    a_sum_s1_acc = _mm256_hadd_pd(a_sum_s1_acc, a_sum_s1_acc);
    b_sum_s1_acc = _mm256_hadd_pd(b_sum_s1_acc, b_sum_s1_acc);

    // 1. any way to hadd??
    // 2. _mm256_extractf128_pd and _mm256_cvtsd_f64: get 64:0

    let a_sum_s0: f64 =
        _mm256_cvtsd_f64(a_sum_s0_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(a_sum_s0_acc, 1));
    let b_sum_s0: f64 =
        _mm256_cvtsd_f64(b_sum_s0_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(b_sum_s0_acc, 1));
    let a_sum_s1: f64 =
        _mm256_cvtsd_f64(a_sum_s1_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(a_sum_s1_acc, 1));
    let b_sum_s1: f64 =
        _mm256_cvtsd_f64(b_sum_s1_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(b_sum_s1_acc, 1));

    let (abcd_s0, _) = table::adjust_eps_table4((a_sum_s0, b_sum_s0), ef_, epsilons, eps);
    let (abcd_s1, _) = table::adjust_eps_table4((a_sum_s1, b_sum_s1), ef_, epsilons, eps);

    let loss_s0 = calculate_loss_table4(abcd_s0);
    let loss_s1 = calculate_loss_table4(abcd_s1);

    (loss_s0, loss_s1)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
/// losss[..,2m,2mi+1..] 2mi: dom, 2mi+1: rec
pub unsafe fn calculate_loss_gt_constada_simd(
    losss: &mut LossStruct,
    genot: &Genot,
    sample_weight: &SampleWeight,
    //ps_pad: &[f64],
    phe: &Phe,
    boost_param: BoostParam,
    skip_snv: &HashSet<usize>,
) {
    assert_eq!(losss.inner_mut().len(), genot.m() * 2);
    assert_eq!(phe.n(), genot.n());

    //let n = phe.n();
    let epsilons = epsilon::calculate_epsilons(sample_weight.ps().unwrap(), phe, boost_param.eps());
    //let epsilons = epsilon::calculate_epsilons(&ps_pad[..n], phe, boost_param.eps());
    log::debug!("epsilon case, cont: {:.4e},{:.4e}", epsilons.0, epsilons.1);

    // first calculate E:=case_ws_sum and F:=control_ws_sum
    let ef_ = table::calculate_ef_sum(sample_weight.ps().unwrap(), phe);

    losss
        .inner_mut()
        .par_chunks_mut(2)
        .enumerate()
        .for_each(|(mi, loss)| {
            if skip_snv.contains(&mi) {
                loss[0] = f64::MAX;
                loss[1] = f64::MAX;
            } else {
                let loss_ = calculate_loss_gt_constada_simd_sm(
                    &genot.to_genot_snv(mi),
                    sample_weight.ps_pad().unwrap(),
                    //ps_pad,
                    phe,
                    ef_,
                    epsilons,
                    boost_param.eps(),
                );
                loss[0] = loss_.0;
                loss[1] = loss_.1;
            }
        });
}

fn calculate_loss_gt_constada_nosimd_sm(
    //pred_s: &GenotSnvRef,
    pred_sm: &[u8],
    //pred_sm: &[B8],
    ps: &[f64],
    phe: &Phe,
    ef_: (f64, f64),
    epsilons: (f64, f64),
    eps: Option<Eps>,
) -> f64 {
    let (table4, _) = table::calculate_table4_epsilons(pred_sm, ps, phe, ef_, epsilons, eps);

    calculate_loss_table4(table4)
}

pub fn calculate_loss_gt_constada_nosimd(
    losss: &mut LossStruct,
    genot: &Genot,
    //predictions: &GenotBi<Vec<u8>>,
    //predictions: &[B8],
    sample_weight: &SampleWeight,
    //ps: &[f64],
    phe: &Phe,
    boost_param: BoostParam,
) {
    assert_eq!(losss.inner_mut().len(), genot.m() * 2);

    let epsilons = epsilon::calculate_epsilons(sample_weight.ps().unwrap(), phe, boost_param.eps());
    //log::debug!("epsilon case, cont: {:.2e},{:.2e}", epsilons.0, epsilons.1);

    // first calculate E:=case_ws_sum and F:=control_ws_sum
    let ef_ = table::calculate_ef_sum(sample_weight.ps().unwrap(), phe);

    losss
        .inner_mut()
        .par_iter_mut()
        .enumerate()
        .for_each(|(li, loss)| {
            *loss = calculate_loss_gt_constada_nosimd_sm(
                &genot.vals_snv_s_u8(li / 2, li % 2),
                //predictions.predict_snv_s(li / 2, li % 2),
                sample_weight.ps().unwrap(),
                //ps,
                phe,
                ef_,
                epsilons,
                boost_param.eps(),
            );
        })
    //log::debug!("losss {:?}", losss);
    //losss
}

/*
/// calculate loss for FreeModelMissing
///  loss = M + 2*(sqrt(D2 * N2)+sqrt(D1 * N1)+sqrt(D0 *N0))
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn calculate_table7_sum_simd(
    predicts: &GenotSnvRef,
    ps: &[f64],
    phe: &Phe,
) -> ContingencyTable {
    // aarch64 is completely different: https://moshg.github.io/rust-std-ja/core/arch/aarch64/index.html
    //#[cfg(target_arch = "aarch64")]
    //use std::arch::aarch64::*;
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    use std::convert::TryInto;

    let n = phe.n();
    let ys = phe.inner();
    let pred_s0m = predicts.predict_s(0);
    let pred_s1m = predicts.predict_s(1);

    // b: 0x80808080 = -x: value as i32
    // ~b = x - 1
    // x = ~b + 1
    //log::debug!("0x7f7f7f7f + 1: {:#010x}", 0x7f7f7f7f + 1 as i32);
    //log::debug!("0x7f7f7f7f + 1: {}", 0x7f7f7f7f + 1 as i32);
    //log::debug!("0x7f7f7f80: {}", 0x7f7f7f80 as i32);

    let bit_ext_mask: __m256i = _mm256_set_epi32(
        0x10101010,
        0x01010101,
        0x20202020,
        0x02020202,
        0x40404040,
        0x04040404,
        -0x7f7f7f80, // =0x80808080
        0x08080808,
    );
    let zerod: __m256d = _mm256_setzero_pd();
    let zeros: __m256i = _mm256_setzero_si256();
    let ones: __m256i = _mm256_cmpeq_epi32(zeros, zeros);

    let mut d2_sum_acc = _mm256_setzero_pd();
    let mut n2_sum_acc = _mm256_setzero_pd();
    let mut d1_sum_acc = _mm256_setzero_pd();
    let mut n1_sum_acc = _mm256_setzero_pd();
    let mut d0_sum_acc = _mm256_setzero_pd();
    let mut n0_sum_acc = _mm256_setzero_pd();

    // bi=0-3, shift=24,16,8,0
    //let shift_v = _mm256_set1_epi32(shift);
    let shifts: [__m256i; 4] = [
        _mm256_set1_epi32(24),
        _mm256_set1_epi32(16),
        _mm256_set1_epi32(8),
        _mm256_set1_epi32(0),
    ];

    // when n=33, (n/32+1)=2
    // sample index to look through is [0..64)
    // size of y: n/8+5=9: [0..64), which is included
    for ni in 0..(n / 32 + 1) {
        //log::debug!("ni {}", ni);

        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        // 1. use _mm256_set_epi8
        // o2. use from_be() and use set1 <- assume to be fast??
        let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let mis_s0_b32 = u32::from_be_bytes(mis_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let mis_s1_b32 = u32::from_be_bytes(mis_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        //log::debug!("mis 0: {:#010b}", mis_s0m[4 * ni]);
        //log::debug!("mis 1: {:#010b}", mis_s0m[4 * ni + 1]);
        //log::debug!("mis 2: {:#010b}", mis_s0m[4 * ni + 2]);
        //log::debug!("mis 3: {:#010b}", mis_s0m[4 * ni + 3]);
        //log::debug!("mis b32: {:#034b}", mis_s0_b32);
        //log::debug!("misv {:?}", misv_s0_32);

        let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let yv_32 = _mm256_set1_epi32(ys_b32 as i32);

        // bitwise not is fastest in xor by xor(v, ones)
        // D2sum = y & pred0 & pred1 = y & ( pred0 & pred1 )
        // N2sum = !y & pred0 & pred1 = !y & ( pred0 & pred1 )
        // D1sum = y & pred0 & !pred1 = y & ( !pred1 & pred0 )
        // N1sum = !y & pred0 & !pred1 = !y & ( !pred1 & pred0 )
        // D0sum = y & !pred0 & !pred1 = !pred0 & ( !pred1 & y )
        // N0sum = !y & !pred0 & !pred1 = !y & ( !pred0 & (!pred1) )
        let flagv_d2_32 = _mm256_and_si256(yv_32, _mm256_and_si256(predv_s0_32, predv_s1_32));
        let flagv_n2_32 = _mm256_andnot_si256(yv_32, _mm256_and_si256(predv_s0_32, predv_s1_32));
        let flagv_d1_32 = _mm256_and_si256(yv_32, _mm256_andnot_si256(predv_s1_32, predv_s0_32));
        let flagv_n1_32 = _mm256_andnot_si256(yv_32, _mm256_andnot_si256(predv_s1_32, predv_s0_32));
        let flagv_d0_32 = _mm256_andnot_si256(predv_s0_32, _mm256_andnot_si256(predv_s1_32, yv_32));
        let flagv_n0_32 = _mm256_andnot_si256(
            yv_32,
            _mm256_andnot_si256(predv_s0_32, _mm256_xor_si256(predv_s1_32, ones)),
        );

        // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
        // abcdefgh(x4)|...
        // -> extracted (at highest position)
        // 000d0000(x4)|...
        // -> mask
        // 11111111|00000000|00000000|11111111|...
        let flagv_d2_32_ext = _mm256_and_si256(flagv_d2_32, bit_ext_mask);
        let flagv_n2_32_ext = _mm256_and_si256(flagv_n2_32, bit_ext_mask);
        let flagv_d1_32_ext = _mm256_and_si256(flagv_d1_32, bit_ext_mask);
        let flagv_n1_32_ext = _mm256_and_si256(flagv_n1_32, bit_ext_mask);
        let flagv_d0_32_ext = _mm256_and_si256(flagv_d0_32, bit_ext_mask);
        let flagv_n0_32_ext = _mm256_and_si256(flagv_n0_32, bit_ext_mask);

        let take_mask_d2_32 = _mm256_cmpeq_epi8(flagv_d2_32_ext, bit_ext_mask);
        let take_mask_n2_32 = _mm256_cmpeq_epi8(flagv_n2_32_ext, bit_ext_mask);
        let take_mask_d1_32 = _mm256_cmpeq_epi8(flagv_d1_32_ext, bit_ext_mask);
        let take_mask_n1_32 = _mm256_cmpeq_epi8(flagv_n1_32_ext, bit_ext_mask);
        let take_mask_d0_32 = _mm256_cmpeq_epi8(flagv_d0_32_ext, bit_ext_mask);
        let take_mask_n0_32 = _mm256_cmpeq_epi8(flagv_n0_32_ext, bit_ext_mask);

        // bi=0-3, shift=24,16,8,0
        //const SHIFTS: &'static [i32] = &[24, 16, 8, 0];

        for bi in 0usize..4 {
            // DCBADCBA...DCBA
            // -> b=1: for B
            // BA00BA00...BA00

            let shift_v = shifts[bi];
            let take_mask_d2 = _mm256_sllv_epi32(take_mask_d2_32, shift_v);
            let take_mask_n2 = _mm256_sllv_epi32(take_mask_n2_32, shift_v);
            let take_mask_d1 = _mm256_sllv_epi32(take_mask_d1_32, shift_v);
            let take_mask_n1 = _mm256_sllv_epi32(take_mask_n1_32, shift_v);
            let take_mask_d0 = _mm256_sllv_epi32(take_mask_d0_32, shift_v);
            let take_mask_n0 = _mm256_sllv_epi32(take_mask_n0_32, shift_v);

            //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
            //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
            //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
            //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);

            let psv_lo_ptr = ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let psv_hi_ptr = ps[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            let psv_lo: __m256d = _mm256_load_pd(psv_lo_ptr as *const _);
            let psv_hi: __m256d = _mm256_load_pd(psv_hi_ptr as *const _);

            //log::debug!("ps lo {:?}", psv_lo);
            //log::debug!("ps hi {:?}", psv_hi);

            // first for low
            let ps_masked_d2_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_d2));
            let ps_masked_n2_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_n2));
            let ps_masked_d1_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_d1));
            let ps_masked_n1_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_n1));
            let ps_masked_d0_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_d0));
            let ps_masked_n0_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_n0));

            d2_sum_acc = _mm256_add_pd(d2_sum_acc, ps_masked_d2_lo);
            n2_sum_acc = _mm256_add_pd(n2_sum_acc, ps_masked_n2_lo);
            d1_sum_acc = _mm256_add_pd(d1_sum_acc, ps_masked_d1_lo);
            n1_sum_acc = _mm256_add_pd(n1_sum_acc, ps_masked_n1_lo);
            d0_sum_acc = _mm256_add_pd(d0_sum_acc, ps_masked_d0_lo);
            n0_sum_acc = _mm256_add_pd(n0_sum_acc, ps_masked_n0_lo);

            //log::debug!("ps a s0 lo {:?}", ps_masked_a_s0_lo);
            //log::debug!("ps a s1 lo {:?}", ps_masked_a_s1_lo);
            //log::debug!("ps b s0 lo {:?}", ps_masked_b_s0_lo);
            //log::debug!("ps b s1 lo {:?}", ps_masked_b_s1_lo);

            // for high
            let take_mask_d2_hi = _mm256_slli_epi64(take_mask_d2, 32);
            let take_mask_n2_hi = _mm256_slli_epi64(take_mask_n2, 32);
            let take_mask_d1_hi = _mm256_slli_epi64(take_mask_d1, 32);
            let take_mask_n1_hi = _mm256_slli_epi64(take_mask_n1, 32);
            let take_mask_d0_hi = _mm256_slli_epi64(take_mask_d0, 32);
            let take_mask_n0_hi = _mm256_slli_epi64(take_mask_n0, 32);

            let ps_masked_d2_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_d2_hi));
            let ps_masked_n2_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_n2_hi));
            let ps_masked_d1_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_d1_hi));
            let ps_masked_n1_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_n1_hi));
            let ps_masked_d0_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_d0_hi));
            let ps_masked_n0_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_n0_hi));

            //log::debug!("a s0 hi {:?}", ps_masked_a_s0_hi);

            d2_sum_acc = _mm256_add_pd(d2_sum_acc, ps_masked_d2_hi);
            n2_sum_acc = _mm256_add_pd(n2_sum_acc, ps_masked_n2_hi);
            d1_sum_acc = _mm256_add_pd(d1_sum_acc, ps_masked_d1_hi);
            n1_sum_acc = _mm256_add_pd(n1_sum_acc, ps_masked_n1_hi);
            d0_sum_acc = _mm256_add_pd(d0_sum_acc, ps_masked_d0_hi);
            n0_sum_acc = _mm256_add_pd(n0_sum_acc, ps_masked_n0_hi);
        }
    }

    // sum 4 double horizontally to get the whole sum
    d2_sum_acc = _mm256_hadd_pd(d2_sum_acc, d2_sum_acc);
    n2_sum_acc = _mm256_hadd_pd(n2_sum_acc, n2_sum_acc);
    d1_sum_acc = _mm256_hadd_pd(d1_sum_acc, d1_sum_acc);
    n1_sum_acc = _mm256_hadd_pd(n1_sum_acc, n1_sum_acc);
    d0_sum_acc = _mm256_hadd_pd(d0_sum_acc, d0_sum_acc);
    n0_sum_acc = _mm256_hadd_pd(n0_sum_acc, n0_sum_acc);

    // 1. any way to hadd??
    // 2. _mm256_extractf128_pd and _mm256_cvtsd_f64: get 64:0

    let d2_sum: f64 =
        _mm256_cvtsd_f64(d2_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(d2_sum_acc, 1));
    let n2_sum: f64 =
        _mm256_cvtsd_f64(n2_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(n2_sum_acc, 1));
    let d1_sum: f64 =
        _mm256_cvtsd_f64(d1_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(d1_sum_acc, 1));
    let n1_sum: f64 =
        _mm256_cvtsd_f64(n1_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(n1_sum_acc, 1));
    let d0_sum: f64 =
        _mm256_cvtsd_f64(d0_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(d0_sum_acc, 1));
    let n0_sum: f64 =
        _mm256_cvtsd_f64(n0_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(n0_sum_acc, 1));

    let m_sum: f64 = 1.0 - (d2_sum + n2_sum + d1_sum + n1_sum + d0_sum + n0_sum);
    // TODO: check m_sum>-MACHNE_EPSILON
    let m_sum = m_sum.max(0.0);

    ContingencyTable::new_seven((d2_sum, n2_sum, d1_sum, n1_sum, d0_sum, n0_sum, m_sum))
}
*/

/// calculate loss for FreeModelMissing
///  loss = M + 2*(sqrt(D2 * N2)+sqrt(D1 * N1)+sqrt(D0 *N0))
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn calculate_loss_gt_freemodelmissing_simd_sm(
    gsnv: &GenotSnvRef,
    ps_pad: &[f64],
    phe: &Phe,
    epsilons: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    eps: Option<Eps>,
) -> f64 {
    let table7_ori = table::calculate_table7_sum_simd(gsnv, ps_pad, phe);

    if eps.is_none() {
        let loss_ = calculate_loss_table7(table7_ori);
        return loss_;
    }

    // TODO: cleaner
    if eps.is_some() && eps.unwrap().dom() {
        let (table7_or_5, _) = table::adjust_eps_table7_dom(
            table7_ori, //(d2_sum, n2_sum, d1_sum, n1_sum, d0_sum, n0_sum),
            eps.unwrap(),
        );

        let loss_ = calculate_loss_table7_or_5(table7_or_5);

        loss_
    } else {
        let (table7, _) = table::adjust_eps_table7_nondom(
            table7_ori, //(d2_sum, n2_sum, d1_sum, n1_sum, d0_sum, n0_sum),
            epsilons,
            eps.unwrap(),
        );

        //if mi % 20000 == 0 {
        //    log::debug!("table7 {:?}", table7);
        //    log::debug!("cont table {:?}", predicts.stat_contingency_table(phe))
        //}
        //debug
        //log::debug!("table7 {:?}", table7);

        let loss_ = calculate_loss_table7(table7);

        loss_
    }
}

// unsafe is necessary with #[target_feature]
//pub unsafe fn calculate_loss_gt_freemodelmissing_simd(
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn calculate_loss_gt_freemodelmissing_simd(
    losss: &mut LossStruct,
    genot: &Genot,
    sample_weight: &SampleWeight,
    //ps_pad: &[f64],
    phe: &Phe,
    boost_param: BoostParam,
    skip_snv: &HashSet<usize>,
) {
    //let n = phe.n();
    let epsilons = epsilon::calculate_epsilons(sample_weight.ps().unwrap(), phe, boost_param.eps());
    log::debug!("epsilon case, cont: {:.4e},{:.4e}", epsilons.0, epsilons.1);

    assert_eq!(losss.inner_mut().len(), genot.m());

    //unsafe {
    losss
        .inner_mut()
        .par_iter_mut()
        .enumerate()
        .for_each(|(mi, loss)| {
            if skip_snv.contains(&mi) {
                *loss = f64::MAX;
            } else {
                *loss = calculate_loss_gt_freemodelmissing_simd_sm(
                    &genot.to_genot_snv(mi),
                    sample_weight.ps_pad().unwrap(),
                    //ps_pad,
                    phe,
                    epsilons,
                    boost_param.eps(),
                )
            }
        });
    //}
}

/// compute decreasing loss of logistic loss function not least-square loss
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn calculate_loss_gt_logit_simd_sm(
    gsnv: &GenotSnvRef,
    coef: &Coef,
    wls_pad: &[f64],
    zs_pad: &[f64],
    use_adjloss: bool,
    loss_max_theory: f64,
) -> f64 {
    let (s0, s1, s2) = match coef {
        Coef::Score3(_) => {
            let (s0, s1, s2) = coef.score3_f64();
            (s0, s1, s2)
        }
        Coef::Score4(_) => {
            let (s0, s1, s2, sm) = coef.score4_f64();
            if sm.abs() > 1e-10 {
                panic!("sm should be 0.0");
            }
            (s0, s1, s2)
        }
        _ => panic!("wrong"),
    };

    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    //use std::convert::TryInto;

    //let n = phe.n();
    //let n_pad=wls.len();
    let n = wls_pad.len() - 32;
    //let ys = phe.inner();
    let pred_s0m = gsnv.predict_s(0);
    let pred_s1m = gsnv.predict_s(1);

    let bit_ext_mask: __m256i = _mm256_set_epi32(
        0x10101010,
        0x01010101,
        0x20202020,
        0x02020202,
        0x40404040,
        0x04040404,
        -0x7f7f7f80, // =0x80808080
        0x08080808,
    );
    let zerod: __m256d = _mm256_setzero_pd();
    let zeros: __m256i = _mm256_setzero_si256();
    let ones: __m256i = _mm256_cmpeq_epi32(zeros, zeros);

    let s2d: __m256d = _mm256_set1_pd(s2);
    let s1d: __m256d = _mm256_set1_pd(s1);
    let s0d: __m256d = _mm256_set1_pd(s0);

    let mut l_acc: __m256d = _mm256_setzero_pd();

    //let mut wzs_sum2_acc = _mm256_setzero_pd();
    //let mut wzs_sum1_acc = _mm256_setzero_pd();
    //let mut wzs_sum0_acc = _mm256_setzero_pd();
    //let mut wls_sum2_acc = _mm256_setzero_pd();
    //let mut wls_sum1_acc = _mm256_setzero_pd();
    //let mut wls_sum0_acc = _mm256_setzero_pd();

    let shifts: [__m256i; 4] = [
        _mm256_set1_epi32(24),
        _mm256_set1_epi32(16),
        _mm256_set1_epi32(8),
        _mm256_set1_epi32(0),
    ];

    for ni in 0..(n / 32 + 1) {
        //log::debug!("ni {}", ni);

        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let yv_32 = _mm256_set1_epi32(ys_b32 as i32);

        // bitwise not is fastest in xor by xor(v, ones)
        // sum2=( pred0 & pred1 )
        // sum1= ( !pred1 & pred0 )
        // sum0= ( !pred0 & (!pred1) )
        let flagv_2_32 = _mm256_and_si256(predv_s0_32, predv_s1_32);
        let flagv_1_32 = _mm256_andnot_si256(predv_s1_32, predv_s0_32);
        let flagv_0_32 = _mm256_andnot_si256(predv_s0_32, _mm256_xor_si256(predv_s1_32, ones));

        // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
        // abcdefgh(x4)|...
        // -> extracted (at highest position)
        // 000d0000(x4)|...
        // -> mask
        // 11111111|00000000|00000000|11111111|...
        let flagv_2_32_ext = _mm256_and_si256(flagv_2_32, bit_ext_mask);
        let flagv_1_32_ext = _mm256_and_si256(flagv_1_32, bit_ext_mask);
        let flagv_0_32_ext = _mm256_and_si256(flagv_0_32, bit_ext_mask);

        let take_mask_2_32 = _mm256_cmpeq_epi8(flagv_2_32_ext, bit_ext_mask);
        let take_mask_1_32 = _mm256_cmpeq_epi8(flagv_1_32_ext, bit_ext_mask);
        let take_mask_0_32 = _mm256_cmpeq_epi8(flagv_0_32_ext, bit_ext_mask);

        // bi=0-3, shift=24,16,8,0
        //const SHIFTS: &'static [i32] = &[24, 16, 8, 0];

        for bi in 0usize..4 {
            // DCBADCBA...DCBA
            // -> b=1: for B
            // BA00BA00...BA00

            let shift_v = shifts[bi];
            let take_mask_2 = _mm256_sllv_epi32(take_mask_2_32, shift_v);
            let take_mask_1 = _mm256_sllv_epi32(take_mask_1_32, shift_v);
            let take_mask_0 = _mm256_sllv_epi32(take_mask_0_32, shift_v);

            //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
            //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
            //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
            //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);

            let wlsv_lo_ptr = wls_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let wlsv_hi_ptr = wls_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            let zsv_lo_ptr = zs_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let zsv_hi_ptr = zs_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            //let psv_lo_ptr = ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            //let psv_hi_ptr = ps[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            let wlsv_lo: __m256d = _mm256_load_pd(wlsv_lo_ptr as *const _);
            let wlsv_hi: __m256d = _mm256_load_pd(wlsv_hi_ptr as *const _);
            let zsv_lo: __m256d = _mm256_load_pd(zsv_lo_ptr as *const _);
            let zsv_hi: __m256d = _mm256_load_pd(zsv_hi_ptr as *const _);

            //log::debug!("ps lo {:?}", psv_lo);
            //log::debug!("ps hi {:?}", psv_hi);

            // to create f_lo, first set score0 and score1, and then set score2
            //let f_0_lo: __m256d= _mm256_blendv_pd(s0d, s1d, _mm256_castsi256_pd(take_mask_1));
            // missing=0.0
            // TODO: if missing!=0.0, make zerod -> smd
            let f_0_lo: __m256d = _mm256_blendv_pd(zerod, s0d, _mm256_castsi256_pd(take_mask_0));
            let f_01_lo: __m256d = _mm256_blendv_pd(f_0_lo, s1d, _mm256_castsi256_pd(take_mask_1));
            let f_lo: __m256d = _mm256_blendv_pd(f_01_lo, s2d, _mm256_castsi256_pd(take_mask_2));
            // create f-z
            let f_z_lo: __m256d = _mm256_sub_pd(f_lo, zsv_lo);
            // (f-z)**2
            let f_z_2_lo: __m256d = _mm256_mul_pd(f_z_lo, f_z_lo);
            // w * (f-z)**2
            let w_f_z_2_lo: __m256d = _mm256_mul_pd(wlsv_lo, f_z_2_lo);

            l_acc = _mm256_add_pd(l_acc, w_f_z_2_lo);

            // for high
            let take_mask_2_hi = _mm256_slli_epi64(take_mask_2, 32);
            let take_mask_1_hi = _mm256_slli_epi64(take_mask_1, 32);
            let take_mask_0_hi = _mm256_slli_epi64(take_mask_0, 32);

            let f_0_hi: __m256d = _mm256_blendv_pd(zerod, s0d, _mm256_castsi256_pd(take_mask_0_hi));
            let f_01_hi: __m256d =
                _mm256_blendv_pd(f_0_hi, s1d, _mm256_castsi256_pd(take_mask_1_hi));
            //let f_01_hi: __m256d= _mm256_blendv_pd(s0d, s1d, _mm256_castsi256_pd(take_mask_1_hi));
            let f_hi: __m256d = _mm256_blendv_pd(f_01_hi, s2d, _mm256_castsi256_pd(take_mask_2_hi));
            // create f-z
            let f_z_hi: __m256d = _mm256_sub_pd(f_hi, zsv_hi);
            // (f-z)**2
            let f_z_2_hi: __m256d = _mm256_mul_pd(f_z_hi, f_z_hi);
            // w * (f-z)**2
            let w_f_z_2_hi: __m256d = _mm256_mul_pd(wlsv_hi, f_z_2_hi);

            l_acc = _mm256_add_pd(l_acc, w_f_z_2_hi);
        }
    }

    // sum 4 double horizontally to get the whole sum
    l_acc = _mm256_hadd_pd(l_acc, l_acc);

    // 1. any way to hadd??
    // 2. _mm256_extractf128_pd and _mm256_cvtsd_f64: get 64:0

    let loss_least_square: f64 =
        _mm256_cvtsd_f64(l_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(l_acc, 1));

    let loss: f64;
    if use_adjloss {
        loss = 0.5 * (loss_least_square - loss_max_theory);
    } else {
        loss = loss_least_square;
    }
    //let loss: f64 = 0.5 * (loss_least_square - loss_max_theory);

    return loss;
}

/// only difference from logit is
/// 1. coef is linearconst
/// 2. assume no missing
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn calculate_loss_gt_logit_add_simd_sm(
    predicts: &GenotSnvRef,
    coef: &Coef,
    wls: &[f64],
    zs: &[f64],
    use_adjloss: bool,
    loss_max_theory: f64,
) -> f64 {
    let (constt, alphat) = coef.linearconst_f64();

    if constt.is_nan() || alphat.is_nan() {
        //log::debug!("coef is nan {}, {}", c, a);
        return f64::NAN;
    };

    //let (s0, s1, s2, sm) = coef.score4_f64();
    let s0 = constt;
    let s1 = constt + alphat;
    let s2 = constt + 2.0 * alphat;

    //if sm.abs() > 1e-10 {
    //    panic!("sm should be 0.0");
    //}

    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    //use std::convert::TryInto;

    //let n = phe.n();
    //let n_pad=wls.len();
    let n = wls.len() - 32;
    //let ys = phe.inner();
    let pred_s0m = predicts.predict_s(0);
    let pred_s1m = predicts.predict_s(1);

    let bit_ext_mask: __m256i = _mm256_set_epi32(
        0x10101010,
        0x01010101,
        0x20202020,
        0x02020202,
        0x40404040,
        0x04040404,
        -0x7f7f7f80, // =0x80808080
        0x08080808,
    );
    //let zerod: __m256d = _mm256_setzero_pd();
    //let zeros: __m256i = _mm256_setzero_si256();
    //let ones: __m256i = _mm256_cmpeq_epi32(zeros, zeros);

    let s2d: __m256d = _mm256_set1_pd(s2);
    let s1d: __m256d = _mm256_set1_pd(s1);
    let s0d: __m256d = _mm256_set1_pd(s0);

    let mut l_acc: __m256d = _mm256_setzero_pd();

    //let mut wzs_sum2_acc = _mm256_setzero_pd();
    //let mut wzs_sum1_acc = _mm256_setzero_pd();
    //let mut wzs_sum0_acc = _mm256_setzero_pd();
    //let mut wls_sum2_acc = _mm256_setzero_pd();
    //let mut wls_sum1_acc = _mm256_setzero_pd();
    //let mut wls_sum0_acc = _mm256_setzero_pd();

    let shifts: [__m256i; 4] = [
        _mm256_set1_epi32(24),
        _mm256_set1_epi32(16),
        _mm256_set1_epi32(8),
        _mm256_set1_epi32(0),
    ];

    for ni in 0..(n / 32 + 1) {
        //log::debug!("ni {}", ni);

        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let yv_32 = _mm256_set1_epi32(ys_b32 as i32);

        // bitwise not is fastest in xor by xor(v, ones)
        // sum2=( pred0 & pred1 )
        // sum1= ( !pred1 & pred0 )
        // sum0= ( !pred0 & (!pred1) )
        let flagv_2_32 = _mm256_and_si256(predv_s0_32, predv_s1_32);
        let flagv_1_32 = _mm256_andnot_si256(predv_s1_32, predv_s0_32);
        //let flagv_0_32 = _mm256_andnot_si256(predv_s0_32, _mm256_xor_si256(predv_s1_32, ones));

        // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
        // abcdefgh(x4)|...
        // -> extracted (at highest position)
        // 000d0000(x4)|...
        // -> mask
        // 11111111|00000000|00000000|11111111|...
        let flagv_2_32_ext = _mm256_and_si256(flagv_2_32, bit_ext_mask);
        let flagv_1_32_ext = _mm256_and_si256(flagv_1_32, bit_ext_mask);
        //let flagv_0_32_ext = _mm256_and_si256(flagv_0_32, bit_ext_mask);

        let take_mask_2_32 = _mm256_cmpeq_epi8(flagv_2_32_ext, bit_ext_mask);
        let take_mask_1_32 = _mm256_cmpeq_epi8(flagv_1_32_ext, bit_ext_mask);
        //let take_mask_0_32 = _mm256_cmpeq_epi8(flagv_0_32_ext, bit_ext_mask);

        // bi=0-3, shift=24,16,8,0
        //const SHIFTS: &'static [i32] = &[24, 16, 8, 0];

        for bi in 0usize..4 {
            // DCBADCBA...DCBA
            // -> b=1: for B
            // BA00BA00...BA00

            let shift_v = shifts[bi];
            let take_mask_2 = _mm256_sllv_epi32(take_mask_2_32, shift_v);
            let take_mask_1 = _mm256_sllv_epi32(take_mask_1_32, shift_v);
            //let take_mask_0 = _mm256_sllv_epi32(take_mask_0_32, shift_v);

            //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
            //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
            //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
            //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);

            let wlsv_lo_ptr = wls[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let wlsv_hi_ptr = wls[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            let zsv_lo_ptr = zs[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let zsv_hi_ptr = zs[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            //let psv_lo_ptr = ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            //let psv_hi_ptr = ps[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            let wlsv_lo: __m256d = _mm256_load_pd(wlsv_lo_ptr as *const _);
            let wlsv_hi: __m256d = _mm256_load_pd(wlsv_hi_ptr as *const _);
            let zsv_lo: __m256d = _mm256_load_pd(zsv_lo_ptr as *const _);
            let zsv_hi: __m256d = _mm256_load_pd(zsv_hi_ptr as *const _);

            //log::debug!("ps lo {:?}", psv_lo);
            //log::debug!("ps hi {:?}", psv_hi);

            // to create f_lo, first set score0 and score1, and then set score2
            let f_01_lo: __m256d = _mm256_blendv_pd(s0d, s1d, _mm256_castsi256_pd(take_mask_1));
            // no missing here
            //let f_0_lo: __m256d = _mm256_blendv_pd(zerod, s0d, _mm256_castsi256_pd(take_mask_0));
            //let f_01_lo: __m256d = _mm256_blendv_pd(f_0_lo, s1d, _mm256_castsi256_pd(take_mask_1));
            let f_lo: __m256d = _mm256_blendv_pd(f_01_lo, s2d, _mm256_castsi256_pd(take_mask_2));
            // create f-z
            let f_z_lo: __m256d = _mm256_sub_pd(f_lo, zsv_lo);
            // (f-z)**2
            let f_z_2_lo: __m256d = _mm256_mul_pd(f_z_lo, f_z_lo);
            // w * (f-z)**2
            let w_f_z_2_lo: __m256d = _mm256_mul_pd(wlsv_lo, f_z_2_lo);

            l_acc = _mm256_add_pd(l_acc, w_f_z_2_lo);

            // for high
            let take_mask_2_hi = _mm256_slli_epi64(take_mask_2, 32);
            let take_mask_1_hi = _mm256_slli_epi64(take_mask_1, 32);
            //let take_mask_0_hi = _mm256_slli_epi64(take_mask_0, 32);

            let f_01_hi: __m256d = _mm256_blendv_pd(s0d, s1d, _mm256_castsi256_pd(take_mask_1_hi));
            //let f_0_hi: __m256d = _mm256_blendv_pd(zerod, s0d, _mm256_castsi256_pd(take_mask_0_hi));
            //let f_01_hi: __m256d =
            //_mm256_blendv_pd(f_0_hi, s1d, _mm256_castsi256_pd(take_mask_1_hi));
            //let f_01_hi: __m256d= _mm256_blendv_pd(s0d, s1d, _mm256_castsi256_pd(take_mask_1_hi));
            let f_hi: __m256d = _mm256_blendv_pd(f_01_hi, s2d, _mm256_castsi256_pd(take_mask_2_hi));
            // create f-z
            let f_z_hi: __m256d = _mm256_sub_pd(f_hi, zsv_hi);
            // (f-z)**2
            let f_z_2_hi: __m256d = _mm256_mul_pd(f_z_hi, f_z_hi);
            // w * (f-z)**2
            let w_f_z_2_hi: __m256d = _mm256_mul_pd(wlsv_hi, f_z_2_hi);

            l_acc = _mm256_add_pd(l_acc, w_f_z_2_hi);
        }
    }

    // sum 4 double horizontally to get the whole sum
    l_acc = _mm256_hadd_pd(l_acc, l_acc);

    // 1. any way to hadd??
    // 2. _mm256_extractf128_pd and _mm256_cvtsd_f64: get 64:0

    let loss_least_square: f64 =
        _mm256_cvtsd_f64(l_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(l_acc, 1));

    let loss: f64;
    if use_adjloss {
        loss = 0.5 * (loss_least_square - loss_max_theory);
    } else {
        loss = loss_least_square;
    }
    //let loss: f64 = 0.5 * (loss_least_square - loss_max_theory);

    return loss;
}

//fn compute_loss_least_square_max_theory(zs: &[f64], wls: &[f64], n: usize) -> f64 {
fn compute_loss_least_square_max_theory(zs: &[f64], wls: &[f64]) -> f64 {
    zs.iter().zip(wls.iter()).map(|(z, w)| w * z * z).sum()
}

// unsafe is necessary with #[target_feature]
//pub unsafe fn calculate_loss_gt_freemodelmissing_simd(
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn calculate_loss_gt_logit_simd(
    losss: &mut LossStruct,
    genot: &Genot,
    sample_weight: &SampleWeight,
    //wzs_pad: &[f64],
    //wls_pad: &[f64],
    //zs_pad: &[f64],
    phe: &Phe,
    boost_param: BoostParam,
    skip_snv: &HashSet<usize>,
    use_adjloss: bool,
) {
    let epsilons_wzs =
        epsilon::calculate_epsilons_logit_wzs(sample_weight.wzs().unwrap(), phe, boost_param.eps());
    //let n=phe.n();
    //let epsilons_wzs = epsilon::calculate_epsilons_logit_wzs(&wzs_pad[..n], phe, boost_param.eps());
    //let epsilons_wzs = epsilon::calculate_epsilons_logit_wzs(wzs_pad, phe, boost_param.eps());
    log::debug!(
        "epsilon_wzs case, cont: {:.4e},{:.4e}",
        epsilons_wzs.0,
        epsilons_wzs.1
    );

    let epsilons_wls =
        epsilon::calculate_epsilons_logit_wls(sample_weight.wls().unwrap(), phe, boost_param.eps());
    //let epsilons_wls = epsilon::calculate_epsilons_logit_wls(&wls_pad[..n], phe, boost_param.eps());
    //let epsilons_wls = epsilon::calculate_epsilons_logit_wls(wls_pad, phe, boost_param.eps());
    log::debug!(
        "epsilon_wls case, cont: {:.4e},{:.4e}",
        epsilons_wls.0,
        epsilons_wls.1
    );

    assert_eq!(losss.inner_mut().len(), genot.m());

    let loss_max_theory = compute_loss_least_square_max_theory(
        sample_weight.zs().unwrap(),
        sample_weight.wls().unwrap(),
    );
    //let loss_max_theory = compute_loss_least_square_max_theory(&zs_pad[..n], &wls_pad[..n]);
    log::debug!("loss_max_theory {}", loss_max_theory);

    //unsafe {
    losss
        .inner_mut()
        .par_iter_mut()
        .enumerate()
        .for_each(|(mi, loss)| {
            if skip_snv.contains(&mi) {
                *loss = f64::MAX;
            } else {
                let (coef, _, _) = coefficient::calculate_coef_logit_eps(
                    &genot.to_genot_snv(mi),
                    sample_weight.wzs_pad().unwrap(),
                    sample_weight.wls_pad().unwrap(),
                    //wzs_pad,
                    //wls_pad,
                    phe,
                    epsilons_wzs,
                    epsilons_wls,
                    boost_param.eps(),
                    //boost_param.learning_rate(),
                    boost_param.boost_type(),
                    boost_param.eff_eps(),
                    true,
                    false,
                );
                *loss = calculate_loss_gt_logit_simd_sm(
                    &genot.to_genot_snv(mi),
                    &coef,
                    sample_weight.wls_pad().unwrap(),
                    sample_weight.zs_pad().unwrap(),
                    //wls_pad,
                    //zs_pad,
                    use_adjloss,
                    loss_max_theory,
                );
            }
        });
    //}
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn calculate_loss_gt_logit_add_simd(
    losss: &mut LossStruct,
    genot: &Genot,
    sample_weight: &SampleWeight,
    //wzs_pad: &[f64],
    //wls_pad: &[f64],
    //zs_pad: &[f64],
    //phe: &Phe,
    _boost_param: BoostParam,
    skip_snv: &HashSet<usize>,
    use_adjloss: bool,
) {
    assert_eq!(losss.inner_mut().len(), genot.m());

    let loss_max_theory = compute_loss_least_square_max_theory(
        sample_weight.zs().unwrap(),
        sample_weight.wls().unwrap(),
    );
    //let n = genot.n();
    //let loss_max_theory = compute_loss_least_square_max_theory(&zs_pad[..n], &wls_pad[..n]);
    log::debug!("loss_max_theory {}", loss_max_theory);

    //unsafe {
    losss
        .inner_mut()
        .par_iter_mut()
        .enumerate()
        .for_each(|(mi, loss)| {
            if skip_snv.contains(&mi) {
                *loss = f64::MAX;
            } else {
                let coef = coefficient::calculate_coef_logit_add_on_loss(
                    &genot.to_genot_snv(mi),
                    sample_weight.wzs_pad().unwrap(),
                    sample_weight.wls_pad().unwrap(),
                    /////////////////////////
                    /////////////////////// TMP
                    //None,
                    //boost_param.learning_rate(),
                );
                // debug
                //let (c, a) = coef.linearconst_f64();
                //if c.is_nan() || a.is_nan() {
                //    log::debug!("coef is nan {}, {}", c, a);
                //};
                *loss = calculate_loss_gt_logit_add_simd_sm(
                    &genot.to_genot_snv(mi),
                    &coef,
                    sample_weight.wls_pad().unwrap(),
                    sample_weight.zs_pad().unwrap(),
                    //wls_pad,
                    //zs_pad,
                    use_adjloss,
                    loss_max_theory,
                );
                //if loss.is_nan() {
                //    log::debug!("loss is nan {}", loss);
                //};
            }
        });
    //}
}

/* pub fn calculate_loss_gt_freemodelmissing_simd(
    losss: &mut [f64],
    predictions: &Genot,
    //predictions: &GenotBi<Vec<u8>>,
    //predictions: &[B8],
    ps: &[f64],
    ys: &Phe,
    //m: usize,
    boost_param: BoostParam,
) {
    //unimplemented!();

    let epsilons = epsilon::calculate_epsilons(ps, ys, boost_param.eps());
    //log::debug!("epsilon case, cont: {:.2e},{:.2e}", epsilons.0, epsilons.1);

    // first calculate E:=case_ws_sum and F:=control_ws_sum
    //let ef_ = calculate_ef_sum(ps, ys, n);

    fn loss_iter(
        loss: &mut f64,
        predicts: &GenotSnvRef,
        //predicts: &[B8],
        ps: &[f64],
        ys: &Phe,
        epsilons: (f64, f64),
        eps: Eps,
    ) {
        //log::debug!("predicts len {}", predicts.len());
        unsafe {
            let loss_ = calculate_loss_gt_freemodelmissing_simd_sm(predicts, ps, ys, epsilons, eps);
            *loss = loss_
        }
    }

    //losss.par_chunks_mut(2).enumerate().for_each(|(mi, loss)| {
    losss.par_iter_mut().enumerate().for_each(|(mi, loss)| {
        //log::debug!("mi {}", mi);
        loss_iter(
            loss,
            &predictions.to_genot_snv(mi),
            //predictions.predictions_snv_s(mi),
            //predict::predictions_snv_s(predictions, mi, n),
            //predict::get_predictions_snv(predictions, 0, mi, n),
            //predict::get_predictions_snv(predictions, 1, mi, n),
            ps,
            ys,
            epsilons,
            boost_param.eps(),
        )
    });
}
 */

fn calculate_loss_gt_freemodelmissing_nosimd_sm(
    pred_s: &[u8],
    //pred_s: &GenotSnvRef,
    //pred_s: &[B8],
    ps: &[f64],
    phe: &Phe,
    epsilons: (f64, f64),
    eps: Option<Eps>,
) -> f64 {
    //let pred = pred_s.vals();
    let (table, _) = table::calculate_table7_epsilons(&pred_s, ps, phe, epsilons, eps);

    calculate_loss_table7(table)
}

// TODO: this can be integrated to freemodelmissing_simd()
pub fn calculate_loss_gt_freemodelmissing_nosimd(
    losss: &mut LossStruct,
    genot: &Genot,
    //predictions: &GenotBi<Vec<u8>>,
    //predictions: &[B8],
    //ps: &[f64],
    sample_weight: &SampleWeight,
    phe: &Phe,
    boost_param: BoostParam,
) {
    assert_eq!(losss.inner_mut().len(), genot.m());

    let epsilons = epsilon::calculate_epsilons(sample_weight.ps().unwrap(), phe, boost_param.eps());
    //log::debug!("epsilon case, cont: {:.2e},{:.2e}", epsilons.0, epsilons.1);

    losss
        .inner_mut()
        .par_iter_mut()
        .enumerate()
        .for_each(|(mi, loss)| {
            *loss = calculate_loss_gt_freemodelmissing_nosimd_sm(
                &genot.vals_snv(mi),
                //&predictions.to_genot_snv(mi).vals(),
                //&predictions.vals_snv(mi),
                sample_weight.ps().unwrap(),
                //ps,
                phe,
                epsilons,
                boost_param.eps(),
            );
        })
    //log::debug!("losss {:?}", losss);
    //losss
}

// TODO: since I have integrated loss_gt_simd and loss_gt_nosimd into one func. many tests cannot be tested...
// TODO; should integrate setup_test, setup2_test...
#[cfg(test)]
mod tests {
    use super::*;

    use crate::boosting_train::sample_weight::SampleWeight;
    use genetics::alloc;

    //fn is_eq_f64(v: f64, w: f64, e: f64) -> bool {
    //    (v - w).abs() < e
    //}

    #[test]
    fn test_cast() {
        // 192u8
        let a: u8 = 0b1100_0000;
        //log::debug!("a: {}", a);
        //println!("a: {:#010b}", a);
        // -64i8
        let b: i8 = a as i8;
        //println!("b: {}", b);
        //println!("b: {:#010b}", b);
        assert_eq!(format!("{:b}", a), format!("{:b}", b));
        assert_eq!(a as u16 + b.abs() as u16, 256u16);

        let c: i32 = -0x7f7f7f80;
        //println!("c: {}", c);
        //println!("c: {:#034b}", c);
        let d: u32 = 0x80808080;
        //println!("d: {}", d);
        //println!("d: {:#034b}", d);
        // both are the same
        assert_eq!(format!("{:b}", c), format!("{:b}", d));
    }

    // using outside of file should be avoided here
    /*     fn load_dataset(fin: &str) -> (Dataset, Vec<f64>, LossStruct) {
           let dataset: Dataset = Dataset::new(&fin, None, None, None, false);
           let n = dataset.samples().samples_n();
           let m = dataset.snvs().snvs_n();

           let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
           ps.resize(n + 32, 1.0f64 / (n as f64));
           //let ps: Vec<f64> = vec![1.0f64 / (n as f64); n + 64];
           //let ps: Vec<f64> = vec![1.0f64 / (n as f64); n + 32];

           //let losss = vec![0.0f64; 2 * m];
           let losss = LossStruct::new(crate::BoostType::ConstAda, m);

           (dataset, ps, losss)
       }

       //fn setup_test2() -> (DatasetTwin, Vec<f64>, Vec<f64>) {
       fn setup_test2() -> (Dataset, Vec<f64>, LossStruct) {
           let fin = String::from("../../test/data/1kg_n10000/1kg_n10000");
           load_dataset(&fin)
       }

       fn setup_test3() -> (Dataset, Vec<f64>, LossStruct) {
           let fin = String::from("../../test/data/1kg_n100000/1kg_n100000");
           load_dataset(&fin)
       }
    */
    /*     #[test]
       fn test_calculate_loss_gt_simd2() {
           let (dataset, ps, mut losss) = setup_test2();

           println!("done setup");

           assert_eq!(losss.inner().len(), dataset.snvs().snvs_n() * 2);

           calculate_loss_gt_constada_simd(
               &mut losss,
               &dataset.genot(),
               &ps,
               &dataset.samples().phe(),
               // ConstAda
               BoostParam::new_type1(),
           );
       }
    */
    /*     // 1. using Dataset, Phe raised error (above)
    // 2. using inner
    //      -> ng
    // 3. using new vec (not aligned ) with the same length
    // Oh, ps not aligned
    #[test]
    fn test_calculate_loss_gt_simd_old() {
        let (dataset, ps, mut losss) = setup_test2();

        println!("done setup");

        assert_eq!(losss.inner().len(), dataset.snvs().snvs_n() * 2);

        //let ys = vec![0u8; dataset.samples().phe().inner().len()];
        let n = dataset.genot().n();
        let len_n = n / 8 + 5;
        let ys = vec![0u8; len_n];
        //let mut ys: Vec<B8> = alloc::with_capacity_align_u8(len_n);
        //ys.resize(len_n, 0x00);

        let len_g = dataset.genot().genot_inner().inner().len();
        let g = vec![0u8; len_g];
        //let mut g: Vec<B8> = alloc::with_capacity_align_u8(len_g);
        //g.resize(len_g, 0x00);
        let m = dataset.genot().m();
        assert_eq!(len_g, len_n * m * 2);

        calculate_loss_gt_simd_old(
            &mut losss.inner(),
            &g,
            //&dataset.genot().genot_inner().inner(),
            &ps,
            &ys,
            //&dataset.samples().phe().inner(),
            // ConstAda
            dataset.genot().n(),
            BoostParam::new_type1(),
        );
    } */

    /*     #[test]
       fn test_calculate_loss_gt_nosimd2() {
           let (dataset, ps, mut losss) = setup_test2();

           calculate_loss_gt_constada_simd(
               &mut losss,
               &dataset.genot(),
               &ps,
               &dataset.samples().phe(),
               BoostParam::new_type1(),
           );
       }

       #[test]
       fn test_calculate_loss_gt_simd3() {
           let (dataset, ps, mut losss) = setup_test3();

           calculate_loss_gt_constada_simd(
               &mut losss,
               &dataset.genot(),
               &ps,
               &dataset.samples().phe(),
               BoostParam::new_type1(),
           );
       }

       #[test]
       fn test_calculate_loss_gt_nosimd3() {
           let (dataset, ps, mut losss) = setup_test3();

           calculate_loss_gt_constada_nosimd(
               &mut losss,
               &dataset.genot(),
               &ps,
               &dataset.samples().phe(),
               BoostParam::new_type1(),
           );
       }
    */

    // TODO: larger test
    // with eps
    //fn setup_test2() -> (Genot, Phe, Vec<f64>) {
    fn setup_test2() -> (Genot, Phe, SampleWeight) {
        let x = vec![0u8, 1, 1, 2];
        let n = x.len();
        let genot = Genot::new(1, n, &x);
        let y = vec![false, true, false, true];
        let phe = Phe::new(&y);

        let mut ps_pad: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        for p in [0.1, 0.2, 0.3, 0.4].iter() {
            ps_pad.push(*p);
        }
        for _ in 0..32 {
            ps_pad.push(0.0);
        }
        // ps:  [0.1,0.2,0.3,0.4]
        // ys:  [0,1,0,1]
        // predict
        // dom: [0,1,1,1]
        // rec: [0,0,0,1]
        //
        // abcd
        // dom: [0.6,0.3,0.0,0.1]
        // rec: [0.4,0.0,0.2,0.4]

        let sw = SampleWeight::_new_test(n, ps_pad);

        //let sw = SampleWeight {
        //    n,
        //    len_n: n + 32,
        //    ps_pad: Some(ps_pad),
        //    ..Default::default()
        //};
        (genot, phe, sw)
        //(genot, phe, ps_pad)
    }

    // with eps
    //fn setup_test3() -> (Genot, Phe, Vec<f64>) {
    fn setup_test3() -> (Genot, Phe, SampleWeight) {
        let x = vec![0u8, 0, 0, 0, 1, 1, 1, 2, 2, 3];
        let n = x.len();
        let genot = Genot::new(1, n, &x);
        let y = vec![
            true, true, false, false, true, false, false, false, false, true,
        ];
        let phe = Phe::new(&y);

        let mut ps_pad: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        for _ in 0..n {
            ps_pad.push(0.1);
        }
        for _ in 0..32 {
            ps_pad.push(0.0);
        }
        // ps:  [0.1] * 10
        // ys:
        // table7: [0.0,0.2,0.1,0.2,0.2,0.2,0.1]

        let sw = SampleWeight::_new_test(n, ps_pad);

        (genot, phe, sw)
        //(genot, phe, ps)
    }

    #[test]
    fn test_calculate_loss_gt_constada_simd_2() {
        // This error below means SIMD memory is not aligned.
        // "process didn't exit successfully: (signal: 11, SIGSEGV: invalid memory reference)"

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                let (genot, phe, sample_weight) = setup_test2();

                let m = genot.m();
                let losss = vec![0.0f64; 2 * m];
                let mut losss = LossStruct::ConstAda(losss, m, 2);

                unsafe {
                    calculate_loss_gt_constada_simd(
                        &mut losss,
                        &genot,
                        &sample_weight,
                        //&ps,
                        &phe,
                        BoostParam::new_type1(),
                        &HashSet::<usize>::new(),
                    );
                }
                //calculate_loss_gt_simd(&mut losss, &mistakes, &ps, &ys, m, n);
                //println!("losss: {:?}", losss);

                // abcd
                // ps:  [0.1,0.2,0.3,0.4]
                // ys:  [0,1,0,1]
                // eps: (0.2, 0.2)
                // dom: [0.6,0.3,0.0,0.1]
                //      -> [0.8,0.5,0.2,0.3]
                // rec: [0.4,0.0,0.2,0.4]
                //      ->[0.6,0.2,0.4,0.6]
                // dom: 1.7548090
                // rec: 1.67261622
                // -> not realistic due to small n
                assert_float_absolute_eq!(losss.inner_mut()[0], 1.7548090);
                assert_float_absolute_eq!(losss.inner_mut()[0], 1.7548090);
                //assert!(is_eq_f64(losss.inner_mut()[1], 1.67261622, 1e-7));
                //assert!(is_eq_f64(losss.inner_mut()[0], 1.7548090, 1e-7));
                //assert!(is_eq_f64(losss.inner_mut()[1], 1.67261622, 1e-7));
            }
        }
    }

    #[test]
    fn test_calculate_loss_gt_constada_nosimd_2() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                let (genot, phe, sw) = setup_test2();

                let m = genot.m();
                let losss = vec![0.0f64; 2 * m];
                let mut losss = LossStruct::ConstAda(losss, m, 2);
                unsafe {
                    calculate_loss_gt_constada_simd(
                        &mut losss,
                        &genot,
                        &sw,
                        &phe,
                        BoostParam::new_type1(),
                        &HashSet::<usize>::new(),
                    );
                }

                let losss_nosimd = vec![0.0f64; 2 * m];
                let mut losss_nosimd = LossStruct::ConstAda(losss_nosimd, m, 2);
                calculate_loss_gt_constada_nosimd(
                    &mut losss_nosimd,
                    &genot,
                    &sw,
                    &phe,
                    BoostParam::new_type1(),
                );

                assert_float_absolute_eq!(losss_nosimd.inner_mut()[0], losss.inner_mut()[0]);
                assert_float_absolute_eq!(losss_nosimd.inner_mut()[1], losss.inner_mut()[1]);
                //assert!(is_eq_f64(
                //    losss_nosimd.inner_mut()[0],
                //    losss.inner_mut()[0],
                //    1e-7
                //));
                //assert!(is_eq_f64(
                //    losss_nosimd.inner_mut()[1],
                //    losss.inner_mut()[1],
                //    1e-7
                //));
            }
        }
    }

    #[test]
    fn test_calculate_loss_gt_freemodelmissing_simd_3() {
        // This error below is due to that SIMD memory is not aligned.
        // "process didn't exit successfully: (signal: 11, SIGSEGV: invalid memory reference)"

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                let (genot, phe, ps) = setup_test3();

                let m = genot.m();
                let losss = vec![0.0f64; m];
                let mut losss = LossStruct::LossOne(losss, m);

                unsafe {
                    calculate_loss_gt_freemodelmissing_simd(
                        &mut losss,
                        &genot,
                        &ps,
                        &phe,
                        BoostParam::new_type1(),
                        &HashSet::<usize>::new(),
                    );
                }

                // abcd
                // eps: (0.05, 0.05)
                // table7: [0.0,0.2,0.1,0.2,0.2,0.2,0.1]
                //         ->[0.05,0.25,0.15,0.25,0.25,0.25,0.1]
                // -> not realistic due to small n
                assert_float_absolute_eq!(losss.inner_mut()[0], 1.21090513);
                //assert!(is_eq_f64(losss.inner_mut()[0], 1.21090513, 1e-7));
            }
        }
    }

    #[test]
    fn test_calculate_loss_gt_freemodelmissing_nosimd_3() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                let (genot, phe, sw) = setup_test2();

                let m = genot.m();
                let losss = vec![0.0f64; m];
                let mut losss = LossStruct::LossOne(losss, m);
                unsafe {
                    calculate_loss_gt_freemodelmissing_simd(
                        &mut losss,
                        &genot,
                        &sw,
                        &phe,
                        BoostParam::new_type1(),
                        &HashSet::<usize>::new(),
                    );
                }

                let losss_nosimd = vec![0.0f64; m];
                let mut losss_nosimd = LossStruct::LossOne(losss_nosimd, m);
                calculate_loss_gt_freemodelmissing_nosimd(
                    &mut losss_nosimd,
                    &genot,
                    &sw,
                    &phe,
                    BoostParam::new_type1(),
                );

                assert_float_absolute_eq!(losss_nosimd.inner_mut()[0], losss.inner_mut()[0]);

                //assert!(is_eq_f64(
                //    losss_nosimd.inner_mut()[0],
                //    losss.inner_mut()[0],
                //    1e-7
                //));
            }
        }
    }
}
