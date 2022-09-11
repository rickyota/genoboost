use super::epsilon;
use super::table;
use super::LossStruct;
use crate::{BoostParam, ContingencyTable, Eps};
use genetics::genot::prelude::*;
use genetics::samples::prelude::*;
use rayon::prelude::*;

pub fn calculate_loss_table4(table: ContingencyTable) -> f64 {
    let (a, b, c, d) = table.four();
    //let (a, b, c, d) = table.four();
    2.0 * ((a * b).sqrt() + (c * d).sqrt())
}

pub fn calculate_loss_table7(table: ContingencyTable) -> f64 {
    let (d2, n2, d1, n1, d0, n0, m) = table.seven();
    m + 2.0 * ((d2 * n2).sqrt() + (d1 * n1).sqrt() + (d0 * n0).sqrt())
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
    predicts: &GenotSnvRef,
    ps: &[f64],
    phe: &Phe,
    ef_: (f64, f64),
    epsilons: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    eps: Eps,
) -> (f64, f64) {
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
    //println!("0x7f7f7f7f + 1: {:#010x}", 0x7f7f7f7f + 1 as i32);
    //println!("0x7f7f7f7f + 1: {}", 0x7f7f7f7f + 1 as i32);
    //println!("0x7f7f7f80: {}", 0x7f7f7f80 as i32);

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
        //println!("ni {}", ni);

        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        // 1. use _mm256_set_epi8
        // o2. use from_be() and use set1 <- assume to be fast??
        let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        //println!("mis 0: {:#010b}", mis_s0m[4 * ni]);
        //println!("mis 1: {:#010b}", mis_s0m[4 * ni + 1]);
        //println!("mis 2: {:#010b}", mis_s0m[4 * ni + 2]);
        //println!("mis 3: {:#010b}", mis_s0m[4 * ni + 3]);
        //println!("mis b32: {:#034b}", mis_s0_b32);
        //println!("misv {:?}", misv_s0_32);

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

            //println!("take_mask a s0 {:?}", take_mask_a_s0);
            //println!("take_mask b s0 {:?}", take_mask_b_s0);
            //println!("take_mask a s1 {:?}", take_mask_a_s1);
            //println!("take_mask b s1 {:?}", take_mask_b_s1);

            // ps[i].as_ptr() should be enough
            let psv_lo_ptr = ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let psv_hi_ptr = ps[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();

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

            //println!("ps a s0 lo {:?}", ps_masked_a_s0_lo);
            //println!("ps a s1 lo {:?}", ps_masked_a_s1_lo);
            //println!("ps b s0 lo {:?}", ps_masked_b_s0_lo);
            //println!("ps b s1 lo {:?}", ps_masked_b_s1_lo);

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

            //println!("a s0 hi {:?}", ps_masked_a_s0_hi);

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

    let (abcd_s0, _) = table::fixed_abcd_sum((a_sum_s0, b_sum_s0), ef_, epsilons, eps);
    let (abcd_s1, _) = table::fixed_abcd_sum((a_sum_s1, b_sum_s1), ef_, epsilons, eps);

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
    ps: &[f64],
    phe: &Phe,
    boost_param: BoostParam,
) {
    assert_eq!(losss.inner().len(), genot.m() * 2);
    assert_eq!(phe.n(), genot.n());

    let epsilons = epsilon::calculate_epsilons(ps, phe, boost_param.eps());
    println!("epsilon case, cont: {:.4e},{:.4e}", epsilons.0, epsilons.1);

    // first calculate E:=case_ws_sum and F:=control_ws_sum
    let ef_ = table::calculate_ef_sum(ps, phe);

    losss
        .inner()
        .par_chunks_mut(2)
        .enumerate()
        .for_each(|(mi, loss)| {
            let loss_ = calculate_loss_gt_constada_simd_sm(
                &genot.to_genot_snv(mi),
                ps,
                phe,
                ef_,
                epsilons,
                boost_param.eps(),
            );
            loss[0] = loss_.0;
            loss[1] = loss_.1;
        });
}

fn calculate_loss_gt_constada_nosimd_sm(
    //pred_s: &GenotSnvRef,
    pred_sm: &[u8],
    //pred_sm: &[B8],
    ps: &[f64],
    ys: &Phe,
    ef_: (f64, f64),
    epsilons: (f64, f64),
    eps: Eps,
) -> f64 {
    let (abcd, _) = table::calculate_abcd_sum_sm(pred_sm, ps, ys, ef_, epsilons, eps);

    calculate_loss_table4(abcd)
}

// TODO: use rayon::par_chunks_mut ?
/// losss[..,2m,2mi+1..] 2mi: dom, 2mi+1: rec
pub fn calculate_loss_gt_constada_nosimd(
    losss: &mut LossStruct,
    predictions: &Genot,
    //predictions: &GenotBi<Vec<u8>>,
    //predictions: &[B8],
    ps: &[f64],
    ys: &Phe,
    boost_param: BoostParam,
) {
    assert_eq!(losss.inner().len(), predictions.m() * 2);

    let epsilons = epsilon::calculate_epsilons(ps, ys, boost_param.eps());
    //println!("epsilon case, cont: {:.2e},{:.2e}", epsilons.0, epsilons.1);

    // first calculate E:=case_ws_sum and F:=control_ws_sum
    let ef_ = table::calculate_ef_sum(ps, ys);

    losss
        .inner()
        .par_iter_mut()
        .enumerate()
        .for_each(|(li, loss)| {
            *loss = calculate_loss_gt_constada_nosimd_sm(
                &predictions.vals_snv_s_u8(li / 2, li % 2),
                //predictions.predict_snv_s(li / 2, li % 2),
                ps,
                ys,
                ef_,
                epsilons,
                boost_param.eps(),
            );
        })
    //println!("losss {:?}", losss);
    //losss
}

/// calculate loss for FreeModelMissing
///  loss = M + 2*(sqrt(D2 * N2)+sqrt(D1 * N1)+sqrt(D0 *N0))
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn calculate_loss_gt_freemodelmissing_simd_sm(
    predicts: &GenotSnvRef,
    ps: &[f64],
    phe: &Phe,
    epsilons: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    eps: Eps,
) -> f64 {
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
    //println!("0x7f7f7f7f + 1: {:#010x}", 0x7f7f7f7f + 1 as i32);
    //println!("0x7f7f7f7f + 1: {}", 0x7f7f7f7f + 1 as i32);
    //println!("0x7f7f7f80: {}", 0x7f7f7f80 as i32);

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
        //println!("ni {}", ni);

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

        //println!("mis 0: {:#010b}", mis_s0m[4 * ni]);
        //println!("mis 1: {:#010b}", mis_s0m[4 * ni + 1]);
        //println!("mis 2: {:#010b}", mis_s0m[4 * ni + 2]);
        //println!("mis 3: {:#010b}", mis_s0m[4 * ni + 3]);
        //println!("mis b32: {:#034b}", mis_s0_b32);
        //println!("misv {:?}", misv_s0_32);

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

            //println!("take_mask a s0 {:?}", take_mask_a_s0);
            //println!("take_mask b s0 {:?}", take_mask_b_s0);
            //println!("take_mask a s1 {:?}", take_mask_a_s1);
            //println!("take_mask b s1 {:?}", take_mask_b_s1);

            let psv_lo_ptr = ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let psv_hi_ptr = ps[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();

            //println!("ps ind {}", 32 * ni + 8 * bi);
            //println!("ps ind {}", 32 * ni + 8 * bi + 4);
            //println!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //println!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            let psv_lo: __m256d = _mm256_load_pd(psv_lo_ptr as *const _);
            let psv_hi: __m256d = _mm256_load_pd(psv_hi_ptr as *const _);

            //println!("ps lo {:?}", psv_lo);
            //println!("ps hi {:?}", psv_hi);

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

            //println!("ps a s0 lo {:?}", ps_masked_a_s0_lo);
            //println!("ps a s1 lo {:?}", ps_masked_a_s1_lo);
            //println!("ps b s0 lo {:?}", ps_masked_b_s0_lo);
            //println!("ps b s1 lo {:?}", ps_masked_b_s1_lo);

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

            //println!("a s0 hi {:?}", ps_masked_a_s0_hi);

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

    let (table7, _) = table::fixed_seven_sum(
        (d2_sum, n2_sum, d1_sum, n1_sum, d0_sum, n0_sum),
        epsilons,
        eps,
    );

    //if mi % 20000 == 0 {
    //    println!("table7 {:?}", table7);
    //    println!("cont table {:?}", predicts.stat_contingency_table(phe))
    //}
    //debug
    //println!("table7 {:?}", table7);

    let loss_ = calculate_loss_table7(table7);

    loss_
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn calculate_loss_gt_freemodelmissing_simd(
    losss: &mut LossStruct,
    genot: &Genot,
    ps: &[f64],
    phe: &Phe,
    boost_param: BoostParam,
) {
    let epsilons = epsilon::calculate_epsilons(ps, phe, boost_param.eps());
    println!("epsilon case, cont: {:.4e},{:.4e}", epsilons.0, epsilons.1);

    assert_eq!(losss.inner().len(), genot.m());

    losss
        .inner()
        .par_iter_mut()
        .enumerate()
        .for_each(|(mi, loss)| {
            *loss = calculate_loss_gt_freemodelmissing_simd_sm(
                &genot.to_genot_snv(mi),
                ps,
                phe,
                epsilons,
                boost_param.eps(),
            )
        });
    /*     losss
    .inner()
    .par_iter_mut()
    .enumerate()
    .for_each(|(mi, loss)| unsafe {
        //if mi % 100000 == 0 {
        //    println!("loss mi {}", mi);
        //}
        *loss = calculate_loss_gt_freemodelmissing_simd_sm(
            &genot.to_genot_snv(mi),
            ps,
            phe,
            epsilons,
            boost_param.eps(),
        );
    }); */
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
    //println!("epsilon case, cont: {:.2e},{:.2e}", epsilons.0, epsilons.1);

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
        //println!("predicts len {}", predicts.len());
        unsafe {
            let loss_ = calculate_loss_gt_freemodelmissing_simd_sm(predicts, ps, ys, epsilons, eps);
            *loss = loss_
        }
    }

    //losss.par_chunks_mut(2).enumerate().for_each(|(mi, loss)| {
    losss.par_iter_mut().enumerate().for_each(|(mi, loss)| {
        //println!("mi {}", mi);
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
    ys: &Phe,
    epsilons: (f64, f64),
    eps: Eps,
) -> f64 {
    //let pred = pred_s.vals();
    let (table, _) = table::calculate_table7_sum_sm(&pred_s, ps, ys, epsilons, eps);

    calculate_loss_table7(table)
}

// TODO: this can be integrated to freemodelmissing_simd()
pub fn calculate_loss_gt_freemodelmissing_nosimd(
    losss: &mut LossStruct,
    predictions: &Genot,
    //predictions: &GenotBi<Vec<u8>>,
    //predictions: &[B8],
    ps: &[f64],
    ys: &Phe,
    boost_param: BoostParam,
) {
    assert_eq!(losss.inner().len(), predictions.m());

    let epsilons = epsilon::calculate_epsilons(ps, ys, boost_param.eps());
    //println!("epsilon case, cont: {:.2e},{:.2e}", epsilons.0, epsilons.1);

    losss
        .inner()
        .par_iter_mut()
        .enumerate()
        .for_each(|(mi, loss)| {
            *loss = calculate_loss_gt_freemodelmissing_nosimd_sm(
                &predictions.vals_snv(mi),
                //&predictions.to_genot_snv(mi).vals(),
                //&predictions.vals_snv(mi),
                ps,
                ys,
                epsilons,
                boost_param.eps(),
            );
        })
    //println!("losss {:?}", losss);
    //losss
}

// TODO: since I have integrated loss_gt_simd and loss_gt_nosimd into one func. many tests cannot be tested...
// TODO; should integrate setup_test, setup2_test...
#[cfg(test)]
mod tests {
    use super::*;

    use genetics::alloc;

    fn is_eq_f64(v: f64, w: f64, e: f64) -> bool {
        (v - w).abs() < e
    }

    #[test]
    fn test_cast() {
        // 192u8
        let a: u8 = 0b1100_0000;
        //println!("a: {}", a);
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
    fn setup_test2() -> (Genot, Phe, Vec<f64>) {
        let x = vec![0u8, 1, 1, 2];
        let n = x.len();
        let genot = Genot::new(1, n, x);
        let y = vec![false, true, false, true];
        let phe = Phe::new(y);

        let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        for p in [0.1, 0.2, 0.3, 0.4].iter() {
            ps.push(*p);
        }
        for _ in 0..32 {
            ps.push(0.0);
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

        (genot, phe, ps)
    }

    // with eps
    fn setup_test3() -> (Genot, Phe, Vec<f64>) {
        let x = vec![0u8, 0, 0, 0, 1, 1, 1, 2, 2, 3];
        let n = x.len();
        let genot = Genot::new(1, n, x);
        let y = vec![
            true, true, false, false, true, false, false, false, false, true,
        ];
        let phe = Phe::new(y);

        let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        for _ in 0..n {
            ps.push(0.1);
        }
        for _ in 0..32 {
            ps.push(0.0);
        }
        // ps:  [0.1] * 10
        // ys:
        // table7: [0.0,0.2,0.1,0.2,0.2,0.2,0.1]

        (genot, phe, ps)
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx2")]
    #[test]
    fn test_calculate_loss_gt_constada_simd_2() {
        // This error below means SIMD memory is not aligned.
        // "process didn't exit successfully: (signal: 11, SIGSEGV: invalid memory reference)"
        let (predictions, phe, ps) = setup_test2();

        let m = predictions.m();
        let losss = vec![0.0f64; 2 * m];
        let mut losss = LossStruct::ConstAda(losss, m, 2);

        calculate_loss_gt_constada_simd(
            &mut losss,
            &predictions,
            &ps,
            &phe,
            BoostParam::new_type1(),
        );
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
        assert!(is_eq_f64(losss.inner()[0], 1.7548090, 1e-7));
        assert!(is_eq_f64(losss.inner()[1], 1.67261622, 1e-7));
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx2")]
    #[test]
    fn test_calculate_loss_gt_constada_nosimd_2() {
        let (predictions, phe, ps) = setup_test2();

        let m = predictions.m();
        let losss = vec![0.0f64; 2 * m];
        let mut losss = LossStruct::ConstAda(losss, m, 2);
        calculate_loss_gt_constada_simd(
            &mut losss,
            &predictions,
            &ps,
            &phe,
            BoostParam::new_type1(),
        );

        let losss_nosimd = vec![0.0f64; 2 * m];
        let mut losss_nosimd = LossStruct::ConstAda(losss_nosimd, m, 2);
        calculate_loss_gt_constada_nosimd(
            &mut losss_nosimd,
            &predictions,
            &ps,
            &phe,
            BoostParam::new_type1(),
        );

        assert!(is_eq_f64(losss_nosimd.inner()[0], losss.inner()[0], 1e-7));
        assert!(is_eq_f64(losss_nosimd.inner()[1], losss.inner()[1], 1e-7));
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx2")]
    #[test]
    fn test_calculate_loss_gt_freemodelmissing_simd_3() {
        // This error below means SIMD memory is not aligned.
        // "process didn't exit successfully: (signal: 11, SIGSEGV: invalid memory reference)"
        let (predictions, phe, ps) = setup_test3();

        let m = predictions.m();
        let losss = vec![0.0f64; m];
        let mut losss = LossStruct::ModelFree(losss, m);

        calculate_loss_gt_freemodelmissing_simd(
            &mut losss,
            &predictions,
            &ps,
            &phe,
            BoostParam::new_type1(),
        );

        // abcd
        // eps: (0.05, 0.05)
        // table7: [0.0,0.2,0.1,0.2,0.2,0.2,0.1]
        //         ->[0.05,0.25,0.15,0.25,0.25,0.25,0.1]
        // -> not realistic due to small n
        assert!(is_eq_f64(losss.inner()[0], 1.21090513, 1e-7));
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx2")]
    #[test]
    fn test_calculate_loss_gt_freemodelmissing_nosimd_3() {
        let (predictions, phe, ps) = setup_test2();

        let m = predictions.m();
        let losss = vec![0.0f64; m];
        let mut losss = LossStruct::ModelFree(losss, m);
        calculate_loss_gt_freemodelmissing_simd(
            &mut losss,
            &predictions,
            &ps,
            &phe,
            BoostParam::new_type1(),
        );

        let losss_nosimd = vec![0.0f64; m];
        let mut losss_nosimd = LossStruct::ModelFree(losss_nosimd, m);
        calculate_loss_gt_freemodelmissing_nosimd(
            &mut losss_nosimd,
            &predictions,
            &ps,
            &phe,
            BoostParam::new_type1(),
        );

        assert!(is_eq_f64(losss_nosimd.inner()[0], losss.inner()[0], 1e-7));
    }
}
