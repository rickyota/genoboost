use genetics::genot::prelude::*;

use super::Sum3by3;

pub fn calc_ws_sum_logit_mi(
    gsnv: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
) -> ((f64, f64, f64), (f64, f64, f64)) {
    unsafe {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                //let (wzs_sum, wls_sum) = calc::calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
                return calc_ws_sum_logit_simd_mi(gsnv, wzs_pad, wls_pad);
            }
        }
        return calc_ws_sum_logit_nosimd_mi(gsnv, wzs_pad, wls_pad);
    }
}

// unnecessary: implement no missing ver.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn calc_ws_sum_logit_simd_mi(
    gsnv: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
    //phe: &Phe,
    //epsilons_wls: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //epsilons_wzs: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //eps: Eps,
) -> ((f64, f64, f64), (f64, f64, f64)) {
    //) -> Coef {

    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    //use std::convert::TryInto;

    // (x) n + 32
    let n = gsnv.n();

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

    let mut wzs_sum2_acc = _mm256_setzero_pd();
    let mut wzs_sum1_acc = _mm256_setzero_pd();
    let mut wzs_sum0_acc = _mm256_setzero_pd();
    let mut wls_sum2_acc = _mm256_setzero_pd();
    let mut wls_sum1_acc = _mm256_setzero_pd();
    let mut wls_sum0_acc = _mm256_setzero_pd();

    let shifts: [__m256i; 4] = [
        _mm256_set1_epi32(24),
        _mm256_set1_epi32(16),
        _mm256_set1_epi32(8),
        _mm256_set1_epi32(0),
    ];

    for ni in 0..(n / 32 + 1) {
        //log::debug!("ni {}", ni);
        //println!("ni {}", ni);

        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let yv_32 = _mm256_set1_epi32(ys_b32 as i32);

        // bitwise not is fastest using xor by xor(v, ones)
        // sum2=( pred0 & pred1 )
        // sum1= ( !pred1 & pred0 )
        // sum0= ( !pred0 & (!pred1) )
        // D2sum = y & pred0 & pred1 = y & ( pred0 & pred1 )
        // N2sum = !y & pred0 & pred1 = !y & ( pred0 & pred1 )
        // D1sum = y & pred0 & !pred1 = y & ( !pred1 & pred0 )
        // N1sum = !y & pred0 & !pred1 = !y & ( !pred1 & pred0 )
        // D0sum = y & !pred0 & !pred1 = !pred0 & ( !pred1 & y )
        // N0sum = !y & !pred0 & !pred1 = !y & ( !pred0 & (!pred1) )
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

            //log::debug!("bi {}", bi);

            let shift_v = shifts[bi];
            let take_mask_2 = _mm256_sllv_epi32(take_mask_2_32, shift_v);
            let take_mask_1 = _mm256_sllv_epi32(take_mask_1_32, shift_v);
            let take_mask_0 = _mm256_sllv_epi32(take_mask_0_32, shift_v);

            //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
            //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
            //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
            //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);

            let wzsv_lo_ptr = wzs_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let wzsv_hi_ptr = wzs_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            let wlsv_lo_ptr = wls_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let wlsv_hi_ptr = wls_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            //let psv_lo_ptr = ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            //let psv_hi_ptr = ps[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            //println!(
            //    "wzsv_lo_ptr  val {:?}",
            //    32 * ni + 8 * bi..32 * ni + 8 * bi + 4
            //);

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            //println!("wzsv_lo_ptr {:?}", wzsv_lo_ptr);
            let wzsv_lo: __m256d = _mm256_load_pd(wzsv_lo_ptr as *const _);
            //println!("point 5.1");
            let wzsv_hi: __m256d = _mm256_load_pd(wzsv_hi_ptr as *const _);
            //println!("point 5.2");
            let wlsv_lo: __m256d = _mm256_load_pd(wlsv_lo_ptr as *const _);
            //println!("point 5.3");
            let wlsv_hi: __m256d = _mm256_load_pd(wlsv_hi_ptr as *const _);

            //log::debug!("ps lo {:?}", psv_lo);
            //log::debug!("ps hi {:?}", psv_hi);

            // first for low
            let wzs_masked_2_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_2));
            let wzs_masked_1_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_1));
            let wzs_masked_0_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_0));

            let wls_masked_2_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_2));
            let wls_masked_1_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_1));
            let wls_masked_0_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_0));

            wzs_sum2_acc = _mm256_add_pd(wzs_sum2_acc, wzs_masked_2_lo);
            wzs_sum1_acc = _mm256_add_pd(wzs_sum1_acc, wzs_masked_1_lo);
            wzs_sum0_acc = _mm256_add_pd(wzs_sum0_acc, wzs_masked_0_lo);
            wls_sum2_acc = _mm256_add_pd(wls_sum2_acc, wls_masked_2_lo);
            wls_sum1_acc = _mm256_add_pd(wls_sum1_acc, wls_masked_1_lo);
            wls_sum0_acc = _mm256_add_pd(wls_sum0_acc, wls_masked_0_lo);

            //log::debug!("ps a s0 lo {:?}", ps_masked_a_s0_lo);
            //log::debug!("ps a s1 lo {:?}", ps_masked_a_s1_lo);
            //log::debug!("ps b s0 lo {:?}", ps_masked_b_s0_lo);
            //log::debug!("ps b s1 lo {:?}", ps_masked_b_s1_lo);

            // for high
            let take_mask_2_hi = _mm256_slli_epi64(take_mask_2, 32);
            let take_mask_1_hi = _mm256_slli_epi64(take_mask_1, 32);
            let take_mask_0_hi = _mm256_slli_epi64(take_mask_0, 32);

            let wzs_masked_2_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_2_hi));
            let wzs_masked_1_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_1_hi));
            let wzs_masked_0_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_0_hi));
            let wls_masked_2_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_2_hi));
            let wls_masked_1_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_1_hi));
            let wls_masked_0_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_0_hi));

            //log::debug!("a s0 hi {:?}", ps_masked_a_s0_hi);

            wzs_sum2_acc = _mm256_add_pd(wzs_sum2_acc, wzs_masked_2_hi);
            wzs_sum1_acc = _mm256_add_pd(wzs_sum1_acc, wzs_masked_1_hi);
            wzs_sum0_acc = _mm256_add_pd(wzs_sum0_acc, wzs_masked_0_hi);
            wls_sum2_acc = _mm256_add_pd(wls_sum2_acc, wls_masked_2_hi);
            wls_sum1_acc = _mm256_add_pd(wls_sum1_acc, wls_masked_1_hi);
            wls_sum0_acc = _mm256_add_pd(wls_sum0_acc, wls_masked_0_hi);
        }
    }

    // sum 4 double horizontally to get the whole sum
    wzs_sum2_acc = _mm256_hadd_pd(wzs_sum2_acc, wzs_sum2_acc);
    wzs_sum1_acc = _mm256_hadd_pd(wzs_sum1_acc, wzs_sum1_acc);
    wzs_sum0_acc = _mm256_hadd_pd(wzs_sum0_acc, wzs_sum0_acc);
    wls_sum2_acc = _mm256_hadd_pd(wls_sum2_acc, wls_sum2_acc);
    wls_sum1_acc = _mm256_hadd_pd(wls_sum1_acc, wls_sum1_acc);
    wls_sum0_acc = _mm256_hadd_pd(wls_sum0_acc, wls_sum0_acc);

    // 1. any way to hadd??
    // 2. _mm256_extractf128_pd and _mm256_cvtsd_f64: get 64:0

    let wzs_sum2: f64 =
        _mm256_cvtsd_f64(wzs_sum2_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum2_acc, 1));
    let wzs_sum1: f64 =
        _mm256_cvtsd_f64(wzs_sum1_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum1_acc, 1));
    let wzs_sum0: f64 =
        _mm256_cvtsd_f64(wzs_sum0_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum0_acc, 1));
    let wls_sum2: f64 =
        _mm256_cvtsd_f64(wls_sum2_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum2_acc, 1));
    let wls_sum1: f64 =
        _mm256_cvtsd_f64(wls_sum1_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum1_acc, 1));
    let wls_sum0: f64 =
        _mm256_cvtsd_f64(wls_sum0_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum0_acc, 1));

    (
        (wzs_sum2, wzs_sum1, wzs_sum0),
        (wls_sum2, wls_sum1, wls_sum0),
    )

    //let s2=wzs_sum2/wls_sum2;
    //let s1=wzs_sum1/wls_sum1;
    //let s0=wzs_sum0/wls_sum0;
    //let sm=0.0f64;

    //Coef::Score4((s0, s1, s2, sm))
}

/// wzs_sum.0 = \sum_{g=2} U
/// wls_sum.0 = \sum_{g=2} W
pub fn calc_ws_sum_logit_nosimd_mi(
    gsnv: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
) -> ((f64, f64, f64), (f64, f64, f64)) {
    // TODO: reverse order to (0,1,2)
    // g*g 2,1,0
    let mut wzs_sum = (0.0f64, 0.0f64, 0.0f64);
    let mut wls_sum = (0.0f64, 0.0f64, 0.0f64);

    gsnv.iter().enumerate().for_each(|(mi, g)| {
        let wzs = wzs_pad[mi];
        let wls = wls_pad[mi];
        match g {
            0 => {
                wzs_sum.2 += wzs;
                wls_sum.2 += wls;
                //wzs_sum.0 += wzs;
                //wls_sum.0 += wls;
            }
            1 => {
                wzs_sum.1 += wzs;
                wls_sum.1 += wls;
            }
            2 => {
                wzs_sum.0 += wzs;
                wls_sum.0 += wls;
                //wzs_sum.2 += wzs;
                //wls_sum.2 += wls;
            }
            _ => panic!("genot must be 0,1,2."),
        }
    });

    (wzs_sum, wls_sum)
}

//((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)),
//((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)),
pub fn calc_ws_sum_logit_interaction_mi(
    gsnv_1: &GenotSnvRef,
    gsnv_2: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
) -> (Sum3by3, Sum3by3) {
    unsafe {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                //let (wzs_sum, wls_sum) = calc::calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
                return calc_ws_sum_logit_interaction_simd_mi(gsnv_1, gsnv_2, wzs_pad, wls_pad);
            }
        }
        return calc_ws_sum_logit_interaction_nosimd_mi(gsnv_1, gsnv_2, wzs_pad, wls_pad);
    }
}

//pub fn calc_ws_sum_logit_interaction_mi(
//    gsnv_1: &GenotSnvRef,
//    gsnv_2: &GenotSnvRef,
//    wzs_pad: &[f64],
//    wls_pad: &[f64],
//) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
//    unsafe {
//        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
//        {
//            if is_x86_feature_detected!("avx2") {
//                //let (wzs_sum, wls_sum) = calc::calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
//                return calc_ws_sum_logit_interaction_simd_mi(gsnv_1, gsnv_2, wzs_pad, wls_pad);
//            }
//        }
//        return calc_ws_sum_logit_interaction_nosimd_mi(gsnv_1, gsnv_2, wzs_pad, wls_pad);
//    }
//}

pub fn calc_ws_sum_logit_interaction_nosimd_mi(
    gsnv_1: &GenotSnvRef,
    gsnv_2: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
) -> (Sum3by3, Sum3by3) {
    //) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {

    // For g1, g2
    // ((x00, x01, x02),(x10,x11,x12),(x20,x21,x22))
    let init3 = (0.0f64, 0.0f64, 0.0f64);
    let mut wzs_sum = (init3, init3, init3);
    let mut wls_sum = (init3, init3, init3);

    gsnv_1
        .iter()
        .zip(gsnv_2.iter())
        .enumerate()
        .for_each(|(mi, (g1, g2))| {
            let wzs = wzs_pad[mi];
            let wls = wls_pad[mi];
            match (g1, g2) {
                (0, 0) => {
                    wzs_sum.0 .0 += wzs;
                    wls_sum.0 .0 += wls;
                }
                (0, 1) => {
                    wzs_sum.0 .1 += wzs;
                    wls_sum.0 .1 += wls;
                }
                (0, 2) => {
                    wzs_sum.0 .2 += wzs;
                    wls_sum.0 .2 += wls;
                }
                (1, 0) => {
                    wzs_sum.1 .0 += wzs;
                    wls_sum.1 .0 += wls;
                }
                (1, 1) => {
                    wzs_sum.1 .1 += wzs;
                    wls_sum.1 .1 += wls;
                }
                (1, 2) => {
                    wzs_sum.1 .2 += wzs;
                    wls_sum.1 .2 += wls;
                }
                (2, 0) => {
                    wzs_sum.2 .0 += wzs;
                    wls_sum.2 .0 += wls;
                }
                (2, 1) => {
                    wzs_sum.2 .1 += wzs;
                    wls_sum.2 .1 += wls;
                }
                (2, 2) => {
                    wzs_sum.2 .2 += wzs;
                    wls_sum.2 .2 += wls;
                }
                _ => panic!("Wrong genot."),
            }
        });

    (Sum3by3::new(wzs_sum), Sum3by3::new(wls_sum))
}

//pub fn calc_ws_sum_logit_interaction_nosimd_mi(
//    gsnv_1: &GenotSnvRef,
//    gsnv_2: &GenotSnvRef,
//    wzs_pad: &[f64],
//    wls_pad: &[f64],
//) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
//    // TODO: reverse order to (0,1,2,4)
//    // g*g 4,2,1,0
//    let mut wzs_sum = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
//    let mut wls_sum = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
//
//    gsnv_1
//        .iter()
//        .zip(gsnv_2.iter())
//        .enumerate()
//        .for_each(|(mi, (g1, g2))| {
//            let wzs = wzs_pad[mi];
//            let wls = wls_pad[mi];
//            match g1 * g2 {
//                0 => {
//                    wzs_sum.3 += wzs;
//                    wls_sum.3 += wls;
//                    //wzs_sum.0 += wzs;
//                    //wls_sum.0 += wls;
//                }
//                1 => {
//                    wzs_sum.2 += wzs;
//                    wls_sum.2 += wls;
//                    //wzs_sum.1 += wzs;
//                    //wls_sum.1 += wls;
//                }
//                2 => {
//                    wzs_sum.1 += wzs;
//                    wls_sum.1 += wls;
//                    //wzs_sum.2 += wzs;
//                    //wls_sum.2 += wls;
//                }
//                4 => {
//                    wzs_sum.0 += wzs;
//                    wls_sum.0 += wls;
//                    //wzs_sum.3 += wzs;
//                    //wls_sum.3 += wls;
//                }
//                _ => panic!("genot*genot must be 0,1,2,4."),
//            }
//        });
//
//    (wzs_sum, wls_sum)
//}

// TODO: calculate time vs nosimd
// assume no missing in gsnv
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn calc_ws_sum_logit_interaction_simd_mi(
    gsnv_1: &GenotSnvRef,
    gsnv_2: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
) -> (Sum3by3, Sum3by3) {
    //) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    //use std::convert::TryInto;

    // n + 32
    let n = gsnv_1.n();

    let pred1_s0m = gsnv_1.predict_s(0);
    let pred1_s1m = gsnv_1.predict_s(1);
    let pred2_s0m = gsnv_2.predict_s(0);
    let pred2_s1m = gsnv_2.predict_s(1);

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

    let mut wzs_sum_0_0_acc = _mm256_setzero_pd();
    let mut wzs_sum_0_1_acc = _mm256_setzero_pd();
    let mut wzs_sum_0_2_acc = _mm256_setzero_pd();
    let mut wzs_sum_1_0_acc = _mm256_setzero_pd();
    let mut wzs_sum_1_1_acc = _mm256_setzero_pd();
    let mut wzs_sum_1_2_acc = _mm256_setzero_pd();
    let mut wzs_sum_2_0_acc = _mm256_setzero_pd();
    let mut wzs_sum_2_1_acc = _mm256_setzero_pd();
    let mut wzs_sum_2_2_acc = _mm256_setzero_pd();

    let mut wls_sum_0_0_acc = _mm256_setzero_pd();
    let mut wls_sum_0_1_acc = _mm256_setzero_pd();
    let mut wls_sum_0_2_acc = _mm256_setzero_pd();
    let mut wls_sum_1_0_acc = _mm256_setzero_pd();
    let mut wls_sum_1_1_acc = _mm256_setzero_pd();
    let mut wls_sum_1_2_acc = _mm256_setzero_pd();
    let mut wls_sum_2_0_acc = _mm256_setzero_pd();
    let mut wls_sum_2_1_acc = _mm256_setzero_pd();
    let mut wls_sum_2_2_acc = _mm256_setzero_pd();

    //let mut wls_sum4_acc = _mm256_setzero_pd();
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
        //println!("ni {}", ni);

        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        let pred1_s0_b32 = u32::from_le_bytes(pred1_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred1_s1_b32 = u32::from_le_bytes(pred1_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred2_s0_b32 = u32::from_le_bytes(pred2_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred2_s1_b32 = u32::from_le_bytes(pred2_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let predv1_s0_32 = _mm256_set1_epi32(pred1_s0_b32 as i32);
        let predv1_s1_32 = _mm256_set1_epi32(pred1_s1_b32 as i32);
        let predv2_s0_32 = _mm256_set1_epi32(pred2_s0_b32 as i32);
        let predv2_s1_32 = _mm256_set1_epi32(pred2_s1_b32 as i32);

        //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let yv_32 = _mm256_set1_epi32(ys_b32 as i32);

        // bitwise not is fastest using xor by not(v) = xor(v, ones)
        // sum2=( pred0 & pred1 )
        // sum1= ( !pred1 & pred0 )
        // sum0= ( !pred0 & (!pred1) )
        let flagv1_2_32 = _mm256_and_si256(predv1_s0_32, predv1_s1_32);
        let flagv1_1_32 = _mm256_andnot_si256(predv1_s1_32, predv1_s0_32);
        let flagv1_0_32 = _mm256_andnot_si256(predv1_s0_32, _mm256_xor_si256(predv1_s1_32, ones));

        let flagv2_2_32 = _mm256_and_si256(predv2_s0_32, predv2_s1_32);
        let flagv2_1_32 = _mm256_andnot_si256(predv2_s1_32, predv2_s0_32);
        let flagv2_0_32 = _mm256_andnot_si256(predv2_s0_32, _mm256_xor_si256(predv2_s1_32, ones));

        // interaction
        // (g1, g2)=(0,0)
        let flagv_by_0_0_32 = _mm256_and_si256(flagv1_0_32, flagv2_0_32);
        // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
        // abcdefgh(x4)|...
        // -> extracted (at highest position)
        // 000d0000(x4)|...
        // -> mask
        // 11111111|00000000|00000000|11111111|...
        let flagv_by_0_0_32_ext = _mm256_and_si256(flagv_by_0_0_32, bit_ext_mask);
        let take_mask_0_0_32 = _mm256_cmpeq_epi8(flagv_by_0_0_32_ext, bit_ext_mask);

        let flagv_by_0_1_32 = _mm256_and_si256(flagv1_0_32, flagv2_1_32);
        let flagv_by_0_1_32_ext = _mm256_and_si256(flagv_by_0_1_32, bit_ext_mask);
        let take_mask_0_1_32 = _mm256_cmpeq_epi8(flagv_by_0_1_32_ext, bit_ext_mask);
        let flagv_by_0_2_32 = _mm256_and_si256(flagv1_0_32, flagv2_2_32);
        let flagv_by_0_2_32_ext = _mm256_and_si256(flagv_by_0_2_32, bit_ext_mask);
        let take_mask_0_2_32 = _mm256_cmpeq_epi8(flagv_by_0_2_32_ext, bit_ext_mask);
        let flagv_by_1_0_32 = _mm256_and_si256(flagv1_1_32, flagv2_0_32);
        let flagv_by_1_0_32_ext = _mm256_and_si256(flagv_by_1_0_32, bit_ext_mask);
        let take_mask_1_0_32 = _mm256_cmpeq_epi8(flagv_by_1_0_32_ext, bit_ext_mask);
        let flagv_by_1_1_32 = _mm256_and_si256(flagv1_1_32, flagv2_1_32);
        let flagv_by_1_1_32_ext = _mm256_and_si256(flagv_by_1_1_32, bit_ext_mask);
        let take_mask_1_1_32 = _mm256_cmpeq_epi8(flagv_by_1_1_32_ext, bit_ext_mask);
        let flagv_by_1_2_32 = _mm256_and_si256(flagv1_1_32, flagv2_2_32);
        let flagv_by_1_2_32_ext = _mm256_and_si256(flagv_by_1_2_32, bit_ext_mask);
        let take_mask_1_2_32 = _mm256_cmpeq_epi8(flagv_by_1_2_32_ext, bit_ext_mask);
        let flagv_by_2_0_32 = _mm256_and_si256(flagv1_2_32, flagv2_0_32);
        let flagv_by_2_0_32_ext = _mm256_and_si256(flagv_by_2_0_32, bit_ext_mask);
        let take_mask_2_0_32 = _mm256_cmpeq_epi8(flagv_by_2_0_32_ext, bit_ext_mask);
        let flagv_by_2_1_32 = _mm256_and_si256(flagv1_2_32, flagv2_1_32);
        let flagv_by_2_1_32_ext = _mm256_and_si256(flagv_by_2_1_32, bit_ext_mask);
        let take_mask_2_1_32 = _mm256_cmpeq_epi8(flagv_by_2_1_32_ext, bit_ext_mask);
        let flagv_by_2_2_32 = _mm256_and_si256(flagv1_2_32, flagv2_2_32);
        let flagv_by_2_2_32_ext = _mm256_and_si256(flagv_by_2_2_32, bit_ext_mask);
        let take_mask_2_2_32 = _mm256_cmpeq_epi8(flagv_by_2_2_32_ext, bit_ext_mask);

        // bi=0-3, shift=24,16,8,0
        //const SHIFTS: &'static [i32] = &[24, 16, 8, 0];

        for bi in 0usize..4 {
            // DCBADCBA...DCBA
            // -> b=1: for B
            // BA00BA00...BA00

            //log::debug!("bi {}", bi);
            //println!("bi {}", bi);

            let shift_v = shifts[bi];
            let take_mask_0_0 = _mm256_sllv_epi32(take_mask_0_0_32, shift_v);
            let take_mask_0_1 = _mm256_sllv_epi32(take_mask_0_1_32, shift_v);
            let take_mask_0_2 = _mm256_sllv_epi32(take_mask_0_2_32, shift_v);
            let take_mask_1_0 = _mm256_sllv_epi32(take_mask_1_0_32, shift_v);
            let take_mask_1_1 = _mm256_sllv_epi32(take_mask_1_1_32, shift_v);
            let take_mask_1_2 = _mm256_sllv_epi32(take_mask_1_2_32, shift_v);
            let take_mask_2_0 = _mm256_sllv_epi32(take_mask_2_0_32, shift_v);
            let take_mask_2_1 = _mm256_sllv_epi32(take_mask_2_1_32, shift_v);
            let take_mask_2_2 = _mm256_sllv_epi32(take_mask_2_2_32, shift_v);

            let wzsv_lo_ptr = wzs_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let wzsv_hi_ptr = wzs_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            let wlsv_lo_ptr = wls_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let wlsv_hi_ptr = wls_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();

            //println!(
            //    "wzsv_lo_ptr  val {:?}",
            //    32 * ni + 8 * bi..32 * ni + 8 * bi + 4
            //);

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            //println!("wzsv_lo_ptr {:?}", wzsv_lo_ptr);
            let wzsv_lo: __m256d = _mm256_load_pd(wzsv_lo_ptr as *const _);
            let wzsv_hi: __m256d = _mm256_load_pd(wzsv_hi_ptr as *const _);
            let wlsv_lo: __m256d = _mm256_load_pd(wlsv_lo_ptr as *const _);
            let wlsv_hi: __m256d = _mm256_load_pd(wlsv_hi_ptr as *const _);

            //log::debug!("ps lo {:?}", psv_lo);
            //log::debug!("ps hi {:?}", psv_hi);

            // first for low
            let wzs_masked_0_0_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_0_0));
            wzs_sum_0_0_acc = _mm256_add_pd(wzs_sum_0_0_acc, wzs_masked_0_0_lo);
            let wzs_masked_0_1_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_0_1));
            wzs_sum_0_1_acc = _mm256_add_pd(wzs_sum_0_1_acc, wzs_masked_0_1_lo);
            let wzs_masked_0_2_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_0_2));
            wzs_sum_0_2_acc = _mm256_add_pd(wzs_sum_0_2_acc, wzs_masked_0_2_lo);
            let wzs_masked_1_0_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_1_0));
            wzs_sum_1_0_acc = _mm256_add_pd(wzs_sum_1_0_acc, wzs_masked_1_0_lo);
            let wzs_masked_1_1_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_1_1));
            wzs_sum_1_1_acc = _mm256_add_pd(wzs_sum_1_1_acc, wzs_masked_1_1_lo);
            let wzs_masked_1_2_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_1_2));
            wzs_sum_1_2_acc = _mm256_add_pd(wzs_sum_1_2_acc, wzs_masked_1_2_lo);
            let wzs_masked_2_0_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_2_0));
            wzs_sum_2_0_acc = _mm256_add_pd(wzs_sum_2_0_acc, wzs_masked_2_0_lo);
            let wzs_masked_2_1_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_2_1));
            wzs_sum_2_1_acc = _mm256_add_pd(wzs_sum_2_1_acc, wzs_masked_2_1_lo);
            let wzs_masked_2_2_lo =
                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_2_2));
            wzs_sum_2_2_acc = _mm256_add_pd(wzs_sum_2_2_acc, wzs_masked_2_2_lo);

            let wls_masked_0_0_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_0_0));
            wls_sum_0_0_acc = _mm256_add_pd(wls_sum_0_0_acc, wls_masked_0_0_lo);
            let wls_masked_0_1_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_0_1));
            wls_sum_0_1_acc = _mm256_add_pd(wls_sum_0_1_acc, wls_masked_0_1_lo);
            let wls_masked_0_2_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_0_2));
            wls_sum_0_2_acc = _mm256_add_pd(wls_sum_0_2_acc, wls_masked_0_2_lo);
            let wls_masked_1_0_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_1_0));
            wls_sum_1_0_acc = _mm256_add_pd(wls_sum_1_0_acc, wls_masked_1_0_lo);
            let wls_masked_1_1_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_1_1));
            wls_sum_1_1_acc = _mm256_add_pd(wls_sum_1_1_acc, wls_masked_1_1_lo);
            let wls_masked_1_2_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_1_2));
            wls_sum_1_2_acc = _mm256_add_pd(wls_sum_1_2_acc, wls_masked_1_2_lo);
            let wls_masked_2_0_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_2_0));
            wls_sum_2_0_acc = _mm256_add_pd(wls_sum_2_0_acc, wls_masked_2_0_lo);
            let wls_masked_2_1_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_2_1));
            wls_sum_2_1_acc = _mm256_add_pd(wls_sum_2_1_acc, wls_masked_2_1_lo);
            let wls_masked_2_2_lo =
                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_2_2));
            wls_sum_2_2_acc = _mm256_add_pd(wls_sum_2_2_acc, wls_masked_2_2_lo);

            // for high
            let take_mask_0_0_hi = _mm256_slli_epi64(take_mask_0_0, 32);
            let take_mask_0_1_hi = _mm256_slli_epi64(take_mask_0_1, 32);
            let take_mask_0_2_hi = _mm256_slli_epi64(take_mask_0_2, 32);
            let take_mask_1_0_hi = _mm256_slli_epi64(take_mask_1_0, 32);
            let take_mask_1_1_hi = _mm256_slli_epi64(take_mask_1_1, 32);
            let take_mask_1_2_hi = _mm256_slli_epi64(take_mask_1_2, 32);
            let take_mask_2_0_hi = _mm256_slli_epi64(take_mask_2_0, 32);
            let take_mask_2_1_hi = _mm256_slli_epi64(take_mask_2_1, 32);
            let take_mask_2_2_hi = _mm256_slli_epi64(take_mask_2_2, 32);

            //let take_mask_4_hi = _mm256_slli_epi64(take_mask_4, 32);
            //let take_mask_2_hi = _mm256_slli_epi64(take_mask_2, 32);
            //let take_mask_1_hi = _mm256_slli_epi64(take_mask_1, 32);
            //let take_mask_0_hi = _mm256_slli_epi64(take_mask_0, 32);

            let wzs_masked_0_0_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_0_0_hi));
            wzs_sum_0_0_acc = _mm256_add_pd(wzs_sum_0_0_acc, wzs_masked_0_0_hi);
            let wzs_masked_0_1_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_0_1_hi));
            wzs_sum_0_1_acc = _mm256_add_pd(wzs_sum_0_1_acc, wzs_masked_0_1_hi);
            let wzs_masked_0_2_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_0_2_hi));
            wzs_sum_0_2_acc = _mm256_add_pd(wzs_sum_0_2_acc, wzs_masked_0_2_hi);
            let wzs_masked_1_0_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_1_0_hi));
            wzs_sum_1_0_acc = _mm256_add_pd(wzs_sum_1_0_acc, wzs_masked_1_0_hi);
            let wzs_masked_1_1_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_1_1_hi));
            wzs_sum_1_1_acc = _mm256_add_pd(wzs_sum_1_1_acc, wzs_masked_1_1_hi);
            let wzs_masked_1_2_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_1_2_hi));
            wzs_sum_1_2_acc = _mm256_add_pd(wzs_sum_1_2_acc, wzs_masked_1_2_hi);
            let wzs_masked_2_0_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_2_0_hi));
            wzs_sum_2_0_acc = _mm256_add_pd(wzs_sum_2_0_acc, wzs_masked_2_0_hi);
            let wzs_masked_2_1_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_2_1_hi));
            wzs_sum_2_1_acc = _mm256_add_pd(wzs_sum_2_1_acc, wzs_masked_2_1_hi);
            let wzs_masked_2_2_hi =
                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_2_2_hi));
            wzs_sum_2_2_acc = _mm256_add_pd(wzs_sum_2_2_acc, wzs_masked_2_2_hi);

            let wls_masked_0_0_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_0_0_hi));
            wls_sum_0_0_acc = _mm256_add_pd(wls_sum_0_0_acc, wls_masked_0_0_hi);
            let wls_masked_0_1_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_0_1_hi));
            wls_sum_0_1_acc = _mm256_add_pd(wls_sum_0_1_acc, wls_masked_0_1_hi);
            let wls_masked_0_2_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_0_2_hi));
            wls_sum_0_2_acc = _mm256_add_pd(wls_sum_0_2_acc, wls_masked_0_2_hi);
            let wls_masked_1_0_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_1_0_hi));
            wls_sum_1_0_acc = _mm256_add_pd(wls_sum_1_0_acc, wls_masked_1_0_hi);
            let wls_masked_1_1_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_1_1_hi));
            wls_sum_1_1_acc = _mm256_add_pd(wls_sum_1_1_acc, wls_masked_1_1_hi);
            let wls_masked_1_2_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_1_2_hi));
            wls_sum_1_2_acc = _mm256_add_pd(wls_sum_1_2_acc, wls_masked_1_2_hi);
            let wls_masked_2_0_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_2_0_hi));
            wls_sum_2_0_acc = _mm256_add_pd(wls_sum_2_0_acc, wls_masked_2_0_hi);
            let wls_masked_2_1_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_2_1_hi));
            wls_sum_2_1_acc = _mm256_add_pd(wls_sum_2_1_acc, wls_masked_2_1_hi);
            let wls_masked_2_2_hi =
                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_2_2_hi));
            wls_sum_2_2_acc = _mm256_add_pd(wls_sum_2_2_acc, wls_masked_2_2_hi);
        }
    }

    // sum 4 double horizontally to get the whole sum
    wzs_sum_0_0_acc = _mm256_hadd_pd(wzs_sum_0_0_acc, wzs_sum_0_0_acc);
    wzs_sum_0_1_acc = _mm256_hadd_pd(wzs_sum_0_1_acc, wzs_sum_0_1_acc);
    wzs_sum_0_2_acc = _mm256_hadd_pd(wzs_sum_0_2_acc, wzs_sum_0_2_acc);
    wzs_sum_1_0_acc = _mm256_hadd_pd(wzs_sum_1_0_acc, wzs_sum_1_0_acc);
    wzs_sum_1_1_acc = _mm256_hadd_pd(wzs_sum_1_1_acc, wzs_sum_1_1_acc);
    wzs_sum_1_2_acc = _mm256_hadd_pd(wzs_sum_1_2_acc, wzs_sum_1_2_acc);
    wzs_sum_2_0_acc = _mm256_hadd_pd(wzs_sum_2_0_acc, wzs_sum_2_0_acc);
    wzs_sum_2_1_acc = _mm256_hadd_pd(wzs_sum_2_1_acc, wzs_sum_2_1_acc);
    wzs_sum_2_2_acc = _mm256_hadd_pd(wzs_sum_2_2_acc, wzs_sum_2_2_acc);

    wls_sum_0_0_acc = _mm256_hadd_pd(wls_sum_0_0_acc, wls_sum_0_0_acc);
    wls_sum_0_1_acc = _mm256_hadd_pd(wls_sum_0_1_acc, wls_sum_0_1_acc);
    wls_sum_0_2_acc = _mm256_hadd_pd(wls_sum_0_2_acc, wls_sum_0_2_acc);
    wls_sum_1_0_acc = _mm256_hadd_pd(wls_sum_1_0_acc, wls_sum_1_0_acc);
    wls_sum_1_1_acc = _mm256_hadd_pd(wls_sum_1_1_acc, wls_sum_1_1_acc);
    wls_sum_1_2_acc = _mm256_hadd_pd(wls_sum_1_2_acc, wls_sum_1_2_acc);
    wls_sum_2_0_acc = _mm256_hadd_pd(wls_sum_2_0_acc, wls_sum_2_0_acc);
    wls_sum_2_1_acc = _mm256_hadd_pd(wls_sum_2_1_acc, wls_sum_2_1_acc);
    wls_sum_2_2_acc = _mm256_hadd_pd(wls_sum_2_2_acc, wls_sum_2_2_acc);

    //wzs_sum4_acc = _mm256_hadd_pd(wzs_sum4_acc, wzs_sum4_acc);
    //wzs_sum2_acc = _mm256_hadd_pd(wzs_sum2_acc, wzs_sum2_acc);
    //wzs_sum1_acc = _mm256_hadd_pd(wzs_sum1_acc, wzs_sum1_acc);
    //wzs_sum0_acc = _mm256_hadd_pd(wzs_sum0_acc, wzs_sum0_acc);
    //wls_sum4_acc = _mm256_hadd_pd(wls_sum4_acc, wls_sum4_acc);
    //wls_sum2_acc = _mm256_hadd_pd(wls_sum2_acc, wls_sum2_acc);
    //wls_sum1_acc = _mm256_hadd_pd(wls_sum1_acc, wls_sum1_acc);
    //wls_sum0_acc = _mm256_hadd_pd(wls_sum0_acc, wls_sum0_acc);

    // 1. any way to hadd??
    // 2. _mm256_extractf128_pd and _mm256_cvtsd_f64: get 64:0

    let wzs_sum_0_0: f64 = _mm256_cvtsd_f64(wzs_sum_0_0_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_0_0_acc, 1));
    let wzs_sum_0_1: f64 = _mm256_cvtsd_f64(wzs_sum_0_1_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_0_1_acc, 1));
    let wzs_sum_0_2: f64 = _mm256_cvtsd_f64(wzs_sum_0_2_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_0_2_acc, 1));
    let wzs_sum_1_0: f64 = _mm256_cvtsd_f64(wzs_sum_1_0_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_1_0_acc, 1));
    let wzs_sum_1_1: f64 = _mm256_cvtsd_f64(wzs_sum_1_1_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_1_1_acc, 1));
    let wzs_sum_1_2: f64 = _mm256_cvtsd_f64(wzs_sum_1_2_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_1_2_acc, 1));
    let wzs_sum_2_0: f64 = _mm256_cvtsd_f64(wzs_sum_2_0_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_2_0_acc, 1));
    let wzs_sum_2_1: f64 = _mm256_cvtsd_f64(wzs_sum_2_1_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_2_1_acc, 1));
    let wzs_sum_2_2: f64 = _mm256_cvtsd_f64(wzs_sum_2_2_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum_2_2_acc, 1));

    let wls_sum_0_0: f64 = _mm256_cvtsd_f64(wls_sum_0_0_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_0_0_acc, 1));
    let wls_sum_0_1: f64 = _mm256_cvtsd_f64(wls_sum_0_1_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_0_1_acc, 1));
    let wls_sum_0_2: f64 = _mm256_cvtsd_f64(wls_sum_0_2_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_0_2_acc, 1));
    let wls_sum_1_0: f64 = _mm256_cvtsd_f64(wls_sum_1_0_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_1_0_acc, 1));
    let wls_sum_1_1: f64 = _mm256_cvtsd_f64(wls_sum_1_1_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_1_1_acc, 1));
    let wls_sum_1_2: f64 = _mm256_cvtsd_f64(wls_sum_1_2_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_1_2_acc, 1));
    let wls_sum_2_0: f64 = _mm256_cvtsd_f64(wls_sum_2_0_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_2_0_acc, 1));
    let wls_sum_2_1: f64 = _mm256_cvtsd_f64(wls_sum_2_1_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_2_1_acc, 1));
    let wls_sum_2_2: f64 = _mm256_cvtsd_f64(wls_sum_2_2_acc)
        + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum_2_2_acc, 1));

    //let wzs_sum4: f64 =
    //    _mm256_cvtsd_f64(wzs_sum4_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum4_acc, 1));
    //let wzs_sum2: f64 =
    //    _mm256_cvtsd_f64(wzs_sum2_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum2_acc, 1));
    //let wzs_sum1: f64 =
    //    _mm256_cvtsd_f64(wzs_sum1_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum1_acc, 1));
    //let wzs_sum0: f64 =
    //    _mm256_cvtsd_f64(wzs_sum0_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum0_acc, 1));
    //let wls_sum4: f64 =
    //    _mm256_cvtsd_f64(wls_sum4_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum4_acc, 1));
    //let wls_sum2: f64 =
    //    _mm256_cvtsd_f64(wls_sum2_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum2_acc, 1));
    //let wls_sum1: f64 =
    //    _mm256_cvtsd_f64(wls_sum1_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum1_acc, 1));
    //let wls_sum0: f64 =
    //    _mm256_cvtsd_f64(wls_sum0_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum0_acc, 1));

    //(
    //    (wzs_sum4, wzs_sum2, wzs_sum1, wzs_sum0),
    //    (wls_sum4, wls_sum2, wls_sum1, wls_sum0),
    //)

    (
        Sum3by3::new((
            (wzs_sum_0_0, wzs_sum_0_1, wzs_sum_0_2),
            (wzs_sum_1_0, wzs_sum_1_1, wzs_sum_1_2),
            (wzs_sum_2_0, wzs_sum_2_1, wzs_sum_2_2),
        )),
        Sum3by3::new((
            (wls_sum_0_0, wls_sum_0_1, wls_sum_0_2),
            (wls_sum_1_0, wls_sum_1_1, wls_sum_1_2),
            (wls_sum_2_0, wls_sum_2_1, wls_sum_2_2),
        )),
    )
}

//// assume no missing in gsnv
//#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
//#[target_feature(enable = "avx2")]
//pub unsafe fn calc_ws_sum_logit_interaction_simd_mi(
//    gsnv_1: &GenotSnvRef,
//    gsnv_2: &GenotSnvRef,
//    wzs_pad: &[f64],
//    wls_pad: &[f64],
//) -> ((f64, f64, f64, f64), (f64, f64, f64, f64)) {
//    #[cfg(target_arch = "x86")]
//    use std::arch::x86::*;
//    #[cfg(target_arch = "x86_64")]
//    use std::arch::x86_64::*;
//    //use std::convert::TryInto;
//
//    // n + 32
//    let n = gsnv_1.n();
//
//    let pred1_s0m = gsnv_1.predict_s(0);
//    let pred1_s1m = gsnv_1.predict_s(1);
//    let pred2_s0m = gsnv_2.predict_s(0);
//    let pred2_s1m = gsnv_2.predict_s(1);
//
//    let bit_ext_mask: __m256i = _mm256_set_epi32(
//        0x10101010,
//        0x01010101,
//        0x20202020,
//        0x02020202,
//        0x40404040,
//        0x04040404,
//        -0x7f7f7f80, // =0x80808080
//        0x08080808,
//    );
//    let zerod: __m256d = _mm256_setzero_pd();
//    let zeros: __m256i = _mm256_setzero_si256();
//    let ones: __m256i = _mm256_cmpeq_epi32(zeros, zeros);
//
//    let mut wzs_sum4_acc = _mm256_setzero_pd();
//    let mut wzs_sum2_acc = _mm256_setzero_pd();
//    let mut wzs_sum1_acc = _mm256_setzero_pd();
//    let mut wzs_sum0_acc = _mm256_setzero_pd();
//    let mut wls_sum4_acc = _mm256_setzero_pd();
//    let mut wls_sum2_acc = _mm256_setzero_pd();
//    let mut wls_sum1_acc = _mm256_setzero_pd();
//    let mut wls_sum0_acc = _mm256_setzero_pd();
//
//    let shifts: [__m256i; 4] = [
//        _mm256_set1_epi32(24),
//        _mm256_set1_epi32(16),
//        _mm256_set1_epi32(8),
//        _mm256_set1_epi32(0),
//    ];
//
//    for ni in 0..(n / 32 + 1) {
//        //log::debug!("ni {}", ni);
//        //println!("ni {}", ni);
//
//        // broadcast 32bit int to 256bit
//        // ex. DCBA -> DCBADCBA...DCBA
//        // (D=abcdefgh)
//
//        let pred1_s0_b32 = u32::from_le_bytes(pred1_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
//        let pred1_s1_b32 = u32::from_le_bytes(pred1_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
//        let pred2_s0_b32 = u32::from_le_bytes(pred2_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
//        let pred2_s1_b32 = u32::from_le_bytes(pred2_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
//
//        let predv1_s0_32 = _mm256_set1_epi32(pred1_s0_b32 as i32);
//        let predv1_s1_32 = _mm256_set1_epi32(pred1_s1_b32 as i32);
//        let predv2_s0_32 = _mm256_set1_epi32(pred2_s0_b32 as i32);
//        let predv2_s1_32 = _mm256_set1_epi32(pred2_s1_b32 as i32);
//
//        //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
//        //let yv_32 = _mm256_set1_epi32(ys_b32 as i32);
//
//        // bitwise not is fastest using xor by not(v) = xor(v, ones)
//        // sum2=( pred0 & pred1 )
//        // sum1= ( !pred1 & pred0 )
//        // sum0= ( !pred0 & (!pred1) )
//        let flagv1_2_32 = _mm256_and_si256(predv1_s0_32, predv1_s1_32);
//        let flagv1_1_32 = _mm256_andnot_si256(predv1_s1_32, predv1_s0_32);
//
//        let flagv2_2_32 = _mm256_and_si256(predv2_s0_32, predv2_s1_32);
//        let flagv2_1_32 = _mm256_andnot_si256(predv2_s1_32, predv2_s0_32);
//
//        // g1 * g2 = 1,2,4 // remaining is 0
//        let flagv_by_4_32 = _mm256_and_si256(flagv1_2_32, flagv2_2_32);
//        let flagv_by_1_32 = _mm256_and_si256(flagv1_1_32, flagv2_1_32);
//        let flagv_by_2_32_1 = _mm256_and_si256(flagv1_1_32, flagv2_2_32);
//        let flagv_by_2_32_2 = _mm256_and_si256(flagv1_2_32, flagv2_1_32);
//        let flagv_by_2_32 = _mm256_or_si256(flagv_by_2_32_1, flagv_by_2_32_2);
//        // 1 or 2 or 4
//        let flagv_by_124_32 =
//            _mm256_or_si256(flagv_by_1_32, _mm256_or_si256(flagv_by_2_32, flagv_by_4_32));
//        // 0 = not(1 or 2 or 4)
//        let flagv_by_0_32 = _mm256_xor_si256(flagv_by_124_32, ones);
//
//        // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
//        // abcdefgh(x4)|...
//        // -> extracted (at highest position)
//        // 000d0000(x4)|...
//        // -> mask
//        // 11111111|00000000|00000000|11111111|...
//        let flagv_by_4_32_ext = _mm256_and_si256(flagv_by_4_32, bit_ext_mask);
//        let flagv_by_2_32_ext = _mm256_and_si256(flagv_by_2_32, bit_ext_mask);
//        let flagv_by_1_32_ext = _mm256_and_si256(flagv_by_1_32, bit_ext_mask);
//        let flagv_by_0_32_ext = _mm256_and_si256(flagv_by_0_32, bit_ext_mask);
//        //let flagv_0_32_ext = _mm256_and_si256(flagv_0_32, bit_ext_mask);
//
//        let take_mask_4_32 = _mm256_cmpeq_epi8(flagv_by_4_32_ext, bit_ext_mask);
//        let take_mask_2_32 = _mm256_cmpeq_epi8(flagv_by_2_32_ext, bit_ext_mask);
//        let take_mask_1_32 = _mm256_cmpeq_epi8(flagv_by_1_32_ext, bit_ext_mask);
//        let take_mask_0_32 = _mm256_cmpeq_epi8(flagv_by_0_32_ext, bit_ext_mask);
//        //let take_mask_0_32 = _mm256_cmpeq_epi8(flagv_0_32_ext, bit_ext_mask);
//
//        // bi=0-3, shift=24,16,8,0
//        //const SHIFTS: &'static [i32] = &[24, 16, 8, 0];
//
//        for bi in 0usize..4 {
//            // DCBADCBA...DCBA
//            // -> b=1: for B
//            // BA00BA00...BA00
//
//            //log::debug!("bi {}", bi);
//            //println!("bi {}", bi);
//
//            let shift_v = shifts[bi];
//            let take_mask_4 = _mm256_sllv_epi32(take_mask_4_32, shift_v);
//            let take_mask_2 = _mm256_sllv_epi32(take_mask_2_32, shift_v);
//            let take_mask_1 = _mm256_sllv_epi32(take_mask_1_32, shift_v);
//            let take_mask_0 = _mm256_sllv_epi32(take_mask_0_32, shift_v);
//
//            //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
//            //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
//            //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
//            //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);
//
//            let wzsv_lo_ptr = wzs_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
//            let wzsv_hi_ptr = wzs_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
//            let wlsv_lo_ptr = wls_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
//            let wlsv_hi_ptr = wls_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
//            //let psv_lo_ptr = ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
//            //let psv_hi_ptr = ps[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
//            //println!(
//            //    "wzsv_lo_ptr  val {:?}",
//            //    32 * ni + 8 * bi..32 * ni + 8 * bi + 4
//            //);
//
//            //println!("point 5");
//
//            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
//            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
//            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
//            //log::debug!(
//            //    "ps hi {:?}",
//            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
//            //);
//
//            //println!("wzsv_lo_ptr {:?}", wzsv_lo_ptr);
//            let wzsv_lo: __m256d = _mm256_load_pd(wzsv_lo_ptr as *const _);
//            let wzsv_hi: __m256d = _mm256_load_pd(wzsv_hi_ptr as *const _);
//            let wlsv_lo: __m256d = _mm256_load_pd(wlsv_lo_ptr as *const _);
//            let wlsv_hi: __m256d = _mm256_load_pd(wlsv_hi_ptr as *const _);
//
//            //log::debug!("ps lo {:?}", psv_lo);
//            //log::debug!("ps hi {:?}", psv_hi);
//
//            // first for low
//            let wzs_masked_4_lo =
//                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_4));
//            let wzs_masked_2_lo =
//                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_2));
//            let wzs_masked_1_lo =
//                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_1));
//            let wzs_masked_0_lo =
//                _mm256_blendv_pd(zerod, wzsv_lo, _mm256_castsi256_pd(take_mask_0));
//
//            let wls_masked_4_lo =
//                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_4));
//            let wls_masked_2_lo =
//                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_2));
//            let wls_masked_1_lo =
//                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_1));
//            let wls_masked_0_lo =
//                _mm256_blendv_pd(zerod, wlsv_lo, _mm256_castsi256_pd(take_mask_0));
//
//            wzs_sum4_acc = _mm256_add_pd(wzs_sum4_acc, wzs_masked_4_lo);
//            wzs_sum2_acc = _mm256_add_pd(wzs_sum2_acc, wzs_masked_2_lo);
//            wzs_sum1_acc = _mm256_add_pd(wzs_sum1_acc, wzs_masked_1_lo);
//            wzs_sum0_acc = _mm256_add_pd(wzs_sum0_acc, wzs_masked_0_lo);
//            wls_sum4_acc = _mm256_add_pd(wls_sum4_acc, wls_masked_4_lo);
//            wls_sum2_acc = _mm256_add_pd(wls_sum2_acc, wls_masked_2_lo);
//            wls_sum1_acc = _mm256_add_pd(wls_sum1_acc, wls_masked_1_lo);
//            wls_sum0_acc = _mm256_add_pd(wls_sum0_acc, wls_masked_0_lo);
//
//            //log::debug!("ps a s0 lo {:?}", ps_masked_a_s0_lo);
//            //log::debug!("ps a s1 lo {:?}", ps_masked_a_s1_lo);
//            //log::debug!("ps b s0 lo {:?}", ps_masked_b_s0_lo);
//            //log::debug!("ps b s1 lo {:?}", ps_masked_b_s1_lo);
//
//            // for high
//            let take_mask_4_hi = _mm256_slli_epi64(take_mask_4, 32);
//            let take_mask_2_hi = _mm256_slli_epi64(take_mask_2, 32);
//            let take_mask_1_hi = _mm256_slli_epi64(take_mask_1, 32);
//            let take_mask_0_hi = _mm256_slli_epi64(take_mask_0, 32);
//
//            let wzs_masked_4_hi =
//                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_4_hi));
//            let wzs_masked_2_hi =
//                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_2_hi));
//            let wzs_masked_1_hi =
//                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_1_hi));
//            let wzs_masked_0_hi =
//                _mm256_blendv_pd(zerod, wzsv_hi, _mm256_castsi256_pd(take_mask_0_hi));
//            let wls_masked_4_hi =
//                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_4_hi));
//            let wls_masked_2_hi =
//                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_2_hi));
//            let wls_masked_1_hi =
//                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_1_hi));
//            let wls_masked_0_hi =
//                _mm256_blendv_pd(zerod, wlsv_hi, _mm256_castsi256_pd(take_mask_0_hi));
//
//            //log::debug!("a s0 hi {:?}", ps_masked_a_s0_hi);
//
//            wzs_sum4_acc = _mm256_add_pd(wzs_sum4_acc, wzs_masked_4_hi);
//            wzs_sum2_acc = _mm256_add_pd(wzs_sum2_acc, wzs_masked_2_hi);
//            wzs_sum1_acc = _mm256_add_pd(wzs_sum1_acc, wzs_masked_1_hi);
//            wzs_sum0_acc = _mm256_add_pd(wzs_sum0_acc, wzs_masked_0_hi);
//            wls_sum4_acc = _mm256_add_pd(wls_sum4_acc, wls_masked_4_hi);
//            wls_sum2_acc = _mm256_add_pd(wls_sum2_acc, wls_masked_2_hi);
//            wls_sum1_acc = _mm256_add_pd(wls_sum1_acc, wls_masked_1_hi);
//            wls_sum0_acc = _mm256_add_pd(wls_sum0_acc, wls_masked_0_hi);
//        }
//    }
//
//    // sum 4 double horizontally to get the whole sum
//    wzs_sum4_acc = _mm256_hadd_pd(wzs_sum4_acc, wzs_sum4_acc);
//    wzs_sum2_acc = _mm256_hadd_pd(wzs_sum2_acc, wzs_sum2_acc);
//    wzs_sum1_acc = _mm256_hadd_pd(wzs_sum1_acc, wzs_sum1_acc);
//    wzs_sum0_acc = _mm256_hadd_pd(wzs_sum0_acc, wzs_sum0_acc);
//    wls_sum4_acc = _mm256_hadd_pd(wls_sum4_acc, wls_sum4_acc);
//    wls_sum2_acc = _mm256_hadd_pd(wls_sum2_acc, wls_sum2_acc);
//    wls_sum1_acc = _mm256_hadd_pd(wls_sum1_acc, wls_sum1_acc);
//    wls_sum0_acc = _mm256_hadd_pd(wls_sum0_acc, wls_sum0_acc);
//
//    // 1. any way to hadd??
//    // 2. _mm256_extractf128_pd and _mm256_cvtsd_f64: get 64:0
//
//    let wzs_sum4: f64 =
//        _mm256_cvtsd_f64(wzs_sum4_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum4_acc, 1));
//    let wzs_sum2: f64 =
//        _mm256_cvtsd_f64(wzs_sum2_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum2_acc, 1));
//    let wzs_sum1: f64 =
//        _mm256_cvtsd_f64(wzs_sum1_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum1_acc, 1));
//    let wzs_sum0: f64 =
//        _mm256_cvtsd_f64(wzs_sum0_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wzs_sum0_acc, 1));
//    let wls_sum4: f64 =
//        _mm256_cvtsd_f64(wls_sum4_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum4_acc, 1));
//    let wls_sum2: f64 =
//        _mm256_cvtsd_f64(wls_sum2_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum2_acc, 1));
//    let wls_sum1: f64 =
//        _mm256_cvtsd_f64(wls_sum1_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum1_acc, 1));
//    let wls_sum0: f64 =
//        _mm256_cvtsd_f64(wls_sum0_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(wls_sum0_acc, 1));
//
//    (
//        (wzs_sum4, wzs_sum2, wzs_sum1, wzs_sum0),
//        (wls_sum4, wls_sum2, wls_sum1, wls_sum0),
//    )
//}

#[cfg(test)]
mod tests {
    use super::*;
    use genetics::alloc;

    // Error "SIGSEGV: invalid memory reference (signal 11)" will be raised if wzs_pad is not aligned

    fn initialize_wzs_wls(n: usize, val: f64) -> (Vec<f64>, Vec<f64>) {
        let v = vec![val; n];
        let wzs_pad = alloc::vec_align_f64(v, n + 32);
        let v = vec![val; n];
        let wls_pad = alloc::vec_align_f64(v, n + 32);

        //let mut wzs_pad = genetics::alloc::with_capacity_align::<f64>(n);
        //let mut wls_pad = genetics::alloc::with_capacity_align::<f64>(n);
        //wzs_pad.resize(n + 32, 0.0f64);
        //wls_pad.resize(n + 32, 0.0f64);
        //// fill values only for n
        //(0..n).zip(wzs_pad.iter_mut()).for_each(|(_, x)| *x = val);
        //(0..n).zip(wls_pad.iter_mut()).for_each(|(_, x)| *x = val);

        (wzs_pad, wls_pad)
    }

    fn initialize_wzs_wls_2(n: usize) -> (Vec<f64>, Vec<f64>) {
        let val1 = 0.01;
        let v = (0..n).map(|i| (i as f64) * val1).collect::<Vec<f64>>();
        let wzs_pad = alloc::vec_align_f64(v, n + 32);
        let val2 = 0.005;
        let v = (0..n).map(|i| (i as f64) * val2).collect::<Vec<f64>>();
        let wls_pad = alloc::vec_align_f64(v, n + 32);

        (wzs_pad, wls_pad)
    }

    #[test]
    fn test_calc_ws_sum_logit_nosimd_mi() {
        let vec = vec![0, 1, 1, 2, 1, 0];
        let gsnv = GenotSnv::new(&vec);

        let n = gsnv.n();
        assert_eq!(n, vec.len());

        let (wzs_pad, wls_pad) = initialize_wzs_wls(n, 0.1f64);

        let (wzs_sum, wls_sum) =
            calc_ws_sum_logit_nosimd_mi(&gsnv.as_genot_snv(), &wzs_pad, &wls_pad);

        assert_float_absolute_eq!(wzs_sum.0, 0.1);
        assert_float_absolute_eq!(wzs_sum.1, 0.3);
        assert_float_absolute_eq!(wzs_sum.2, 0.2);
        //assert_eq!(wzs_sum, (0.1, 0.3, 0.2));
        //assert_eq!(wzs_sum, (0.2, 0.3, 0.1));
        assert_eq!(wzs_sum, wls_sum);
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn test_calc_ws_sum_logit_nosimd_mi_vs_simd() {
        // check when n>32
        // n=33
        let mut vec = vec![0u8; 11];
        let vec1 = vec![1u8; 11];
        let vec2 = vec![2u8; 11];
        vec.extend(vec1);
        vec.extend(vec2);
        assert_eq!(vec.len(), 33);
        let gsnv = GenotSnv::new(&vec);

        let n = gsnv.n();
        assert_eq!(n, vec.len());

        let (wzs_pad, wls_pad) = initialize_wzs_wls_2(n);
        //let (wzs_pad, wls_pad) = initialize_wzs_wls(n, 0.1f64);

        if is_x86_feature_detected!("avx2") {
            unsafe {
                let ans_nosimd =
                    calc_ws_sum_logit_nosimd_mi(&gsnv.as_genot_snv(), &wzs_pad, &wls_pad);

                let ans_simd = calc_ws_sum_logit_simd_mi(&gsnv.as_genot_snv(), &wzs_pad, &wls_pad);

                assert_float_absolute_eq!(ans_nosimd.0 .0, ans_simd.0 .0);
                assert_float_absolute_eq!(ans_nosimd.0 .1, ans_simd.0 .1);
                assert_float_absolute_eq!(ans_nosimd.0 .2, ans_simd.0 .2);
                assert_float_absolute_eq!(ans_nosimd.1 .0, ans_simd.1 .0);
                assert_float_absolute_eq!(ans_nosimd.1 .1, ans_simd.1 .1);
                assert_float_absolute_eq!(ans_nosimd.1 .2, ans_simd.1 .2);
                //assert_eq!(ans_nosimd, ans_simd);
            }
        }
        // else do nth
    }

    #[test]
    fn test_calc_ws_sum_logit_interaction_nosimd_mi() {
        // all patterns
        let vec1 = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
        let vec2 = vec![0, 1, 2, 0, 1, 2, 0, 1, 2];
        let gsnv1 = GenotSnv::new(&vec1);
        let gsnv2 = GenotSnv::new(&vec2);

        let n = gsnv1.n();
        assert_eq!(n, vec1.len());
        assert_eq!(n, vec2.len());

        let (wzs_pad, wls_pad) = initialize_wzs_wls(n, 0.1f64);

        let (wzs_sum, wls_sum) = calc_ws_sum_logit_interaction_nosimd_mi(
            &gsnv1.as_genot_snv(),
            &gsnv2.as_genot_snv(),
            &wzs_pad,
            &wls_pad,
        );

        assert_eq!(
            wzs_sum.to_tuple(),
            ((0.1, 0.1, 0.1), (0.1, 0.1, 0.1), (0.1, 0.1, 0.1))
        );
        //assert_eq!(wzs_sum, (0.1, 0.2, 0.1, 0.5));
        //assert_eq!(wzs_sum, (0.5, 0.1, 0.2, 0.1));
        assert_eq!(wzs_sum.to_tuple(), wls_sum.to_tuple());
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn test_calc_ws_sum_logit_interaction_nosimd_mi_vs_simd() {
        // check when n>32
        // n=33
        let mut vec1 = vec![0u8; 11];
        let vec1_1 = vec![1u8; 11];
        let vec1_2 = vec![2u8; 11];
        vec1.extend(vec1_1);
        vec1.extend(vec1_2);

        let vec2 = vec![vec![0, 1, 2]; 11];
        let vec2 = vec2.into_iter().flatten().collect::<Vec<u8>>();

        assert_eq!(vec1.len(), 33);
        assert_eq!(vec2.len(), 33);

        let gsnv1 = GenotSnv::new(&vec1);
        let gsnv2 = GenotSnv::new(&vec2);

        let n = gsnv1.n();
        let (wzs_pad, wls_pad) = initialize_wzs_wls(n, 0.1f64);

        if is_x86_feature_detected!("avx2") {
            unsafe {
                let ans_nosimd = calc_ws_sum_logit_interaction_nosimd_mi(
                    &gsnv1.as_genot_snv(),
                    &gsnv2.as_genot_snv(),
                    &wzs_pad,
                    &wls_pad,
                );

                let ans_simd = calc_ws_sum_logit_interaction_simd_mi(
                    &gsnv1.as_genot_snv(),
                    &gsnv2.as_genot_snv(),
                    &wzs_pad,
                    &wls_pad,
                );

                let ans_nosimd_0 = ans_nosimd.0.to_tuple();
                let ans_nosimd_1 = ans_nosimd.1.to_tuple();

                let ans_simd_0 = ans_simd.0.to_tuple();
                let ans_simd_1 = ans_simd.1.to_tuple();

                assert_float_absolute_eq!(ans_nosimd_0.0 .0, ans_simd_0.0 .0);
                assert_float_absolute_eq!(ans_nosimd_0.0 .1, ans_simd_0.0 .1);
                assert_float_absolute_eq!(ans_nosimd_0.0 .2, ans_simd_0.0 .2);
                assert_float_absolute_eq!(ans_nosimd_0.1 .0, ans_simd_0.1 .0);
                assert_float_absolute_eq!(ans_nosimd_0.1 .1, ans_simd_0.1 .1);
                assert_float_absolute_eq!(ans_nosimd_0.1 .2, ans_simd_0.1 .2);
                assert_float_absolute_eq!(ans_nosimd_0.2 .0, ans_simd_0.2 .0);
                assert_float_absolute_eq!(ans_nosimd_0.2 .1, ans_simd_0.2 .1);
                assert_float_absolute_eq!(ans_nosimd_0.2 .2, ans_simd_0.2 .2);

                assert_float_absolute_eq!(ans_nosimd_1.0 .0, ans_simd_1.0 .0);
                assert_float_absolute_eq!(ans_nosimd_1.0 .1, ans_simd_1.0 .1);
                assert_float_absolute_eq!(ans_nosimd_1.0 .2, ans_simd_1.0 .2);
                assert_float_absolute_eq!(ans_nosimd_1.1 .0, ans_simd_1.1 .0);
                assert_float_absolute_eq!(ans_nosimd_1.1 .1, ans_simd_1.1 .1);
                assert_float_absolute_eq!(ans_nosimd_1.1 .2, ans_simd_1.1 .2);
                assert_float_absolute_eq!(ans_nosimd_1.2 .0, ans_simd_1.2 .0);
                assert_float_absolute_eq!(ans_nosimd_1.2 .1, ans_simd_1.2 .1);
                assert_float_absolute_eq!(ans_nosimd_1.2 .2, ans_simd_1.2 .2);

                //assert_float_absolute_eq!(ans_nosimd.0 .0, ans_simd.0 .0);
                //assert_float_absolute_eq!(ans_nosimd.0 .1, ans_simd.0 .1);
                //assert_float_absolute_eq!(ans_nosimd.0 .2, ans_simd.0 .2);
                //assert_float_absolute_eq!(ans_nosimd.0 .3, ans_simd.0 .3);
                //assert_float_absolute_eq!(ans_nosimd.1 .1, ans_simd.1 .1);
                //assert_float_absolute_eq!(ans_nosimd.1 .2, ans_simd.1 .2);
                //assert_float_absolute_eq!(ans_nosimd.1 .3, ans_simd.1 .3);
                //assert_eq!(ans_nosimd, ans_simd);
            }
        }
        // else do nth
    }

    //#[test]
    //fn test_calc_ws_sum_logit_interaction_nosimd_mi() {
    //    // all patterns
    //    let vec1 = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
    //    let vec2 = vec![0, 1, 2, 0, 1, 2, 0, 1, 2];
    //    let gsnv1 = GenotSnv::new(&vec1);
    //    let gsnv2 = GenotSnv::new(&vec2);

    //    let n = gsnv1.n();
    //    assert_eq!(n, vec1.len());
    //    assert_eq!(n, vec2.len());

    //    let (wzs_pad, wls_pad) = initialize_wzs_wls(n, 0.1f64);

    //    let (wzs_sum, wls_sum) = calc_ws_sum_logit_interaction_nosimd_mi(
    //        &gsnv1.as_genot_snv(),
    //        &gsnv2.as_genot_snv(),
    //        &wzs_pad,
    //        &wls_pad,
    //    );

    //    // g*g | number
    //    // 0    | 5
    //    // 1    | 1
    //    // 2    | 2
    //    // 4    | 1
    //    assert_eq!(wzs_sum, (0.1, 0.2, 0.1, 0.5));
    //    //assert_eq!(wzs_sum, (0.5, 0.1, 0.2, 0.1));
    //    assert_eq!(wzs_sum, wls_sum);
    //}

    //#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    //#[test]
    //fn test_calc_ws_sum_logit_interaction_nosimd_mi_vs_simd() {
    //    // check when n>32
    //    // n=33
    //    let mut vec1 = vec![0u8; 11];
    //    let vec1_1 = vec![1u8; 11];
    //    let vec1_2 = vec![2u8; 11];
    //    vec1.extend(vec1_1);
    //    vec1.extend(vec1_2);

    //    let vec2 = vec![vec![0, 1, 2]; 11];
    //    let vec2 = vec2.into_iter().flatten().collect::<Vec<u8>>();

    //    assert_eq!(vec1.len(), 33);
    //    assert_eq!(vec2.len(), 33);

    //    let gsnv1 = GenotSnv::new(&vec1);
    //    let gsnv2 = GenotSnv::new(&vec2);

    //    let n = gsnv1.n();
    //    let (wzs_pad, wls_pad) = initialize_wzs_wls(n, 0.1f64);

    //    if is_x86_feature_detected!("avx2") {
    //        unsafe {
    //            let ans_nosimd = calc_ws_sum_logit_interaction_nosimd_mi(
    //                &gsnv1.as_genot_snv(),
    //                &gsnv2.as_genot_snv(),
    //                &wzs_pad,
    //                &wls_pad,
    //            );

    //            let ans_simd = calc_ws_sum_logit_interaction_simd_mi(
    //                &gsnv1.as_genot_snv(),
    //                &gsnv2.as_genot_snv(),
    //                &wzs_pad,
    //                &wls_pad,
    //            );

    //            assert_float_absolute_eq!(ans_nosimd.0 .0, ans_simd.0 .0);
    //            assert_float_absolute_eq!(ans_nosimd.0 .1, ans_simd.0 .1);
    //            assert_float_absolute_eq!(ans_nosimd.0 .2, ans_simd.0 .2);
    //            assert_float_absolute_eq!(ans_nosimd.0 .3, ans_simd.0 .3);
    //            assert_float_absolute_eq!(ans_nosimd.1 .1, ans_simd.1 .1);
    //            assert_float_absolute_eq!(ans_nosimd.1 .2, ans_simd.1 .2);
    //            assert_float_absolute_eq!(ans_nosimd.1 .3, ans_simd.1 .3);
    //            //assert_eq!(ans_nosimd, ans_simd);
    //        }
    //    }
    //    // else do nth
    //}
}
