use crate::{genot::prelude::*, SampleScore};

// nomissing will not be used.
//pub fn scores_add_coef_nomissing(
//    scores: &mut SampleScore,
//    //scores_pad: &mut [f64],
//    score_wgt: (f64, f64, f64),
//    gsnv: &GenotSnvRef,
//) {
//    unsafe {
//        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
//        {
//            if is_x86_feature_detected!("avx2") {
//                let score_wgt_missing = (score_wgt.0, score_wgt.1, score_wgt.2, 0.0);
//                // DO NOT check if missing is included
//                scores_add_coef_simd(scores, score_wgt_missing, gsnv);
//                //scores_add_coef_simd(scores_pad, score_wgt_missing, genot_mi);
//                return;
//            }
//        }
//    }
//    scores_add_coef_nomissing_nosimd(scores, score_wgt, gsnv);
//    //scores_add_coef_nomissing_nosimd(scores_pad, score_wgt, genot_mi);
//    return;
//}
//
//// score_wgt: (s0, s1, s2)
//pub fn scores_add_coef_nomissing_nosimd(
//    scores: &mut SampleScore,
//    //scores: &mut [f64],
//    score_wgt: (f64, f64, f64),
//    gsnv: &GenotSnvRef,
//) {
//    let (s0, s1, s2) = score_wgt;
//    scores
//        .scores_mut()
//        .iter_mut()
//        .zip(gsnv.iter())
//        .for_each(|(score, val)| {
//            let score_add = match val {
//                0 => s0,
//                1 => s1,
//                2 => s2,
//                // panic for NA
//                //3 => sm,
//                _ => {
//                    panic!("Wrong genotype. Possibly NA for linear model.")
//                }
//            };
//            *score += score_add;
//        })
//}

pub fn scores_add_coef(
    scores: &mut SampleScore,
    score_wgt: (f64, f64, f64, f64),
    gsnv: &GenotSnvRef,
) {
    unsafe {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                scores_add_coef_simd(scores, score_wgt, gsnv);
                return;
            }
        }
    }
    scores_add_coef_nosimd(scores, score_wgt, gsnv);
    return;
}

// score_wgt: (s0, s1, s2, sm)
/// If sm=NAN but missing values exist, scores would be NAN and panic.
pub fn scores_add_coef_nosimd(
    scores: &mut SampleScore,
    //scores: &mut [f64],
    score_wgt: (f64, f64, f64, f64),
    gsnv: &GenotSnvRef,
) {
    let (s0, s1, s2, sm) = score_wgt;
    scores
        .scores_mut()
        .iter_mut()
        .zip(gsnv.iter())
        .for_each(|(score, val)| {
            let score_add = match val {
                0 => s0,
                1 => s1,
                2 => s2,
                3 => sm,
                _ => {
                    panic!("Wrong genotype. Possibly NA for linear model.")
                }
            };
            *score += score_add;
        });

    // in score.rs
    //scores.check_no_nan();
}

// TODO: common part should be in scores_add_coef()?-> troublesome due to $[cfg()]
/// If sm=NAN but missing values exist, scores would be NAN and panic.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn scores_add_coef_simd(
    scores: &mut SampleScore,
    //scores_pad: &mut [f64],
    score_wgt: (f64, f64, f64, f64),
    gsnv: &GenotSnvRef,
) {
    scores_add_coef_simd_inner(
        scores.scores_pad_mut(),
        score_wgt,
        gsnv.predict_s(0),
        gsnv.predict_s(1),
        gsnv.n(),
    );

    // make n..score.len() zero
    scores.clear_pad();
    //scores_pad[n..].iter_mut().for_each(|x| *x = 0.0f64);

    //scores.scores().iter().enumerate().for_each(|(i, x)| {
    //    if x.is_nan() {
    //        panic!("score contains nan at i: {}", i);
    //    }
    //});

    //scores.check_no_nan();

    //let (s0, s1, s2, sm) = score_wgt;

    //#[cfg(target_arch = "x86")]
    //use std::arch::x86::*;
    //#[cfg(target_arch = "x86_64")]
    //use std::arch::x86_64::*;

    //let n = gsnv.n();
    //let pred_s0m = gsnv.predict_s(0);
    //let pred_s1m = gsnv.predict_s(1);

    //let bit_ext_mask: __m256i = _mm256_set_epi32(
    //    0x10101010,
    //    0x01010101,
    //    0x20202020,
    //    0x02020202,
    //    0x40404040,
    //    0x04040404,
    //    -0x7f7f7f80, // =0x80808080
    //    0x08080808,
    //);
    ////let zerod: __m256d = _mm256_setzero_pd();
    //let zeros: __m256i = _mm256_setzero_si256();
    //let ones: __m256i = _mm256_cmpeq_epi32(zeros, zeros);

    //let s0d: __m256d = _mm256_set1_pd(s0);
    //let s1d: __m256d = _mm256_set1_pd(s1);
    //let s2d: __m256d = _mm256_set1_pd(s2);
    //let smd: __m256d = _mm256_set1_pd(sm);

    ////let mut l_acc: __m256d = _mm256_setzero_pd();

    //let shifts: [__m256i; 4] = [
    //    _mm256_set1_epi32(24),
    //    _mm256_set1_epi32(16),
    //    _mm256_set1_epi32(8),
    //    _mm256_set1_epi32(0),
    //];

    //for ni in 0..(n / 32 + 1) {
    //    //log::debug!("ni {}", ni);

    //    // broadcast 32bit int to 256bit
    //    // ex. DCBA -> DCBADCBA...DCBA
    //    // (D=abcdefgh)

    //    let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
    //    let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
    //    let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
    //    let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

    //    // bitwise not is fastest using xor by xor(v, ones)
    //    // sum2=( pred0 & pred1 )
    //    // sum1= ( !pred1 & pred0 )
    //    // sum0= ( !pred0 & (!pred1) )
    //    let flagv_2_32 = _mm256_and_si256(predv_s0_32, predv_s1_32);
    //    let flagv_1_32 = _mm256_andnot_si256(predv_s1_32, predv_s0_32);
    //    let flagv_0_32 = _mm256_andnot_si256(predv_s0_32, _mm256_xor_si256(predv_s1_32, ones));

    //    // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
    //    // abcdefgh(x4)|...
    //    // -> extracted (at highest position)
    //    // 000d0000(x4)|...
    //    // -> mask
    //    // 11111111|00000000|00000000|11111111|...
    //    let flagv_2_32_ext = _mm256_and_si256(flagv_2_32, bit_ext_mask);
    //    let flagv_1_32_ext = _mm256_and_si256(flagv_1_32, bit_ext_mask);
    //    let flagv_0_32_ext = _mm256_and_si256(flagv_0_32, bit_ext_mask);

    //    let take_mask_2_32 = _mm256_cmpeq_epi8(flagv_2_32_ext, bit_ext_mask);
    //    let take_mask_1_32 = _mm256_cmpeq_epi8(flagv_1_32_ext, bit_ext_mask);
    //    let take_mask_0_32 = _mm256_cmpeq_epi8(flagv_0_32_ext, bit_ext_mask);

    //    // bi=0-3, shift=24,16,8,0
    //    //const SHIFTS: &'static [i32] = &[24, 16, 8, 0];

    //    for bi in 0usize..4 {
    //        // DCBADCBA...DCBA
    //        // -> b=1: for B
    //        // BA00BA00...BA00

    //        let shift_v = shifts[bi];
    //        let take_mask_2 = _mm256_sllv_epi32(take_mask_2_32, shift_v);
    //        let take_mask_1 = _mm256_sllv_epi32(take_mask_1_32, shift_v);
    //        let take_mask_0 = _mm256_sllv_epi32(take_mask_0_32, shift_v);

    //        //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
    //        //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
    //        //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
    //        //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);

    //        let scores_lo_ptr = scores_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
    //        let scores_hi_ptr = scores_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
    //        //let mut scores_lo_ptr = scores[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_mut_ptr();
    //        //let mut scores_hi_ptr =
    //        //    scores[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_mut_ptr();

    //        //log::debug!("ps ind {}", 32 * ni + 8 * bi);
    //        //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
    //        //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
    //        //log::debug!(
    //        //    "ps hi {:?}",
    //        //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
    //        //);

    //        let scores_lo: __m256d = _mm256_load_pd(scores_lo_ptr as *const _);
    //        let scores_hi: __m256d = _mm256_load_pd(scores_hi_ptr as *const _);
    //        //let zsv_lo: __m256d = _mm256_load_pd(zsv_lo_ptr as *const _);
    //        //let zsv_hi: __m256d = _mm256_load_pd(zsv_hi_ptr as *const _);

    //        //log::debug!("ps lo {:?}", psv_lo);
    //        //log::debug!("ps hi {:?}", psv_hi);

    //        // to create f_lo, first set score0 and score1, and then set score2
    //        //let f_0_lo: __m256d= _mm256_blendv_pd(s0d, s1d, _mm256_castsi256_pd(take_mask_1));
    //        let f_0_lo: __m256d = _mm256_blendv_pd(smd, s0d, _mm256_castsi256_pd(take_mask_0));
    //        //let f_0_lo: __m256d = _mm256_blendv_pd(zerod, s0d, _mm256_castsi256_pd(take_mask_0));
    //        let f_01_lo: __m256d = _mm256_blendv_pd(f_0_lo, s1d, _mm256_castsi256_pd(take_mask_1));
    //        let f_lo: __m256d = _mm256_blendv_pd(f_01_lo, s2d, _mm256_castsi256_pd(take_mask_2));

    //        let scores_lo_add: __m256d = _mm256_add_pd(scores_lo, f_lo);
    //        _mm256_store_pd(scores_lo_ptr as *mut _, scores_lo_add);

    //        // for high
    //        let take_mask_2_hi = _mm256_slli_epi64(take_mask_2, 32);
    //        let take_mask_1_hi = _mm256_slli_epi64(take_mask_1, 32);
    //        let take_mask_0_hi = _mm256_slli_epi64(take_mask_0, 32);

    //        let f_0_hi: __m256d = _mm256_blendv_pd(smd, s0d, _mm256_castsi256_pd(take_mask_0_hi));
    //        //let f_0_hi: __m256d = _mm256_blendv_pd(zerod, s0d, _mm256_castsi256_pd(take_mask_0_hi));
    //        let f_01_hi: __m256d =
    //            _mm256_blendv_pd(f_0_hi, s1d, _mm256_castsi256_pd(take_mask_1_hi));
    //        let f_hi: __m256d = _mm256_blendv_pd(f_01_hi, s2d, _mm256_castsi256_pd(take_mask_2_hi));

    //        let scores_hi_add: __m256d = _mm256_add_pd(scores_hi, f_hi);
    //        _mm256_store_pd(scores_hi_ptr as *mut _, scores_hi_add);
    //    }
    //}

    //// make n..score.len() zero
    //scores_pad[n..].iter_mut().for_each(|x| *x = 0.0f64);
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn scores_add_coef_simd_inner(
    scores_pad: &mut [f64],
    score_wgt: (f64, f64, f64, f64),
    gsnv_s0: &[u8],
    gsnv_s1: &[u8],
    n: usize,
    //gsnv: &GenotSnvRef,
) {
    let (s0, s1, s2, sm) = score_wgt;

    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    //let n = gsnv.n();
    //let pred_s0m = gsnv.predict_s(0);
    //let pred_s1m = gsnv.predict_s(1);

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
    let zeros: __m256i = _mm256_setzero_si256();
    let ones: __m256i = _mm256_cmpeq_epi32(zeros, zeros);

    let s0d: __m256d = _mm256_set1_pd(s0);
    let s1d: __m256d = _mm256_set1_pd(s1);
    let s2d: __m256d = _mm256_set1_pd(s2);
    let smd: __m256d = _mm256_set1_pd(sm);

    //let mut l_acc: __m256d = _mm256_setzero_pd();

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

        let pred_s0_b32 = u32::from_le_bytes(gsnv_s0[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(gsnv_s1[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        // bitwise not is fastest using xor by xor(v, ones)
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

            let scores_lo_ptr = scores_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let scores_hi_ptr = scores_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            //let mut scores_lo_ptr = scores[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_mut_ptr();
            //let mut scores_hi_ptr =
            //    scores[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_mut_ptr();

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            let scores_lo: __m256d = _mm256_load_pd(scores_lo_ptr as *const _);
            let scores_hi: __m256d = _mm256_load_pd(scores_hi_ptr as *const _);
            //let zsv_lo: __m256d = _mm256_load_pd(zsv_lo_ptr as *const _);
            //let zsv_hi: __m256d = _mm256_load_pd(zsv_hi_ptr as *const _);

            //log::debug!("ps lo {:?}", psv_lo);
            //log::debug!("ps hi {:?}", psv_hi);

            // first set smd and overwrite with s0d, s1d, s2d
            let f_0_lo: __m256d = _mm256_blendv_pd(smd, s0d, _mm256_castsi256_pd(take_mask_0));
            //let f_0_lo: __m256d= _mm256_blendv_pd(s0d, s1d, _mm256_castsi256_pd(take_mask_1));
            //let f_0_lo: __m256d = _mm256_blendv_pd(zerod, s0d, _mm256_castsi256_pd(take_mask_0));
            let f_01_lo: __m256d = _mm256_blendv_pd(f_0_lo, s1d, _mm256_castsi256_pd(take_mask_1));
            let f_lo: __m256d = _mm256_blendv_pd(f_01_lo, s2d, _mm256_castsi256_pd(take_mask_2));

            let scores_lo_add: __m256d = _mm256_add_pd(scores_lo, f_lo);
            _mm256_store_pd(scores_lo_ptr as *mut _, scores_lo_add);

            // for high
            let take_mask_2_hi = _mm256_slli_epi64(take_mask_2, 32);
            let take_mask_1_hi = _mm256_slli_epi64(take_mask_1, 32);
            let take_mask_0_hi = _mm256_slli_epi64(take_mask_0, 32);

            let f_0_hi: __m256d = _mm256_blendv_pd(smd, s0d, _mm256_castsi256_pd(take_mask_0_hi));
            //let f_0_hi: __m256d = _mm256_blendv_pd(zerod, s0d, _mm256_castsi256_pd(take_mask_0_hi));
            let f_01_hi: __m256d =
                _mm256_blendv_pd(f_0_hi, s1d, _mm256_castsi256_pd(take_mask_1_hi));
            let f_hi: __m256d = _mm256_blendv_pd(f_01_hi, s2d, _mm256_castsi256_pd(take_mask_2_hi));

            let scores_hi_add: __m256d = _mm256_add_pd(scores_hi, f_hi);
            _mm256_store_pd(scores_hi_ptr as *mut _, scores_hi_add);
        }
    }
}

// TODO: test
/// allow missing
pub fn scores_add_coef_interaction_nosimd(
    scores: &mut SampleScore,
    score_wgt: [[f64; 4]; 4],
    //score_wgt: ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)),
    gsnv_1: &GenotSnvRef,
    gsnv_2: &GenotSnvRef,
) {
    //let (constt, alphat) = score_wgt;
    //let (s0, s1, s2, s4) = score_wgt;

    scores
        .scores_mut()
        .iter_mut()
        .zip(gsnv_1.iter().zip(gsnv_2.iter()))
        .for_each(|(score, (g1, g2))| {
            let score_add = match (g1, g2) {
                (0, 0) => score_wgt[0][0],
                (0, 1) => score_wgt[0][1],
                (0, 2) => score_wgt[0][2],
                (0, 3) => score_wgt[0][3],
                (1, 0) => score_wgt[1][0],
                (1, 1) => score_wgt[1][1],
                (1, 2) => score_wgt[1][2],
                (1, 3) => score_wgt[1][3],
                (2, 0) => score_wgt[2][0],
                (2, 1) => score_wgt[2][1],
                (2, 2) => score_wgt[2][2],
                (2, 3) => score_wgt[2][3],
                (3, 0) => score_wgt[3][0],
                (3, 1) => score_wgt[3][1],
                (3, 2) => score_wgt[3][2],
                (3, 3) => score_wgt[3][3],
                _ => panic!("wrong genotype pair: {}, {}", g1, g2),
            };

            *score += score_add;
        })
}

pub fn scores_add_coef_interaction_no_missing(
    scores: &mut SampleScore,
    score_wgt: ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)),
    gsnv_1: &GenotSnvRef,
    gsnv_2: &GenotSnvRef,
) {
    unsafe {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                scores_add_coef_interaction_no_missing_simd(scores, score_wgt, gsnv_1, gsnv_2);
                return;
            }
        }
    }
    scores_add_coef_interaction_no_missing_nosimd(scores, score_wgt, gsnv_1, gsnv_2);
    return;
}

pub fn scores_add_coef_interaction_no_missing_nosimd(
    scores: &mut SampleScore,
    score_wgt: ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)),
    gsnv_1: &GenotSnvRef,
    gsnv_2: &GenotSnvRef,
) {
    //let (constt, alphat) = score_wgt;
    //let (s0, s1, s2, s4) = score_wgt;

    scores
        .scores_mut()
        .iter_mut()
        .zip(gsnv_1.iter().zip(gsnv_2.iter()))
        .for_each(|(score, (g1, g2))| {
            let score_add = match (g1, g2) {
                (0, 0) => score_wgt.0 .0,
                (0, 1) => score_wgt.0 .1,
                (0, 2) => score_wgt.0 .2,
                (1, 0) => score_wgt.1 .0,
                (1, 1) => score_wgt.1 .1,
                (1, 2) => score_wgt.1 .2,
                (2, 0) => score_wgt.2 .0,
                (2, 1) => score_wgt.2 .1,
                (2, 2) => score_wgt.2 .2,
                _ => panic!("wrong genotype pair: {}, {}", g1, g2),
            };

            *score += score_add;
        })
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn scores_add_coef_interaction_no_missing_simd(
    scores: &mut SampleScore,
    //scores_pad: &mut [f64],
    score_wgt: ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)),
    gsnv1: &GenotSnvRef,
    gsnv2: &GenotSnvRef,
) {
    scores_add_coef_simd_interaction_inner(
        scores.scores_pad_mut(),
        score_wgt,
        gsnv1.predict_s(0),
        gsnv1.predict_s(1),
        gsnv2.predict_s(0),
        gsnv2.predict_s(1),
        gsnv1.n(),
    );

    // make n..score.len() zero
    scores.clear_pad();

    //scores.check_no_nan();
}

// score is NAN for missing values
// -> cannot be used for score
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn scores_add_coef_simd_interaction_inner(
    scores_pad: &mut [f64],
    score_wgt: ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)),
    gsnv1_s0: &[u8],
    gsnv1_s1: &[u8],
    gsnv2_s0: &[u8],
    gsnv2_s1: &[u8],
    n: usize,
    //gsnv: &GenotSnvRef,
) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let ((s00, s01, s02), (s10, s11, s12), (s20, s21, s22)) = score_wgt;

    let s_m_d: __m256d = _mm256_set1_pd(f64::NAN);
    let s_0_0_d: __m256d = _mm256_set1_pd(s00);
    let s_0_1_d: __m256d = _mm256_set1_pd(s01);
    let s_0_2_d: __m256d = _mm256_set1_pd(s02);
    let s_1_0_d: __m256d = _mm256_set1_pd(s10);
    let s_1_1_d: __m256d = _mm256_set1_pd(s11);
    let s_1_2_d: __m256d = _mm256_set1_pd(s12);
    let s_2_0_d: __m256d = _mm256_set1_pd(s20);
    let s_2_1_d: __m256d = _mm256_set1_pd(s21);
    let s_2_2_d: __m256d = _mm256_set1_pd(s22);

    //let n = gsnv.n();
    //let pred_s0m = gsnv.predict_s(0);
    //let pred_s1m = gsnv.predict_s(1);

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
    let zeros: __m256i = _mm256_setzero_si256();
    let ones: __m256i = _mm256_cmpeq_epi32(zeros, zeros);

    //let mut l_acc: __m256d = _mm256_setzero_pd();

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

        let pred1_s0_b32 = u32::from_le_bytes(gsnv1_s0[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred1_s1_b32 = u32::from_le_bytes(gsnv1_s1[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred2_s0_b32 = u32::from_le_bytes(gsnv2_s0[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred2_s1_b32 = u32::from_le_bytes(gsnv2_s1[4 * ni..4 * (ni + 1)].try_into().unwrap());

        let predv1_s0_32 = _mm256_set1_epi32(pred1_s0_b32 as i32);
        let predv1_s1_32 = _mm256_set1_epi32(pred1_s1_b32 as i32);
        let predv2_s0_32 = _mm256_set1_epi32(pred2_s0_b32 as i32);
        let predv2_s1_32 = _mm256_set1_epi32(pred2_s1_b32 as i32);

        // bitwise not is fastest using xor by xor(v, ones)
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

            //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
            //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
            //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
            //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);

            let scores_lo_ptr = scores_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let scores_hi_ptr = scores_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();
            //let mut scores_lo_ptr = scores[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_mut_ptr();
            //let mut scores_hi_ptr =
            //    scores[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_mut_ptr();

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            let scores_lo: __m256d = _mm256_load_pd(scores_lo_ptr as *const _);
            let scores_hi: __m256d = _mm256_load_pd(scores_hi_ptr as *const _);
            //let zsv_lo: __m256d = _mm256_load_pd(zsv_lo_ptr as *const _);
            //let zsv_hi: __m256d = _mm256_load_pd(zsv_hi_ptr as *const _);

            //log::debug!("ps lo {:?}", psv_lo);
            //log::debug!("ps hi {:?}", psv_hi);

            // first set smd and overwrite with s0d, s1d, s2d
            let f_0_0_lo: __m256d =
                _mm256_blendv_pd(s_m_d, s_0_0_d, _mm256_castsi256_pd(take_mask_0_0));
            let f_0_1_lo: __m256d =
                _mm256_blendv_pd(f_0_0_lo, s_0_1_d, _mm256_castsi256_pd(take_mask_0_1));
            let f_0_2_lo: __m256d =
                _mm256_blendv_pd(f_0_1_lo, s_0_2_d, _mm256_castsi256_pd(take_mask_0_2));
            let f_1_0_lo: __m256d =
                _mm256_blendv_pd(f_0_2_lo, s_1_0_d, _mm256_castsi256_pd(take_mask_1_0));
            let f_1_1_lo: __m256d =
                _mm256_blendv_pd(f_1_0_lo, s_1_1_d, _mm256_castsi256_pd(take_mask_1_1));
            let f_1_2_lo: __m256d =
                _mm256_blendv_pd(f_1_1_lo, s_1_2_d, _mm256_castsi256_pd(take_mask_1_2));
            let f_2_0_lo: __m256d =
                _mm256_blendv_pd(f_1_2_lo, s_2_0_d, _mm256_castsi256_pd(take_mask_2_0));
            let f_2_1_lo: __m256d =
                _mm256_blendv_pd(f_2_0_lo, s_2_1_d, _mm256_castsi256_pd(take_mask_2_1));
            let f_lo: __m256d =
                _mm256_blendv_pd(f_2_1_lo, s_2_2_d, _mm256_castsi256_pd(take_mask_2_2));

            let scores_lo_add: __m256d = _mm256_add_pd(scores_lo, f_lo);
            _mm256_store_pd(scores_lo_ptr as *mut _, scores_lo_add);

            // for high
            let take_mask_0_0_hi = _mm256_slli_epi64(take_mask_0_0, 32);
            let f_0_0_hi: __m256d =
                _mm256_blendv_pd(s_m_d, s_0_0_d, _mm256_castsi256_pd(take_mask_0_0_hi));
            let take_mask_0_1_hi = _mm256_slli_epi64(take_mask_0_1, 32);
            let f_0_1_hi: __m256d =
                _mm256_blendv_pd(f_0_0_hi, s_0_1_d, _mm256_castsi256_pd(take_mask_0_1_hi));
            let take_mask_0_2_hi = _mm256_slli_epi64(take_mask_0_2, 32);
            let f_0_2_hi: __m256d =
                _mm256_blendv_pd(f_0_1_hi, s_0_2_d, _mm256_castsi256_pd(take_mask_0_2_hi));
            let take_mask_1_0_hi = _mm256_slli_epi64(take_mask_1_0, 32);
            let f_1_0_hi: __m256d =
                _mm256_blendv_pd(f_0_2_hi, s_1_0_d, _mm256_castsi256_pd(take_mask_1_0_hi));
            let take_mask_1_1_hi = _mm256_slli_epi64(take_mask_1_1, 32);
            let f_1_1_hi: __m256d =
                _mm256_blendv_pd(f_1_0_hi, s_1_1_d, _mm256_castsi256_pd(take_mask_1_1_hi));
            let take_mask_1_2_hi = _mm256_slli_epi64(take_mask_1_2, 32);
            let f_1_2_hi: __m256d =
                _mm256_blendv_pd(f_1_1_hi, s_1_2_d, _mm256_castsi256_pd(take_mask_1_2_hi));
            let take_mask_2_0_hi = _mm256_slli_epi64(take_mask_2_0, 32);
            let f_2_0_hi: __m256d =
                _mm256_blendv_pd(f_1_2_hi, s_2_0_d, _mm256_castsi256_pd(take_mask_2_0_hi));
            let take_mask_2_1_hi = _mm256_slli_epi64(take_mask_2_1, 32);
            let f_2_1_hi: __m256d =
                _mm256_blendv_pd(f_2_0_hi, s_2_1_d, _mm256_castsi256_pd(take_mask_2_1_hi));
            let take_mask_2_2_hi = _mm256_slli_epi64(take_mask_2_2, 32);
            let f_hi: __m256d =
                _mm256_blendv_pd(f_2_1_hi, s_2_2_d, _mm256_castsi256_pd(take_mask_2_2_hi));

            let scores_hi_add: __m256d = _mm256_add_pd(scores_hi, f_hi);
            _mm256_store_pd(scores_hi_ptr as *mut _, scores_hi_add);
        }
    }
}

//// score_wgt: (s0, s1, s2, s4)
//pub fn scores_add_coef_interaction(
//    scores: &mut [f64],
//    score_wgt: (f64, f64, f64, f64),
//    gsnv_1: &GenotSnvRef,
//    gsnv_2: &GenotSnvRef,
//) {
//    let (s0, s1, s2, s4) = score_wgt;
//    scores
//        .iter_mut()
//        .zip(gsnv_1.iter().zip(gsnv_2.iter()))
//        .for_each(|(score, (val_1, val_2))| {
//            let score_add = match val_1 * val_2 {
//                0 => s0,
//                1 => s1,
//                2 => s2,
//                4 => s4,
//                // panic for NA
//                //3 => sm,
//                _ => {
//                    panic!("Wrong genotype. Possibly NA for linear model.")
//                }
//            };
//            *score += score_add;
//        })
//}

#[cfg(test)]
mod tests {
    use super::*;

    //#[test]
    //fn test_scores_add_coef_nomissing_nosimd() {
    //    let vec = vec![0, 1, 1, 2, 1, 0];
    //    let g = GenotSnv::new(&vec);
    //    let g_snv_ref = g.as_genot_snv();
    //    let scores = vec![0.0f64, 1.0, 2.0, 3.0, 4.0, 5.0];
    //    let mut scores = SampleScore::new_vec(scores);

    //    let score_wgt = (0.0, 1.0, 3.0);
    //    scores_add_coef_nomissing_nosimd(&mut scores, score_wgt, &g_snv_ref);

    //    let scores_exp = vec![0.0f64, 2.0, 3.0, 6.0, 5.0, 5.0];
    //    assert_eq!(scores.scores(), &scores_exp);
    //}

    //// include missing
    //#[test]
    //#[should_panic]
    //fn test_scores_add_coef_nomissing_nosimd_panic() {
    //    let vec = vec![0, 1, 1, 2, 1, 3];
    //    let g = GenotSnv::new(&vec);
    //    let g_snv_ref = g.as_genot_snv();
    //    let scores = vec![0.0f64, 1.0, 2.0, 3.0, 4.0, 5.0];
    //    let mut scores = SampleScore::new_vec(scores);

    //    let score_wgt = (0.0, 1.0, 3.0);
    //    scores_add_coef_nomissing_nosimd(&mut scores, score_wgt, &g_snv_ref);
    //}

    // include missing
    // -> no check_nan here
    //#[test]
    //#[should_panic]
    //fn test_scores_add_coef_nosimd_panic() {
    //    let vec = vec![0, 1, 1, 2, 1, 3];
    //    let g = GenotSnv::new(&vec);
    //    let g_snv_ref = g.as_genot_snv();
    //    let scores = vec![0.0f64, 1.0, 2.0, 3.0, 4.0, 5.0];
    //    let mut scores = SampleScore::new_vec(scores);

    //    let score_wgt = (0.0, 1.0, 3.0, f64::NAN);
    //    scores_add_coef_nosimd(&mut scores, score_wgt, &g_snv_ref);
    //}

    // -> no check_nan here
    //#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    //#[test]
    //#[should_panic]
    //fn test_scores_add_coef_simd_panic() {
    //    let vec = vec![0, 1, 1, 2, 1, 3];
    //    let g = GenotSnv::new(&vec);
    //    let g_snv_ref = g.as_genot_snv();
    //    let scores = vec![0.0f64, 1.0, 2.0, 3.0, 4.0, 5.0];
    //    let mut scores = SampleScore::new_vec(scores);

    //    let score_wgt = (0.0, 1.0, 3.0, f64::NAN);

    //    if is_x86_feature_detected!("avx2") {
    //        unsafe {
    //            scores_add_coef_simd(&mut scores, score_wgt, &g_snv_ref);
    //        }
    //    }
    //}

    #[test]
    fn test_scores_add_coef_nosimd() {
        let vec = vec![0, 1, 1, 2, 1, 3];
        let g = GenotSnv::new(&vec);
        let g_snv_ref = g.as_genot_snv();
        let scores = vec![0.0f64, 1.0, 2.0, 3.0, 4.0, 5.0];
        let mut scores = SampleScore::new_vec(scores);

        let score_wgt = (0.0, 1.0, 3.0, 5.0);
        scores_add_coef_nosimd(&mut scores, score_wgt, &g_snv_ref);

        let scores_exp = vec![0.0f64, 2.0, 3.0, 6.0, 5.0, 10.0];
        assert_eq!(scores.scores(), &scores_exp);
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn test_scores_add_coef_nosimd_vs_simd() {
        // check when n>32
        // n=33

        let vec = vec![vec![0, 1, 2]; 11];
        let vec = vec.into_iter().flatten().collect::<Vec<u8>>();

        assert_eq!(vec.len(), 33);

        let gsnv = GenotSnv::new(&vec);
        let n = gsnv.n();

        let score_wgt = (0.2, 0.4, 0.8, 0.1);

        if is_x86_feature_detected!("avx2") {
            unsafe {
                let mut scores_nosimd = SampleScore::new(n);
                scores_add_coef_nosimd(&mut scores_nosimd, score_wgt, &gsnv.as_genot_snv());

                let mut scores_simd = SampleScore::new(n);
                scores_add_coef_simd(&mut scores_simd, score_wgt, &gsnv.as_genot_snv());

                let scores_nosimd_pad = scores_nosimd.scores_pad();
                let scores_simd_pad = scores_simd.scores_pad();
                for i in 0..n {
                    assert_float_absolute_eq!(scores_nosimd_pad[i], scores_simd_pad[i]);
                }
                // 0.0 for pad
                for i in n..scores_simd_pad.len() {
                    assert_float_absolute_eq!(scores_simd_pad[i], 0.0f64);
                    assert_float_absolute_eq!(scores_nosimd_pad[i], 0.0f64);
                }
            }
        } //else {
          // panic!("not avx2")
          //}
          //else do nth
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn test_scores_interaction_nosimd_mi_vs_simd() {
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

        let score_wgt = ((0.1, 0.2, 0.3), (0.4, 0.5, 0.6), (0.7, 0.8, 0.9));

        let n = gsnv1.n();

        if is_x86_feature_detected!("avx2") {
            unsafe {
                let mut scores_nosimd = SampleScore::new(n);
                scores_add_coef_interaction_no_missing_nosimd(
                    &mut scores_nosimd,
                    score_wgt,
                    &gsnv1.as_genot_snv(),
                    &gsnv2.as_genot_snv(),
                );

                let mut scores_simd = SampleScore::new(n);
                scores_add_coef_interaction_no_missing_simd(
                    &mut scores_simd,
                    score_wgt,
                    &gsnv1.as_genot_snv(),
                    &gsnv2.as_genot_snv(),
                );

                let scores_nosimd_pad = scores_nosimd.scores_pad();
                let scores_simd_pad = scores_simd.scores_pad();
                for i in 0..n {
                    assert_float_absolute_eq!(scores_nosimd_pad[i], scores_simd_pad[i]);
                }
                // 0.0 for pad
                for i in n..scores_simd_pad.len() {
                    assert_float_absolute_eq!(scores_simd_pad[i], 0.0f64);
                    assert_float_absolute_eq!(scores_nosimd_pad[i], 0.0f64);
                }
            }
        } //else {
          // panic!("not avx2")
    }
}
