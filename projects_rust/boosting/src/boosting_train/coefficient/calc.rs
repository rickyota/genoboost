use genetics::genot::prelude::*;

pub unsafe fn calculate_coef_gt_logit_sm(
    gsnv: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
    //phe: &Phe,
    //epsilons_wls: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //epsilons_wzs: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //eps: Eps,
) -> ((f64, f64, f64), (f64, f64, f64)) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            //let (wzs_sum, wls_sum) = calc::calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
            return calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
        }
    }
    return calculate_coef_gt_logit_nosimd_sm(gsnv, wzs_pad, wls_pad);
}

// TODO: move to coef.rs
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn calculate_coef_gt_logit_simd_sm(
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

    // n + 32
    //let n_pad=wzs.len();
    let n = wzs_pad.len() - 32;
    //let n = phe.n();
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

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            let wzsv_lo: __m256d = _mm256_load_pd(wzsv_lo_ptr as *const _);
            let wzsv_hi: __m256d = _mm256_load_pd(wzsv_hi_ptr as *const _);
            let wlsv_lo: __m256d = _mm256_load_pd(wlsv_lo_ptr as *const _);
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

pub unsafe fn calculate_coef_gt_logit_nosimd_sm(
    _gsnv: &GenotSnvRef,
    _wzs_pad: &[f64],
    _wls_pad: &[f64],
    //phe: &Phe,
    //epsilons_wls: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //epsilons_wzs: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //eps: Eps,
) -> ((f64, f64, f64), (f64, f64, f64)) {
    unimplemented!()
}
