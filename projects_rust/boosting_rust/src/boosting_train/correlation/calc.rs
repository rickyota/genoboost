// TODO: many of commonly used ones should be moved to genetics_rust

#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use std::convert::TryInto;

//use super::super::mistake;
use super::WgtBoost;
use genetics::genot::prelude::*;
use genetics::{alloc, bit, genot_twin, vec, GenotTwin, SnvId, B8};
use rayon::prelude::*;

fn convert_predicts_to_genotype(pred0: bool, pred1: bool) -> u8 {
    unimplemented!("not considering missing");
    pred0 as u8 + pred1 as u8
}

pub fn calculate_mean_m(predicts: &GenotTwinSnvRef) -> f64 {
    let sum = predicts.iter().map(|x| x as f64).sum::<f64>();

    let mean = sum / (predicts.n() as f64);

    mean
}

pub fn calculate_var_m(mean: f64, predicts: &GenotTwinSnvRef) -> f64 {

    let sum2 = predicts
        .iter()
        .map(|x| (x as f64) * (x as f64))
        .sum::<f64>();

    let n = predicts.n();
    let var = sum2 / (n as f64) - mean * mean;

    var
}

/*
fn get_len_n_m256i(n: usize) -> usize {
    n / 32 + 5
}

/// overflowed predict should be 0.
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
unsafe fn create_gs_epi8(predicts: &[B8], n: usize) -> Vec<__m256i> {
    let len_n = get_len_n_m256i(n);
    let mut gs_epi8 = alloc::with_capacity_align_m256i(len_n);

    let len_n_p = genot_twin::len_n(n);
    //let len_n_p = genotype::len_n(n);
    let pred_s0m = &predicts[..len_n_p];
    let pred_s1m = &predicts[len_n_p..];

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
    let ones_8 = _mm256_set1_epi8(1);

    for ni in 0..(n / 32 + 1) {
        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
        // abcdefgh(x4)|...
        // -> extracted (at highest position)
        // 000d0000(x4)|...
        // -> mask
        // 11111111|00000000|00000000|11111111|...
        let flagv_s0_32_ext = _mm256_and_si256(predv_s0_32, bit_ext_mask);
        let flagv_s1_32_ext = _mm256_and_si256(predv_s1_32, bit_ext_mask);

        let take_mask_s0_32 = _mm256_cmpeq_epi8(flagv_s0_32_ext, bit_ext_mask);
        let take_mask_s1_32 = _mm256_cmpeq_epi8(flagv_s1_32_ext, bit_ext_mask);

        let xs_s0_8 = _mm256_and_si256(take_mask_s0_32, ones_8);
        let xs_s1_8 = _mm256_and_si256(take_mask_s1_32, ones_8);
        let gs_8 = _mm256_add_epi8(xs_s0_8, xs_s1_8);

        gs_epi8.push(gs_8);
    }

    gs_epi8
}

#[target_feature(enable = "avx2")]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
unsafe fn calculate_dot(predicts: &[B8], gs_chosen: &[__m256i], n: usize) -> usize {
    let gs: Vec<__m256i> = create_gs_epi8(predicts, n);

    let zeros = _mm256_set1_epi8(0);

    let mut dotsum = _mm256_setzero_si256();
    for ni in 0..(n / 32 + 1) {
        let gs_8 = gs[ni];
        let gs_chosen_8 = gs_chosen[ni];

        let dot_8_lo = _mm256_mullo_epi16(gs_8, gs_chosen_8);
        let dot_8_hi = _mm256_mullo_epi16(
            _mm256_srli_epi16(gs_8, 8),
            _mm256_srli_epi16(gs_chosen_8, 8),
        );
        let dot_8_sum = _mm256_add_epi8(dot_8_lo, dot_8_hi);
        let dot_16_sum = _mm256_and_si256(dot_8_sum, _mm256_set1_epi16(0xFF));
        let dot_32_sum = _mm256_add_epi32(
            _mm256_unpackhi_epi16(dot_16_sum, zeros),
            _mm256_unpacklo_epi16(dot_16_sum, zeros),
        );
        dotsum = _mm256_add_epi32(dotsum, dot_32_sum);
        /*
        __m256i const dot_8_lo = _mm256_mullo_epi16(gs_8, gs_chosen_8);
        __m256i const dot_8_hi = _mm256_mullo_epi16(_mm256_srli_epi16(gs_8, 8), _mm256_srli_epi16(gs_chosen_8, 8));
        __m256i const dot_8_sum = _mm256_add_epi8(dot_8_lo, dot_8_hi);
        __m256i const dot_16_sum = _mm256_and_si256(dot_8_sum, _mm256_set1_epi16(0xFF));
        __m256i const dot_32_sum = _mm256_add_epi32(_mm256_unpackhi_epi16(dot_16_sum, zeros), _mm256_unpacklo_epi16(dot_16_sum, zeros));
        dotsum = _mm256_add_epi32(dotsum, dot_32_sum);
        */
    }

    // sum 4 double horizontally to get the whole sum
    dotsum = _mm256_hadd_epi32(dotsum, dotsum);
    dotsum = _mm256_hadd_epi32(dotsum, dotsum);
    // cannot use since dotsum could be epi32x4 or epi8x32
    // int const dotsum_m=dotsum[0] + dotsum[4];
    //int const dotsum_m = _mm256_extract_epi32(dotsum, 0) + _mm256_extract_epi32(dotsum, 4);
    //return dotsum_m;

    let dotsum_m =
        _mm256_cvtsi256_si32(dotsum) + _mm_cvtsi128_si32(_mm256_extracti128_si256(dotsum, 1));
    //let dotsum_m = _mm256_cvtsd_f64(dotsum) + _mm_cvtsd_f64(_mm256_extractf128_pd(a_sum_s0_acc, 1));

    dotsum_m as usize
}

#[target_feature(enable = "avx2")]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
unsafe fn calculate_r2_simd(
    predicts: &[B8],
    gs_chosen: &[__m256i],
    //predicts_chosen: &[B8],
    gt_mean: f64,
    gt_var: f64,
    gt_mean_chosen: f64,
    gt_var_chosen: f64,
    n: usize,
) -> f64 {
    let dot = calculate_dot(predicts, gs_chosen, n);

    /*
    let mut dot = 0.0f64;
    for ni in 0..n {
        let g = convert_predicts_to_genotype(operate::bget(preds0, ni), operate::bget(preds1, ni));
        let g_c =
            convert_predicts_to_genotype(operate::bget(preds0_c, ni), operate::bget(preds1_c, ni));
        dot += (g as f64) * (g_c as f64);
    }
    */

    let r2 = ((dot as f64) / (n as f64) - gt_mean * gt_mean_chosen)
        / (gt_var.sqrt() * gt_var_chosen.sqrt());
    r2
}

// Now, snv_indexs can be &[SnvIndex], &[Snv] or &[SumStat].
#[target_feature(enable = "avx2")]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub unsafe fn calculate_r2s_simd<T>(
    r2s: &mut [f64],
    predictions: &GenotTwin,
    gt_means: &[f64],
    gt_vars: &[f64],
    is_pruned: &[bool],
    wgt_chosen: &WgtBoost,
    snv_indexs: &[T],
    //snv_indexs: &[SnvIndex],
    n: usize,
) where
    T: AsRef<SnvId> + std::marker::Sync, //Sync is required by rayon
{
    vec::fill(r2s, 0.0);

    let mi_chosen = wgt_chosen.wgt().kind().index_snv().unwrap();
    let chrom_chosen = wgt_chosen.wgt().kind().snv_index().chrom();

    let gs_chosen: Vec<__m256i> = create_gs_epi8(predictions.predictions_snv_s(mi_chosen), n);
    //create_gs_epi8(predict::predictions_snv_s(predictions, mi_chosen, n), n);
    //let gs_chosen: Vec<__m256d> = _mm256_setzero_pd();

    r2s.par_iter_mut()
        .enumerate()
        .filter(|(mi, _)| !is_pruned[*mi] && snv_indexs[*mi].as_ref().chrom() == chrom_chosen)
        .for_each(|(mi, r2)| {
            *r2 = calculate_r2_simd(
                predictions.predictions_snv_s(mi),
                //predict::predictions_snv_s(predictions, mi, n),
                &gs_chosen,
                gt_means[mi],
                gt_vars[mi],
                gt_means[mi_chosen],
                gt_vars[mi_chosen],
                n,
            )
        });
}

fn calculate_r2_nosimd(
    predicts: &[B8],
    predicts_chosen: &[B8],
    gt_mean: f64,
    gt_var: f64,
    gt_mean_chosen: f64,
    gt_var_chosen: f64,
    n: usize,
) -> f64 {
    let len_n = genot_twin::len_n(n);
    //let len_n = genotype::len_n(n);
    let preds0 = &predicts[..len_n];
    let preds1 = &predicts[len_n..];
    let preds0_c = &predicts_chosen[..len_n];
    let preds1_c = &predicts_chosen[len_n..];

    //let mut dot = 0.0f64;
    let mut dot: usize = 0;
    for ni in 0..n {
        let g = convert_predicts_to_genotype(bit::bget(preds0, ni), bit::bget(preds1, ni));
        let g_c = convert_predicts_to_genotype(bit::bget(preds0_c, ni), bit::bget(preds1_c, ni));
        dot += (g as usize) * (g_c as usize);
        //dot += (g as f64) * (g_c as f64);
    }

    let r2 = ((dot as f64) / (n as f64) - gt_mean * gt_mean_chosen)
        / (gt_var.sqrt() * gt_var_chosen.sqrt());
    r2
}

pub fn calculate_r2s_nosimd<T>(
    r2s: &mut [f64],
    predictions: &GenotTwin,
    //predictions: &[B8],
    gt_means: &[f64],
    gt_vars: &[f64],
    is_pruned: &[bool],
    wgt_chosen: &WgtBoost,
    snv_indexs: &[T],
    //snv_indexs: &[SnvIndex],
    n: usize,
) where
    T: AsRef<SnvId> + std::marker::Sync, //Sync is required by rayon
{
    vec::fill(r2s, 0.0);

    let mi_chosen = wgt_chosen.wgt().kind().index_snv().unwrap();
    let chrom_chosen = wgt_chosen.wgt().kind().snv_index().chrom();

    r2s.par_iter_mut()
        .enumerate()
        .filter(|(mi, _)| !is_pruned[*mi] && snv_indexs[*mi].as_ref().chrom() == chrom_chosen)
        .for_each(|(mi, r2)| {
            *r2 = calculate_r2_nosimd(
                predictions.predictions_snv_s(mi),
                predictions.predictions_snv_s(mi_chosen),
                gt_means[mi],
                gt_vars[mi],
                gt_means[mi_chosen],
                gt_vars[mi_chosen],
                n,
            )
        });
}
 */

#[cfg(test)]
mod tests {
    use super::*;

    //use super::mistake;
    //use super::genotype;
    use genetics::{plink, SnvId};

    fn is_eq_f64(v: f64, w: f64, e: f64) -> bool {
        (v - w).abs() < e
    }

    /*
    fn setup_test() -> (Vec<B8>, usize, usize) {
        let fin = String::from("../../test/data/toy1/genot");
        let fin_snv = None;
        let fin_sample = None;

        let m_in: usize = plink::compute_num_snv(&fin).unwrap();
        println!("m_in: {}", m_in);
        let n_in: usize = plink::compute_num_sample(&fin).unwrap();
        println!("n_in: {}", n_in);
        // load snvs
        let snvs_in: Vec<Snv> = plink::load_snvs(&fin, m_in);
        let (m, use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (n, use_samples) = plink::make_use_samples(fin_sample, &fin, n_in);
        //let ys: Vec<B8> = plink::load_ys_b8(&fin, n, &use_samples);

        //(fin, ys, m, n, use_snvs, use_samples)

        let predictions = predict::generate_predictions(&fin, m, n, &use_snvs, &use_samples);

        (predictions, m, n)
    }

    fn setup_test_2() -> (Vec<B8>, Vec<B8>, usize, usize) {
        let (m, n) = (2, 4);
        //let mistakes: [B8; 2] = [0b0011_0111, 0b0000_0000];
        let x0: [u8; 4] = [0, 1, 1, 2];
        let x1: [u8; 4] = [2, 1, 1, 0];

        //let mistakes = mistake::generate_mistakes_from_x_snv(&x, &ys);
        let predict0 = predict::generate_predictions_from_x_snv(&x0);
        let predict1 = predict::generate_predictions_from_x_snv(&x1);
        println!("predict0: {:?}", predict0);

        (predict0, predict1, m, n)
    }

    fn setup_test_3() -> (Vec<B8>, Vec<B8>, usize, usize) {
        // larger test for n>8
        let (m, n) = (2, 9);
        //let mistakes: [B8; 2] = [0b0011_0111, 0b0000_0000];
        let x0: [u8; 9] = [0, 1, 1, 2, 2, 2, 1, 0, 2];
        let x1: [u8; 9] = [2, 0, 1, 0, 1, 2, 0, 2, 1];

        let predict0 = predict::generate_predictions_from_x_snv(&x0);
        let predict1 = predict::generate_predictions_from_x_snv(&x1);
        println!("predict0: {:?}", predict0);

        (predict0, predict1, m, n)
    }

    // test corner len=32 for simd
    fn setup_test_4() -> (Vec<B8>, Vec<B8>, usize, usize) {
        // larger test for n>8
        let (m, n) = (2, 32);
        //let mistakes: [B8; 2] = [0b0011_0111, 0b0000_0000];
        let mut x0 = Vec::with_capacity(32);
        // cannot be all the same values
        x0.push(1u8);
        for _ in 0..31 {
            x0.push(0u8);
        }
        let x1 = x0.clone();

        let predict0 = predict::generate_predictions_from_x_snv(&x0);
        let predict1 = predict::generate_predictions_from_x_snv(&x1);
        println!("predict0: {:?}", predict0);

        (predict0, predict1, m, n)
    }

    // when all vals are 0
    // r2 should be nan
    fn setup_test_5() -> (Vec<B8>, Vec<B8>, usize, usize) {
        // larger test for n>8
        let (m, n) = (2, 4);
        let x0: [u8; 4] = [0, 0, 0, 0];
        let x1: [u8; 4] = [2, 1, 1, 0];

        //let mistakes = mistake::generate_mistakes_from_x_snv(&x, &ys);
        let predict0 = predict::generate_predictions_from_x_snv(&x0);
        let predict1 = predict::generate_predictions_from_x_snv(&x1);
        println!("predict0: {:?}", predict0);

        (predict0, predict1, m, n)
    }

    // when all vals are the same 1
    // r2 should be nan
    fn setup_test_6() -> (Vec<B8>, Vec<B8>, usize, usize) {
        // larger test for n>8
        let (m, n) = (2, 4);
        let x0: [u8; 4] = [1, 1, 1, 1];
        let x1: [u8; 4] = [2, 1, 1, 0];

        //let mistakes = mistake::generate_mistakes_from_x_snv(&x, &ys);
        let predict0 = predict::generate_predictions_from_x_snv(&x0);
        let predict1 = predict::generate_predictions_from_x_snv(&x1);
        println!("predict0: {:?}", predict0);

        (predict0, predict1, m, n)
    }

    #[test]
    fn test_calculate_mean_m_2() {
        //let x0: [u8; 4] = [0, 1, 1, 2]; mean: 1, var: 0.5
        //let x1: [u8; 4] = [2, 1, 1, 0]; mean: 1, var: 0.5
        let mean0_exp = 1.0;
        let mean1_exp = 1.0;

        let (predict0, predict1, _, n) = setup_test_2();

        let mean0 = calculate_mean_m(&predict0, n);
        let mean1 = calculate_mean_m(&predict1, n);

        println!("mean0 {}", mean0);
        assert!(is_eq_f64(mean0, mean0_exp, 1e-7));
        assert!(is_eq_f64(mean1, mean1_exp, 1e-7));
    }
    #[test]
    fn test_calculate_var_m_2() {
        //let x0: [u8; 4] = [0, 1, 1, 2]; mean: 1, var: 0.5
        //let x1: [u8; 4] = [2, 1, 1, 0]; mean: 1, var: 0.5
        let var0_exp = 0.5;
        let var1_exp = 0.5;

        let (predict0, predict1, _, n) = setup_test_2();

        let mean0 = calculate_mean_m(&predict0, n);
        let var0 = calculate_var_m(mean0, &predict0, n);
        let mean1 = calculate_mean_m(&predict1, n);
        let var1 = calculate_var_m(mean1, &predict1, n);

        println!("var0 {}", var0);
        assert!(is_eq_f64(var0, var0_exp, 1e-7));
        assert!(is_eq_f64(var1, var1_exp, 1e-7));
    }

    #[test]
    fn test_calculate_r2_simd() {
        let (predictions, _, n) = setup_test();

        let predicts = predict::predictions_snv_s(&predictions, 0, n);
        let predicts_chosen = predict::predictions_snv_s(&predictions, 1, n);

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        unsafe {
            let gs_chosen: Vec<__m256i> = create_gs_epi8(predicts_chosen, n);

            let r2 = calculate_r2_simd(predicts, &gs_chosen, mean, var, mean_c, var_c, n);
            println!("r2 {}", r2);

            let r2_nosimd =
                calculate_r2_nosimd(predicts, predicts_chosen, mean, var, mean_c, var_c, n);

            assert!(is_eq_f64(r2, r2_nosimd, 1e-7));
        }
    }

    #[test]
    fn test_calculate_r2_nosimd() {
        let (predictions, _, n) = setup_test();

        let predicts = predict::predictions_snv_s(&predictions, 0, n);
        let predicts_chosen = predict::predictions_snv_s(&predictions, 1, n);

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        let r2 = calculate_r2_nosimd(predicts, predicts_chosen, mean, var, mean_c, var_c, n);
        println!("r2 {}", r2);
    }

    #[test]
    fn test_calculate_r2_self_simd() {
        let (predictions, _, n) = setup_test();

        let predicts = predict::predictions_snv_s(&predictions, 0, n);
        let predicts_chosen = predict::predictions_snv_s(&predictions, 0, n);

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        unsafe {
            let gs_chosen: Vec<__m256i> = create_gs_epi8(predicts_chosen, n);

            let r2 = calculate_r2_simd(predicts, &gs_chosen, mean, var, mean_c, var_c, n);
            println!("r2 {}", r2);

            assert!(is_eq_f64(r2, 1.0, 1e-7));
        }
    }

    #[test]
    fn test_calculate_r2_self_nosimd() {
        let (predictions, _, n) = setup_test();

        let predicts = predict::predictions_snv_s(&predictions, 0, n);
        let predicts_chosen = predict::predictions_snv_s(&predictions, 0, n);

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        let r2 = calculate_r2_nosimd(predicts, predicts_chosen, mean, var, mean_c, var_c, n);
        println!("r2 {}", r2);

        assert!(is_eq_f64(r2, 1.0, 1e-7));
    }

    #[test]
    fn test_calculate_r2_nosimd_2() {
        //let x0: [u8; 4] = [0, 1, 1, 2]; mean: 1, var: 0.5
        //let x1: [u8; 4] = [2, 1, 1, 0]; mean: 1, var: 0.5
        let r2_exp = -1.0;

        let (predict0, predict1, _, n) = setup_test_2();

        let predicts = &predict0;
        let predicts_chosen = &predict1;

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        let r2 = calculate_r2_nosimd(predicts, predicts_chosen, mean, var, mean_c, var_c, n);
        println!("r2 {}", r2);
        assert!(is_eq_f64(r2, r2_exp, 1e-7));
    }

    #[test]
    fn test_calculate_r2_simd_2() {
        let r2_exp = -1.0;

        let (predict0, predict1, _, n) = setup_test_2();

        let predicts = &predict0;
        let predicts_chosen = &predict1;

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        unsafe {
            let gs_chosen: Vec<__m256i> = create_gs_epi8(predicts_chosen, n);

            let r2 = calculate_r2_simd(predicts, &gs_chosen, mean, var, mean_c, var_c, n);
            println!("r2 {}", r2);

            assert!(is_eq_f64(r2, r2_exp, 1e-7));
        }
    }
    #[test]
    fn test_calculate_r2_simd_3() {
        let r2_exp = -0.3464;

        let (predict0, predict1, _, n) = setup_test_3();

        let predicts = &predict0;
        let predicts_chosen = &predict1;

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        unsafe {
            let gs_chosen: Vec<__m256i> = create_gs_epi8(predicts_chosen, n);

            let r2 = calculate_r2_simd(predicts, &gs_chosen, mean, var, mean_c, var_c, n);
            println!("r2 {}", r2);

            assert!(is_eq_f64(r2, r2_exp, 1e-4));
        }
    }

    #[test]
    fn test_calculate_r2_simd_4() {
        let (predict0, predict1, _, n) = setup_test_4();

        let predicts = &predict0;
        let predicts_chosen = &predict1;

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        unsafe {
            let gs_chosen: Vec<__m256i> = create_gs_epi8(predicts_chosen, n);

            let r2 = calculate_r2_simd(predicts, &gs_chosen, mean, var, mean_c, var_c, n);
            println!("r2 {}", r2);

            //assert!(is_eq_f64(r2, r2_exp, 1e-4));
        }
    }

    #[test]
    fn test_calculate_r2_nosimd_5() {
        //let r2_exp = -1.0;

        let (predict0, predict1, _, n) = setup_test_5();

        let predicts = &predict0;
        let predicts_chosen = &predict1;

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        let r2 = calculate_r2_nosimd(predicts, predicts_chosen, mean, var, mean_c, var_c, n);
        println!("r2 {}", r2);
        assert!(r2.is_nan());
    }

    #[test]
    fn test_calculate_r2_nosimd_6() {
        //let r2_exp = -1.0;

        let (predict0, predict1, _, n) = setup_test_6();

        let predicts = &predict0;
        let predicts_chosen = &predict1;

        let mean = calculate_mean_m(predicts, n);
        let var = calculate_var_m(mean, predicts, n);
        let mean_c = calculate_mean_m(predicts_chosen, n);
        let var_c = calculate_var_m(mean_c, predicts_chosen, n);

        let r2 = calculate_r2_nosimd(predicts, predicts_chosen, mean, var, mean_c, var_c, n);
        println!("r2 {}", r2);
        assert!(r2.is_nan());
    }
    */
}
