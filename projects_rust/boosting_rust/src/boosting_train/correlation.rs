// move to genetics_rust

mod calculate;

use super::GenotTwin;
use super::WgtBoost;
use genetics::genot::BaseGenotTwin;
use genetics::{vec, SnvId};
use rayon::prelude::*;

//pub fn calculate_means_genotype(predictions: &[B8], m: usize, n: usize) -> Vec<f64> {
pub fn calculate_means_genotype(predictions: &GenotTwin, m: usize, n: usize) -> Vec<f64> {
    let mut means = vec![f64::NAN; m];

    means
        .iter_mut()
        .zip(predictions.iter_snv())
        .par_bridge()
        .for_each(|(mean, gref)| {
            *mean = calculate::calculate_mean_m(&gref);
        });

    // cannot call collect_into_vec()
    /*     predictions
    .iter()
    .par_bridge()
    .map(|pred| calculate::calculate_mean_m(&pred))
    .collect_into_vec(&mut means); */

    /*     means.par_iter_mut().enumerate().for_each(|(mi, mean)| {
        *mean = calculate::calculate_mean_m(&predictions.to_genot_twin_snv(mi));
        // *mean = calculate::calculate_mean_m(predictions.predictions_snv_s(mi), n);
        // *mean = calculate::calculate_mean_m(predict::predictions_snv_s(predictions, mi, n), n);
    }); */

    means
}

pub fn calculate_vars_genotype(
    means: &[f64],
    predictions: &GenotTwin,
    m: usize,
    n: usize,
) -> Vec<f64> {
    let mut vars = vec![f64::NAN; m];

    vars.iter_mut()
        .zip(means.iter())
        .zip(predictions.iter_snv())
        .par_bridge()
        .for_each(|((var, mean), gref)| {
            *var = calculate::calculate_var_m(*mean, &gref);
        });

    /*     vars.par_iter_mut()
    .zip(means.par_iter())
    .enumerate()
    .for_each(|(mi, (var, mean))| {
        *var = calculate::calculate_var_m(*mean, predictions.predictions_snv_s(mi), n);
    }); */

    vars
}

// solved: should rewrite &[SnvIndex] to &[&SnvIndex]
pub fn calculate_r2s<T>(
    r2s: &mut [f64],
    predictions: &GenotTwin,
    //predictions: &[B8],
    gt_means: &[f64],
    gt_vars: &[f64],
    is_pruned: &[bool],
    wgt_chosen: &WgtBoost,
    //wgt_chosen: &Wgt,
    snv_indexs: &[T],
    //snv_indexs: &[SnvIndex],
    //mi_self: usize,
    n: usize,
) where
    T: AsRef<SnvId> + std::marker::Sync,
{
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            log::debug!("Use SIMD.");
            unsafe {
                calculate::calculate_r2s_simd(
                    r2s,
                    predictions,
                    gt_means,
                    gt_vars,
                    is_pruned,
                    wgt_chosen,
                    snv_indexs,
                    n,
                );
            }
            return;
        }
    }

    log::debug!("Do not use SIMD.");
    calculate::calculate_r2s_nosimd(
        r2s,
        predictions,
        gt_means,
        gt_vars,
        is_pruned,
        wgt_chosen,
        snv_indexs,
        n,
    );
}

// rayon
/// when chosen has all the same values, then just make the chosen true
/// when some snv have all the same values, then the snv is false
pub fn update_is_pruned(is_pruned: &mut [bool], r2s: &[f64], clump_r2: f64, mi_chosen: usize) {
    // self
    is_pruned[mi_chosen] = true;
    is_pruned
        .par_iter_mut()
        .zip(r2s.par_iter())
        .filter(|(_, &r2)| (r2 > clump_r2) || (r2 < -clump_r2))
        .for_each(|(v, _)| *v = true);
}

pub fn check_all_pruned(is_pruned: &[bool]) -> bool {
    is_pruned.iter().all(|&v| v)
}

/*
#[cfg(test)]
mod tests {
    use super::*;
gg
    #[test]
    fn test_search_min_loss_gt() {}
}
*/
