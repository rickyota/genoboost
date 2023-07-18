mod adjust_coef;
mod calc;

use super::epsilon;
use super::sample_weight::SampleWeight;
use super::table;
use super::{table::CONTINGENCY_TABLE_FILL, BoostType};
use crate::{ContingencyTable, EffEps, Eps};
use genetics::genot::prelude::*;
use genetics::genot::GenotSnvRef;
use genetics::samples::prelude::*;
use genetics::wgt::Coef;

// TODO: rename
pub fn calculate_coef_root_ada(
    pred_s: &[u8],
    sample_weight: &SampleWeight,
    //ps: &[f64],
    phe: &Phe,
    learning_rate: f64,
    eps: Option<Eps>,
    boost_type: BoostType,
) -> (Coef, bool) {
    let (table_sum, is_eps) =
        table::calculate_table_eps(&pred_s, sample_weight.ps().unwrap(), phe, eps, boost_type);
    let coef_ti = calculate_coefficients_ada(table_sum, boost_type, learning_rate, eps);
    (coef_ti, is_eps)
}

// assume table after eps
pub fn calculate_coefficients_ada(
    table: ContingencyTable,
    boost_type: BoostType,
    lr: f64,
    eps: Option<Eps>,
) -> Coef {
    match boost_type {
        BoostType::Ada => {
            let table2_sum = table.two();
            let (d, n) = table2_sum;
            let alpha = (d / n).ln() / 2.0;

            let alpha = lr * alpha;
            Coef::Binary((0.0, alpha))
        }
        BoostType::ConstAda => {
            let table4_sum = table.four();
            let (d1, n1, d0, n0) = table4_sum;
            let const_ti = ((d1 * d0) / (n1 * n0)).ln() / 4.0;
            let alpha_ti = ((d1 * n0) / (n1 * d0)).ln() / 4.0;

            let const_ti = lr * const_ti;
            let alpha_ti = lr * alpha_ti;
            Coef::Binary((const_ti, alpha_ti))
        }
        BoostType::FreeModelMissing => {
            let table7_sum = table.seven();
            let (d2, n2, d1, n1, d0, n0, _) = table7_sum;

            // TODO: clean do not want to use eps here
            if eps.is_some() && eps.unwrap().dom() {
                // TODO: create table:judge_eps(table)
                let s0;
                let s1;
                let s2;
                if (d2 < CONTINGENCY_TABLE_FILL) || (n2 < CONTINGENCY_TABLE_FILL) {
                    let d1_new = d2 + d1;
                    let n1_new = n2 + n1;
                    s1 = (d1_new / n1_new).ln() / 2.0;
                    s2 = s1;
                    s0 = (d0 / n0).ln() / 2.0;
                } else if (d0 < CONTINGENCY_TABLE_FILL) || (n0 < CONTINGENCY_TABLE_FILL) {
                    let d1_new = d0 + d1;
                    let n1_new = n0 + n1;
                    s1 = (d1_new / n1_new).ln() / 2.0;
                    s0 = s1;
                    s2 = (d2 / n2).ln() / 2.0;
                } else {
                    s0 = (d0 / n0).ln() / 2.0;
                    s1 = (d1 / n1).ln() / 2.0;
                    s2 = (d2 / n2).ln() / 2.0;
                }

                let s0 = lr * s0;
                let s1 = lr * s1;
                let s2 = lr * s2;
                Coef::Score4((s0, s1, s2, 0.0))
            } else {
                let s0 = (d0 / n0).ln() / 2.0;
                let s1 = (d1 / n1).ln() / 2.0;
                let s2 = (d2 / n2).ln() / 2.0;

                let s0 = lr * s0;
                let s1 = lr * s1;
                let s2 = lr * s2;
                Coef::Score4((s0, s1, s2, 0.0))
            }
        }
        _ => panic!(),
    }
}

// DO NOT use predict here: troublesome
// now common fn can be used for loss and wgt
// TODO: integrate to calcualte_coef_root_ada
pub fn calculate_coef_root_logit(
    gsnv: &GenotSnvRef,
    sample_weight: &SampleWeight,
    //wzs_pad: &[f64],
    //wls_pad: &[f64],
    phe: &Phe,
    learning_rate: f64,
    eps: Option<Eps>,
    eff_eps: Option<EffEps>,
    boost_type: BoostType,
) -> (Coef, bool, bool) {
    unsafe {
        match boost_type {
            BoostType::Logit | BoostType::LogitNoMissing => {
                let epsilons_wzs =
                    epsilon::calculate_epsilons_logit_wzs(sample_weight.wzs().unwrap(), phe, eps);
                let epsilons_wls =
                    epsilon::calculate_epsilons_logit_wls(sample_weight.wls().unwrap(), phe, eps);

                let (coef_ti, is_eps, is_eff_eps) = calculate_coef_logit_eps(
                    gsnv,
                    sample_weight.wzs_pad().unwrap(),
                    sample_weight.wls_pad().unwrap(),
                    //wzs_pad,
                    //wls_pad,
                    phe,
                    epsilons_wzs,
                    epsilons_wls,
                    eps,
                    boost_type,
                    eff_eps,
                    false,
                    true,
                );
                let coef_ti = coef_lr(coef_ti, learning_rate, boost_type);

                (coef_ti, is_eps, is_eff_eps)
            }
            BoostType::LogitAdd => {
                let coef_ti = calculate_coef_logit_add(
                    gsnv,
                    //&genot.to_genot_snv(mi),
                    sample_weight.wzs_pad().unwrap(),
                    sample_weight.wls_pad().unwrap(),
                    //wzs_pad,
                    //wls_pad,
                    //phe,
                    learning_rate,
                );
                (coef_ti, false, false)
            }
            _ => panic!(),
        }
    }
}

pub unsafe fn calculate_coef_logit_eps(
    gsnv: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
    phe: &Phe,
    epsilons_wzs: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    epsilons_wls: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    eps: Option<Eps>,
    //learning_rate: Option<f64>,
    boost_type: BoostType,
    eff_eps: Option<EffEps>,
    on_loss: bool,
    verbose: bool,
) -> (Coef, bool, bool) {
    //use std::time::Instant;
    //let start_time = Instant::now();

    let (wzs_sum, wls_sum) = calc::calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
    //calculate_coef_gt_logit_simd_sm(&genot.to_genot_snv(mi), wzs, wls, phe);

    //println!("afr wzs_sum: {} sec",  start_time.elapsed().as_micros());

    //log::debug!("wzs_sum bfr {:?}",wzs_sum);
    //log::debug!("wls_sum bfr {:?}",wls_sum);

    let table8_count = gsnv.stat_contingency_table(phe);
    //let table8_count = predicts.stat_contingency_table(phe);
    //let table8_count = genot.to_genot_snv(mi).stat_contingency_table(phe);

    //println!("afr table8: {} sec",  start_time.elapsed().as_micros());

    let (wzs_sum, wls_sum, is_eps) = table::adjust_eps_logit(
        wzs_sum,
        wls_sum,
        epsilons_wzs,
        epsilons_wls,
        eps,
        table8_count,
    );

    //println!("afr adj: {} sec",  start_time.elapsed().as_micros());

    //log::debug!("wzs_sum afr {:?}",wzs_sum);
    //log::debug!("wls_sum afr {:?}",wls_sum);

    let (coef, is_eff_eps) = match boost_type {
        BoostType::Logit => {
            calculate_coef_logit(wzs_sum, wls_sum, eff_eps, table8_count, on_loss, verbose)
        }
        BoostType::LogitNoMissing => {
            calculate_coef_logit_no_missing(
                wzs_sum,
                wls_sum,
                eff_eps,
                table8_count,
                on_loss,
                verbose,
            )
            //if on_loss{
            //    coefficient::calculate_coef_on_loss_logit_no_missing(wzs_sum, wls_sum, eff_eps, table8_count,verbose)
            //}else{
            //    coefficient::calculate_coef_logit_no_missing(wzs_sum, wls_sum, eff_eps, table8_count,verbose)
            //}
        }
        _ => panic!("wrong"),
    };

    //let coef: Coef = coefficient::calculate_coef_logit(wzs_sum, wls_sum);

    //println!("afr coef: {} sec",  start_time.elapsed().as_micros());

    // CANNOT USE learnin_rate here since lr should not be used when calculating loss function

    /*
    let lr = learning_rate.unwrap_or(1.0);

    let coef: Coef = match boost_type {
        BoostType::Logit => {
            let (s0, s1, s2, _sm) = coef.score4_f64();
            Coef::Score4((lr * s0, lr * s1, lr * s2, 0.0))
        }
        BoostType::LogitNoMissing => {
            let (s0, s1, s2) = coef.score3_f64();
            Coef::Score3((lr * s0, lr * s1, lr * s2))
        }
        _ => panic!("wrong"),
    };
    */

    (coef, is_eps, is_eff_eps)
}

/// deprecated: lr
/* pub unsafe fn calculate_coef_logit_add_lr(
    gsnv: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
    //_phe: &Phe,
    //epsilons_wzs: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //epsilons_wls: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //eps: Eps,
    lr: f64,
) -> Coef {
    //use std::time::Instant;
    //let start_time = Instant::now();

    let (wzs_sum, wls_sum) = calc::calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
    //calculate_coef_gt_logit_simd_sm(&genot.to_genot_snv(mi), wzs, wls, phe);

    //println!("afr wzs_sum: {} sec",  start_time.elapsed().as_micros());

    //log::debug!("wzs_sum bfr {:?}",wzs_sum);
    //log::debug!("wls_sum bfr {:?}",wls_sum);

    //let table8_count = predicts.stat_contingency_table_simd(phe);
    //let table8_count = predicts.stat_contingency_table(phe);
    //let table8_count = genot.to_genot_snv(mi).stat_contingency_table(phe);

    //println!("afr table8: {} sec",  start_time.elapsed().as_micros());

    //let (wzs_sum, wls_sum, is_eps) = table::adjust_eps_logit(
    //    wzs_sum,
    //    wls_sum,
    //    epsilons_wzs,
    //    epsilons_wls,
    //    eps,
    //    table8_count,
    //);

    //println!("afr adj: {} sec",  start_time.elapsed().as_micros());

    //log::debug!("wzs_sum afr {:?}",wzs_sum);
    //log::debug!("wls_sum afr {:?}",wls_sum);

    let coef: Coef = calc_coef_logit_add(wzs_sum, wls_sum);

    //println!("afr coef: {} sec",  start_time.elapsed().as_micros());

    // makd func
    //let lr = lr.unwrap_or(1.0);

    let (c, a) = coef.linearconst_f64();

    let coef = Coef::LinearConst((lr * c, lr * a));

    coef
}
 */

pub unsafe fn calculate_coef_logit_add(
    gsnv: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
    //_phe: &Phe,
    //epsilons_wzs: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //epsilons_wls: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //eps: Eps,
    lr: f64,
) -> Coef {
    /*     //use std::time::Instant;
    //let start_time = Instant::now();

    let (wzs_sum, wls_sum) = calc::calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
    //calculate_coef_gt_logit_simd_sm(&genot.to_genot_snv(mi), wzs, wls, phe);

    let coef: Coef = calc_coef_logit_add(wzs_sum, wls_sum);

    //println!("afr coef: {} sec",  start_time.elapsed().as_micros());

    let (c, a) = coef.linearconst_f64();

    let coef = Coef::LinearConst((c, a));
    //let coef = Coef::LinearConst((lr * c, lr * a)); */

    let coef = calculate_coef_logit_add_on_loss(gsnv, wzs_pad, wls_pad);

    coef.apply_lr(lr)
}

// for coef on calculating loss
pub unsafe fn calculate_coef_logit_add_on_loss(
    gsnv: &GenotSnvRef,
    wzs_pad: &[f64],
    wls_pad: &[f64],
    //_phe: &Phe,
    //epsilons_wzs: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //epsilons_wls: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    //eps: Eps,
    //lr: f64,
) -> Coef {
    //use std::time::Instant;
    //let start_time = Instant::now();

    let (wzs_sum, wls_sum) = calc::calculate_coef_gt_logit_simd_sm(gsnv, wzs_pad, wls_pad);
    //calculate_coef_gt_logit_simd_sm(&genot.to_genot_snv(mi), wzs, wls, phe);

    let coef: Coef = calc_coef_logit_add(wzs_sum, wls_sum);

    //println!("afr coef: {} sec",  start_time.elapsed().as_micros());

    let (c, a) = coef.linearconst_f64();

    let coef = Coef::LinearConst((c, a));
    //let coef = Coef::LinearConst((lr * c, lr * a));

    coef
}

pub fn calculate_coef_logit(
    wzs_sum: (f64, f64, f64),
    wls_sum: (f64, f64, f64),
    eff_eps: Option<EffEps>,
    table8_count: (usize, usize, usize, usize, usize, usize, usize, usize),
    on_loss: bool,
    verbose: bool,
) -> (Coef, bool) {
    let sm = 0.0f64;

    if let (Coef::Score3((s0, s1, s2)), is_eff_eps) =
        calculate_coef_logit_no_missing(wzs_sum, wls_sum, eff_eps, table8_count, on_loss, verbose)
    {
        (Coef::Score4((s0, s1, s2, sm)), is_eff_eps)
    } else {
        panic!()
    }
}

pub fn calculate_coef_logit_no_missing(
    wzs_sum: (f64, f64, f64),
    wls_sum: (f64, f64, f64),
    eff_eps: Option<EffEps>,
    // for eff_eps
    table8_count: (usize, usize, usize, usize, usize, usize, usize, usize),
    on_loss: bool,
    verbose: bool,
) -> (Coef, bool) {
    let (s0, s1, s2) = calculate_coef_from_weights(wzs_sum, wls_sum);

    if verbose {
        log::debug!("Coef2,1,0 bfr EffEps {}, {}, {}", s2, s1, s0);
    }

    if !on_loss {
        adjust_coef::adjust_eff_logit_no_missing((s0, s1, s2), table8_count, eff_eps, verbose)
    } else {
        // calculate coef for loss
        // some eps_eff does not apply on coef for loss
        if eff_eps.is_some() && (!eff_eps.unwrap().is_on_update()) {
            if verbose {
                log::debug!("Adjust eff_eps since for loss.");
            }
            adjust_coef::adjust_eff_logit_no_missing((s0, s1, s2), table8_count, eff_eps, verbose)
        } else {
            // legacy
            if verbose {
                log::debug!("Do not adjust eff_eps since for loss.");
            }
            (Coef::Score3((s0, s1, s2)), false)
        }
    }
}

pub fn calculate_coef_from_weights(
    wzs_sum: (f64, f64, f64),
    wls_sum: (f64, f64, f64),
) -> (f64, f64, f64) {
    // TODO: return score

    let (wzs_sum2, wzs_sum1, wzs_sum0) = wzs_sum;
    let (wls_sum2, wls_sum1, wls_sum0) = wls_sum;

    /*     let s2 = wzs_sum2 / wls_sum2;
    let s1 = wzs_sum1 / wls_sum1;
    let s0 = wzs_sum0 / wls_sum0; */

    let s0 = if wls_sum0 == 0.0f64 {
        0.0f64
    } else {
        wzs_sum0 / wls_sum0
    };

    let s1 = if wls_sum1 == 0.0f64 {
        0.0f64
    } else {
        wzs_sum1 / wls_sum1
    };

    let s2 = if wls_sum2 == 0.0f64 {
        s1
    } else {
        wzs_sum2 / wls_sum2
    };

    (s0, s1, s2)
}

pub fn calc_coef_logit_add(wzs_sum: (f64, f64, f64), wls_sum: (f64, f64, f64)) -> Coef {
    let (u2, u1, u0) = wzs_sum;
    let (w2, w1, w0) = wls_sum;

    let denom = w0 * w1 + w1 * w2 + 4.0 * w2 * w0;

    let c = ((w1 + 4.0 * w2) * u0 + 2.0 * w2 * u1 - w1 * u2) / denom;
    let a = ((-w1 - 2.0 * w2) * u0 + (-w2 + w0) * u1 + (2.0 * w0 + w1) * u2) / denom;

    //when w0=w2=0, denom=0.0
    //if denom==0.0{
    //    log::debug!("denom=0.0");
    //    log::debug!("a {}",a);
    //    log::debug!("c {}",c);
    //    log::debug!("w0 {}",w0);
    //    log::debug!("w1 {}",w1);
    //    log::debug!("w2 {}",w2);
    //}

    /*
    let v0 = 4.0 * w2 + w1;
    let v1 = 2.0 * w2 + w1;
    let v2 = w2 + w1 + w0;
    let v3 = -2.0 * wz2 - wz1;
    let v4 = -(wz2 + wz1 + wz0);

    let c = (v1 * v3 - v0 * v4) / (v0 * v2 - v1 * v1);
    let a = (v1 * v4 - v2 * v3) / (v0 * v2 - v1 * v1);
     */

    Coef::LinearConst((c, a))
}

// boost_type is necessary since how to deal with sm could be different
// TODO: common for not only Logit
pub fn coef_lr(coef: Coef, lr: f64, boost_type: BoostType) -> Coef {
    //let lr = learning_rate.unwrap_or(1.0);

    match boost_type {
        BoostType::Logit => {
            let (s0, s1, s2, _sm) = coef.score4_f64();
            Coef::Score4((lr * s0, lr * s1, lr * s2, 0.0))
        }
        BoostType::LogitNoMissing => {
            let (s0, s1, s2) = coef.score3_f64();
            Coef::Score3((lr * s0, lr * s1, lr * s2))
        }
        _ => panic!("wrong"),
    }
}

pub fn calculate_coef_logit_const(wzs: &[f64], wls: &[f64]) -> (Coef, bool) {
    let wzs_all: f64 = wzs.iter().sum();
    let wls_all: f64 = wls.iter().sum();
    //let wzs_all: f64 = wzs[..n].iter().sum();
    //let wls_all: f64 = wls[..n].iter().sum();

    let alpha = wzs_all / wls_all;

    (Coef::Single(alpha), false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use genetics::alloc;

    //fn is_eq_f64(v: f64, w: f64, e: f64) -> bool {
    //    (v - w).abs() < e
    //}

    #[test]
    fn test_calculate_coefficients_freemodelmissing() {
        let t = ContingencyTable::new_seven((0.02, 0.01, 0.1, 0.2, 0.3, 0.3, 0.07));
        let coef = calculate_coefficients_ada(t, BoostType::FreeModelMissing, 1.0, Some(Eps::Med));
        assert_eq!(
            coef,
            Coef::Score4((0.0, 0.5f64.ln() / 2.0, 2.0f64.ln() / 2.0, 0.0))
        );
    }

    #[test]
    fn test_calculate_coefficients_freemodelmissing_lr() {
        let t = ContingencyTable::new_seven((0.02, 0.01, 0.1, 0.2, 0.3, 0.3, 0.07));
        let coef = calculate_coefficients_ada(t, BoostType::FreeModelMissing, 0.1, Some(Eps::Med));
        assert_eq!(
            coef,
            Coef::Score4((0.0, 0.5f64.ln() / 20.0, 2.0f64.ln() / 20.0, 0.0))
        );
    }
    /// for speed
    #[allow(dead_code)]
    fn setup_test4_logit() -> (Genot, Phe, Vec<f64>, Vec<f64>) {
        let n = 1_000_000;
        let x = vec![0u8; n];
        let genot = Genot::new(1, n, &x);
        let y = vec![true; n];
        let phe = Phe::new(&y);

        let mut wzs: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        wzs.resize(n + 32, 0.0f64);
        let mut wls: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        wls.resize(n + 32, 0.0f64);

        (genot, phe, wzs, wls)
    }

    /*
    // too measure speed
    #[test]
    fn test_calculate_coef_logit_eps_tmp() {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                let (genot, phe, wzs, wls) = setup_test4_logit();

                unsafe {
                    calculate_coef_logit_eps(
                        &genot.to_genot_snv(0),
                        &wzs,
                        &wls,
                        &phe,
                        (0.1, 0.1),
                        (0.1, 0.1),
                        Eps::MedLarge2,
                        //None,
                        BoostType::Logit,
                        None,
                        false,
                        true,
                    );
                }
            }
        }
    }
    */
}
