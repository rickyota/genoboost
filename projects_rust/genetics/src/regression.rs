use crate::samples::prelude::*;
use crate::{wgt::Coef, Wgt};
use crate::{CovId, Covs, CovsTrait, Samples};

//mod logistic;
//mod smartcore;
mod linfa;
#[cfg(feature = "pyo3")]
mod sklearn_py;

use crate::vec;

// linfa
// [ref](https://github.com/zupzup/rust-ml-example/blob/main/src/main.rs)
// [linfa_logistic ref](https://crates.io/crates/linfa-logistic)
// https://docs.rs/linfa-logistic/0.7.0/linfa_logistic/type.LogisticRegression.html
pub fn logistic_regression_covs(samples: &Samples) -> Vec<Wgt> {
    if let Some(covs) = samples.covs() {
        let (coefs, intercept) = linfa::logistic_regression_covs_samples(samples);
        //let (coefs, intercept) = smartcore::logistic_regression_covs_samples(samples);

        assert_eq!(coefs.len(), covs.covs_n());

        wgts_cov_from_coef(&coefs, intercept, covs)

        //let mut wgts_cov: Vec<Wgt> = Vec::new();

        //let cov_id = CovId::new_const();
        //let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
        //wgts_cov.push(wgt_const);

        ////let coefs = lr.params().iter().map(|x| *x).collect::<Vec<f64>>();
        //////let coefs = lr.coefficients().clone().to_row_vector();
        ////assert_eq!(coefs.len(), covs_val.covs_n());
        //let cov_indexs = covs.cov_indexs().unwrap();
        //for wgt_i in 0..covs.covs_n() {
        //    let coef = coefs[wgt_i];
        //    if coef.is_nan() {
        //        panic!("coef is NaN.");
        //    }
        //    let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        //    let wgt_cov = Wgt::new_cov(cov_id, Coef::Linear(coef));
        //    wgts_cov.push(wgt_cov);
        //}

        //wgts_cov
    } else {
        // regression on const

        log::debug!("Regression on const value only.");

        log::debug!("ys true count {}", samples.phe_unwrap().count());
        log::debug!("ys false count {}", samples.phe_unwrap().count_false());

        let intercept = logistic_regression_const_samples(samples);
        wgts_cov_from_intercept(intercept)

        //let y1_n = samples.phe_unwrap().count();
        //let y0_n = samples.phe_unwrap().count_false();

        //let intercept = logistic_regression_const(y0_n, y1_n);

        //let mut wgts_cov: Vec<Wgt> = Vec::new();
        //if intercept.is_nan() {
        //    panic!("Intercept is NaN.");
        //}
        //let cov_id = CovId::new_const();
        //let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
        //wgts_cov.push(wgt_const);

        //wgts_cov
    }
}

fn logistic_regression_const_samples(samples: &Samples) -> f64 {
    let y1_n = samples.phe_unwrap().count();
    let y0_n = samples.phe_unwrap().count_false();

    let intercept = logistic_regression_const(y0_n, y1_n);

    if intercept.is_nan() {
        panic!("Intercept is NaN.");
    }
    intercept
}

fn logistic_regression_const(y0_n: usize, y1_n: usize) -> f64 {
    // regression on const
    let y1_n = y1_n as f64;
    let y0_n = y0_n as f64;

    // TODO: test
    let intercept = (y1_n / y0_n).ln();
    intercept
}

fn wgts_cov_from_coef(coefs: &[f64], intercept: f64, covs: &Covs) -> Vec<Wgt> {
    let mut wgts_cov: Vec<Wgt> = Vec::new();

    let cov_id = CovId::new_const();
    let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
    wgts_cov.push(wgt_const);

    let cov_indexs = covs.cov_indexs().unwrap();
    for wgt_i in 0..covs.covs_n() {
        let coef = coefs[wgt_i];
        // already checked
        //if coef.is_nan() {
        //    panic!("coef is NaN.");
        //}
        let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        let wgt_cov = Wgt::new_cov(cov_id, Coef::Linear(coef));
        wgts_cov.push(wgt_cov);
    }

    wgts_cov
}

fn wgts_cov_from_intercept(intercept: f64) -> Vec<Wgt> {
    let cov_id = CovId::new_const();
    let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
    vec![wgt_const]
}

#[cfg(not(feature = "pyo3"))]
pub fn logistic_regression_covs_sampleweights(_: &Samples, _: &[f64]) -> Vec<Wgt> {
    panic!("Use feature=pyo3");
}

#[cfg(feature = "pyo3")]
pub fn logistic_regression_covs_sampleweights(samples: &Samples, ps: &[f64]) -> Vec<Wgt> {
    if let Some(covs) = samples.covs() {
        let (coefs, intercept) =
            sklearn_py::logistic_regression_covs_sample_sampleweights(samples, ps);

        assert_eq!(coefs.len(), covs.covs_n());

        wgts_cov_from_coef(&coefs, intercept, covs)
    } else {
        // regression on const

        unimplemented!("Not implemented yet.")

        //log::debug!("Regression on const value only.");

        //log::debug!("ys true count {}", samples.phe_unwrap().count());
        //log::debug!("ys false count {}", samples.phe_unwrap().count_false());

        //let y1_n = samples.phe_unwrap().count();
        //let y0_n = samples.phe_unwrap().count_false();

        //let intercept = logistic_regression_const(y0_n, y1_n);

        //let mut wgts_cov: Vec<Wgt> = Vec::new();
        //if intercept.is_nan() {
        //    panic!("Intercept is NaN.");
        //}
        //let cov_id = CovId::new_const();
        //let wgt_const = Wgt::new_cov(cov_id, Coef::Linear(intercept));
        //wgts_cov.push(wgt_const);

        //wgts_cov
    }
}

// for nagelkerke
pub fn logistic_regression_vec2(
    v: Vec<Vec<f64>>, // col major; unnormed 2d
    ys: Vec<bool>,
) -> (Vec<f64>, f64) {
    let means = vec::mean_v2(&v);
    let stds = vec::std_v2(&v);

    // convert to row major
    let v_norm = vec::convert_vec2d_to_row_major(&vec::norm_vec2d(v));
    let (v_norm_1d, (row_n, col_n)) = vec::convert_vec2d_to_vec1(v_norm);
    //let (v_norm_1d, (col_n, row_n)) = vec::convert_vec2d_to_vec1(v_norm);

    assert_eq!(col_n, 1);

    log::debug!("means {:?}", means);
    log::debug!("stds {:?}", stds);

    log::debug!("v_norm_1d[0],[1] {:?}, {:?}", v_norm_1d[0], v_norm_1d[1]);

    linfa::logistic_regression_vec_linfa_normalize(v_norm_1d, col_n, row_n, ys, &means, &stds)
}

// for nagelkerke
// 1d only for now
pub fn logistic_regression_vec_no_intercept_1col(
    //v: Vec<Vec<f64>>, // col major; unnormed 2d
    v: Vec<f64>, // col major; unnormed 2d
    ys: Vec<bool>,
) -> f64 {
    //let means = vec::mean_v2(&v);
    //let stds = vec::std_v2(&v);

    // convert to row major
    //let v_norm = vec::convert_vec2d_to_row_major(&vec::norm_vec2d(v));
    //let (v_norm_1d, (row_n, col_n)) = vec::convert_vec2d_to_vec1(v_norm);
    //let (v_norm_1d, (col_n, row_n)) = vec::convert_vec2d_to_vec1(v_norm);

    //assert_eq!(col_n, 1);

    //log::debug!("means {:?}", means);
    //log::debug!("stds {:?}", stds);

    //log::debug!("v_norm_1d[0],[1] {:?}, {:?}", v_norm_1d[0], v_norm_1d[1]);

    linfa::logistic_regression_vec_linfa_no_intercept_1col(v, ys)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_logisticregression_const() {
        let y1_n = 100;
        let y0_n = 10;
        let intercept = logistic_regression_const(y0_n, y1_n);
        assert_float_absolute_eq!(intercept, (100.0f64 / 10.0).ln());
    }

    #[test]
    fn test_logistic_regression_vec2() {
        // col major
        let v = vec![vec![-1.0f64, -0.01, 0.01, 1.0]];
        let y = vec![false, true, false, true];

        //let v_norm = vec::convert_vec2d_to_row_major(&vec::norm_vec2d(v));
        //let (v_norm_1d, row_n, col_n) = vec::convert_vec2d_to_vec1_row_major(v_norm);
        let (coefs, intercept) = logistic_regression_vec2(v, y);
        //logistic_regression_covs_vec_linfa(v_norm_1d, row_n, col_n, y, &means, &stds);
        assert_eq!(coefs.len(), 1);
        assert_float_absolute_eq!(intercept, 0.0);

        // expected answer was calculated using https://stats.blue/Stats_Suite/logistic_regression_calculator.html
        assert_float_absolute_eq!(coefs[0], 5.2672, 1e-3);
    }
}
