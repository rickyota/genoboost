use crate::WgtBoost;
use genetics::samples::prelude::*;
use genetics::wgt::Coef;
use genetics::{samples::CovsTrait, Samples};
use genetics::{CovId, Wgt};
use smartcore::linalg::naive::dense_matrix::DenseMatrix;
use smartcore::linalg::BaseMatrix;
use smartcore::linear::logistic_regression::LogisticRegression;

// TODO: move to genetic_rust
pub fn logistic_regression_covs(samples: &Samples) -> Vec<WgtBoost> {
    let covs_val = samples.covs();
    if let None = covs_val {
        return vec![];
    }

    let covs_val = covs_val.unwrap();

    let ys = samples.phe().inner_f64();
    // dense_matrix from 2d array
    // https://docs.rs/smartcore/latest/src/smartcore/linalg/naive/dense_matrix.rs.html#248
    let (v, row_n, col_n) = covs_val.vals_column_major_vec();
    let mt = DenseMatrix::from_vec(row_n, col_n, &v);
    //let mt = DenseMatrix::from_2d_vec(&covs_val.vals_column_major());
    // default: alpha=0.0
    // confirmed to output the same answer as sklearn
    let lr = LogisticRegression::fit(&mt, &ys, Default::default()).unwrap();
    println!("logreg");
    println!("intercept {}", lr.intercept());
    println!("coef {}", lr.coefficients());

    let intercept = lr.intercept().clone().to_row_vector()[0];

    let mut wgts_cov: Vec<WgtBoost> = Vec::new();

    let cov_id = CovId::new_const();
    let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));

    let iteration = 0;
    let wgt_boost_const = WgtBoost::construct_wgt_iteration(wgt_const, iteration);

    wgts_cov.push(wgt_boost_const);

    let coefs = lr.coefficients().clone().to_row_vector();
    assert_eq!(coefs.len(), covs_val.covs_n());
    let cov_indexs = covs_val.cov_indexs().unwrap();
    for wgt_i in 0..covs_val.covs_n() {
        let iteration = wgt_i + 1;
        let coef = coefs[wgt_i];
        let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        let wgt_cov = Wgt::construct_cov(cov_id, Coef::Linear(coef));
        let wgt_boost_cov = WgtBoost::construct_wgt_iteration(wgt_cov, iteration);
        wgts_cov.push(wgt_boost_cov);
    }

    wgts_cov
}
