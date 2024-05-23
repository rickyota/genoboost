use crate::WgtBoost;
use genetics::Samples;
//use genetics::{CovId, Wgt};
use genetics::regression;
//use smartcore::linalg::naive::dense_matrix::DenseMatrix;
//use smartcore::linalg::BaseMatrix;
//use smartcore::linear::logistic_regression::LogisticRegression;

//use numpy::PyArray;
//use pyo3::{prelude::*};

/* // moved to genetic_rust
fn logistic_regression_sklearn(x: &Vec<Vec<f64>>, y: &[f64], w: &[f64]) -> Vec<f64> {
    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/src_py/logreg.py"));

    //pyo3::append_to_inittab!();
    // necessary without auto-initialize for static link
    pyo3::prepare_freethreaded_python();

    //Python::with_gil(|py| -> PyArray1<f64> {
    Python::with_gil(|py| -> Vec<f64> {
        PyModule::from_code(py, py_app, "logreg", "logreg").unwrap();
        let app: Py<PyAny> = PyModule::from_code(py, py_app, "", "")
            .unwrap()
            .getattr("logreg_sklearn")
            .unwrap()
            .into();
        //let res = app.call0(py);
        // [ref](https://docs.rs/numpy/latest/numpy/array/struct.PyArray.html#method.from_vec2)
        let args = (
            PyArray::from_vec2(py, x).unwrap(),
            PyArray::from_slice(py, y),
            PyArray::from_slice(py, w),
        );
        //let args=(PyArray::from_vec2(py, x).unwrap(),PyArray::from_vec(py,y),Pyarray::from_vec(py,w));
        //let args=(PyArray::from_vec2(py, &x).unwrap(), y.into_pyarray(py),w.into_pyarray(py));
        let res = app.call1(py, args);
        //let res = app.call1(py,(v ,));

        //let ans=res.unwrap().as_ref(py);
        //let ans: Array1<f64> =res.unwrap().as_ref(py).extract().unwrap();
        let ans: Vec<f64> = res.unwrap().as_ref(py).extract().unwrap();
        //let ans: PyArray1<f64> =res.unwrap().as_ref(py);

        ans
    })
}
 */

/*
// moved to genetic_rust
// deprecated to avoid pyo3
pub fn logistic_regression_covs_sampleweight(
    samples: &Samples,
    ps: &[f64],
    iteration_start: usize,
) -> Vec<WgtBoost> {
    let covs_val = samples.covs();
    if let None = covs_val {
        return vec![];
    }

    let covs_val = covs_val.unwrap();

    let ys = samples.phe().inner_f64();

    let v = covs_val.vals_column_major();

    //let ps_vec=ps.to_vec();

    let intercept_coef = logistic_regression_sklearn(&v, &ys, &ps[..ys.len()]);

    let intercept = intercept_coef[0];
    let coefs = &intercept_coef[1..];

    /*
    // dense_matrix from 2d array
    // cannot set weight
    // https://docs.rs/smartcore/latest/src/smartcore/linalg/naive/dense_matrix.rs.html#248
    let (v, row_n, col_n) = covs_val.vals_column_major_vec();
    let mt = DenseMatrix::from_vec(row_n, col_n, &v);
    //let mt = DenseMatrix::from_2d_vec(&covs_val.vals_column_major());
    // default: alpha=0.0
    // confirmed to output the same answer as sklearn
    let lr = LogisticRegression::fit(&mt, &ys, Default::default()).unwrap();
    log::debug!("logreg");
    log::debug!("intercept {}", lr.intercept());
    log::debug!("coef {}", lr.coefficients());
    let intercept = lr.intercept().clone().to_row_vector()[0];
    */

    let mut wgts_cov: Vec<WgtBoost> = Vec::new();

    let cov_id = CovId::new_const();
    let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));

    //let iteration = 0;
    let iteration = iteration_start;
    let wgt_boost_const = WgtBoost::construct_wgt_iteration(wgt_const, iteration);

    wgts_cov.push(wgt_boost_const);

    //let coefs = lr.coefficients().clone().to_row_vector();
    assert_eq!(coefs.len(), covs_val.covs_n());
    let cov_indexs = covs_val.cov_indexs().unwrap();
    // FIXME: iteration
    //for wgt_i in 0..covs_val.covs_n() {
    for wgt_i in 0..covs_val.covs_n() {
        let iteration = iteration_start + wgt_i + 1;
        let coef = coefs[wgt_i];
        let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        let wgt_cov = Wgt::construct_cov(cov_id, Coef::Linear(coef));
        let wgt_boost_cov = WgtBoost::construct_wgt_iteration(wgt_cov, iteration);
        wgts_cov.push(wgt_boost_cov);
    }

    wgts_cov
}
 */

/* 
// TODO: should pass WgtBoosts and add inside
pub fn logistic_regression_covs_sampleweight(
    samples: &Samples,
    ps: &[f64],
    iteration_start: usize,
) -> Vec<WgtBoost> {
    log::debug!("use sklearn for logreg_cov.");
    let wgts = regression::logistic_regression_covs_sampleweights_sklearn(samples, ps);

    let mut wgts_cov: Vec<WgtBoost> = Vec::new();
    for (wgt_i, wgt) in wgts.into_iter().enumerate() {
        let iteration = iteration_start + wgt_i;
        let wgt_boost_cov = WgtBoost::construct_wgt_iteration(wgt, iteration);
        wgts_cov.push(wgt_boost_cov);
    }

    wgts_cov
} */

/* // TODO: integrate to sapmle_weight Option<ps>
pub fn logistic_regression_covs(samples: &Samples, iteration_start: usize) -> Vec<WgtBoost> {
    // normalized unnecessary for sklearn
    let ps = vec![1.0f64; samples.samples_n()];

    logistic_regression_covs_sampleweight(samples, &ps, iteration_start)
} */

/// regression with smartcore; but not using since cannot change iteration number
// TODO: should pass WgtBoosts and add inside
pub fn logistic_regression_covs(samples: &Samples, iteration_start: usize) -> Vec<WgtBoost> {
    log::debug!("use smartcore for logreg_cov.");
    let wgts = regression::logistic_regression_covs(samples);

    let mut wgts_cov: Vec<WgtBoost> = Vec::new();
    for (wgt_i, wgt) in wgts.into_iter().enumerate() {
        let iteration = iteration_start + wgt_i;
        let wgt_boost_cov = WgtBoost::new_wgt_iteration(wgt, iteration);
        wgts_cov.push(wgt_boost_cov);
    }

    wgts_cov
}
