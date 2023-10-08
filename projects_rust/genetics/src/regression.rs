//use linfa::prelude::*;
//use linfa::traits::Fit;
//use linfa_logistic::LogisticRegression;
//use ndarray::prelude::*;

// mysmartcore v0.3.1
//use mysmartcore::linalg::basic::arrays::Array2;
use mysmartcore::linalg::basic::matrix::DenseMatrix;
use mysmartcore::linear::logistic_regression::LogisticRegression;

// mysmartcore v0.2.1
//use mysmartcore::linalg::naive::dense_matrix::DenseMatrix;
//use mysmartcore::linalg::BaseMatrix;
//use mysmartcore::linear::logistic_regression::LogisticRegression;

// v.0.3.2
//use smartcore::linalg::basic::matrix::DenseMatrix;
//use smartcore::linear::logistic_regression::LogisticRegression;
//use smartcore::linalg::basic::arrays::Array2;

// v0.2.0
//use smartcore::linalg::naive::dense_matrix::DenseMatrix;
//use smartcore::linalg::BaseMatrix;
//use smartcore::linear::logistic_regression::LogisticRegression;

use crate::samples::prelude::*;
use crate::CovId;
use crate::{samples::CovsTrait, Samples};
use crate::{wgt::Coef, Wgt};

//use numpy::PyArray;
//use pyo3::prelude::*;

//mod logistic;

//smartcore v3 or mysmartcore
pub fn logistic_regression_covs(samples: &Samples) -> Vec<Wgt> {
    //let covs_val = samples.covs();
    //// TODO: return None
    //if let None = covs_val {
    //    return vec![];
    //}
    if let Some(covs) = samples.covs() {
        // load normed vals
        let v = covs.vals_row_major_norm();
        let ys = samples.phe_unwrap().inner_i32();

        let cov_n = covs.covs_n();
        //let cov_n=v.len();
        log::debug!("cov_n {}", cov_n);
        log::debug!("ys true count {}", samples.phe_unwrap().count());
        log::debug!("ys false count {}", samples.phe_unwrap().count_false());

        log::debug!("norm v");
        log::debug!("norm mt {}, {}", v[0][0], v[0][1]);

        // dense_matrix from 2d array
        // https://docs.rs/smartcore/latest/src/smartcore/linalg/naive/dense_matrix.rs.html#248
        //let (v, row_n, col_n) = covs_val.vals_row_major_vec();
        //log::debug!("row, col {}, {}",row_n,col_n);
        let mt = DenseMatrix::from_2d_vec(&v);
        //let mt = DenseMatrix::from_2d_vec(&covs_val.vals_row_major());
        // default: alpha=0.0
        // confirmed to output the same result as sklearn
        let lr = LogisticRegression::fit(&mt, &ys, Default::default()).unwrap();
        log::debug!("logreg");
        log::debug!("intercept {:?}", lr.intercept());
        log::debug!("coef {:?}", lr.coefficients());

        //revert normalization
        let means = samples.covs().unwrap().means();
        let stds = samples.covs().unwrap().stds();

        let coefs_norm = lr.coefficients().iter().map(|x| *x).collect::<Vec<f64>>();
        assert_eq!(coefs_norm.len(), covs.covs_n());

        // c=c' / std
        let coefs = coefs_norm
            .iter()
            .zip(stds.iter())
            .map(|(c, s)| c / s)
            .collect::<Vec<f64>>();
        //let coef=coefs_norm.iter().zip(means.iter()).zip(stds.iter()).map(|((c,m),s)| )

        let intercept_norm = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();
        log::debug!("intercept_norm {:?}", intercept_norm);
        assert_eq!(intercept_norm.len(), 1);
        let intercept_norm = intercept_norm[0];
        if intercept_norm.is_nan() {
            panic!("Intercept is NaN.");
        }

        // i = i' - sum( c' * m / std) = i' - sum( c * m )
        let intercept = intercept_norm
            + coefs
                .iter()
                .zip(means.iter())
                .map(|(c, m)| -c * m)
                .sum::<f64>();

        let mut wgts_cov: Vec<Wgt> = Vec::new();

        //let mut intercept=vec![0.0f64;1];
        //let mut intercept=Vec::with_capacity(1);
        //lr.intercept().copy_col_as_vec(0, &mut intercept);
        //let intercept = lr.intercept().clone().to_row_vector()[0];

        //let intercept=lr.intercept().col_iter().collect::<Vec<f64>>();
        //let intercept = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();
        //log::debug!("intercept {:?}", intercept);
        //assert_eq!(intercept.len(), 1);
        //let intercept = intercept[0];
        if intercept.is_nan() {
            panic!("Intercept is NaN.");
        }
        let cov_id = CovId::new_const();
        let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));
        wgts_cov.push(wgt_const);

        let cov_indexs = covs.cov_indexs().unwrap();
        for wgt_i in 0..covs.covs_n() {
            let coef = coefs[wgt_i];
            if coef.is_nan() {
                panic!("coef is NaN.");
            }
            let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
            let wgt_cov = Wgt::construct_cov(cov_id, Coef::Linear(coef));
            wgts_cov.push(wgt_cov);
        }

        wgts_cov
    } else {
        //    return vec![];
        // regress on const

        log::debug!("Regression on const value only.");

        log::debug!("ys true count {}", samples.phe_unwrap().count());
        log::debug!("ys false count {}", samples.phe_unwrap().count_false());

        let y1_n = samples.phe_unwrap().count() as f64;
        let y0_n = samples.phe_unwrap().count_false() as f64;
        // TODO: test
        let intercept = (y1_n / y0_n).ln();

        // cannot run on smartcore since automatically adds const values
        //let ys = samples.phe().inner_i32();
        //let sample_n = samples.samples_n();
        //let v = vec![vec![1.0f64; sample_n]];
        //let mt = DenseMatrix::from_2d_vec(&v);
        //let lr = LogisticRegression::fit(&mt, &ys, Default::default()).unwrap();

        let mut wgts_cov: Vec<Wgt> = Vec::new();
        if intercept.is_nan() {
            panic!("Intercept is NaN.");
        }
        let cov_id = CovId::new_const();
        let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));
        wgts_cov.push(wgt_const);

        wgts_cov
    }
}

/* // smartcore v0.2.0
pub fn logistic_regression_covs(samples: &Samples) -> Vec<Wgt> {
    let covs_val = samples.covs();
    // TODO: return None
    if let None = covs_val {
        return vec![];
    }

    let covs_val = covs_val.unwrap();

    let ys = samples.phe().inner_f64();
    // dense_matrix from 2d array
    // https://docs.rs/smartcore/latest/src/smartcore/linalg/naive/dense_matrix.rs.html#248
    let (v, row_n, col_n) = covs_val.vals_row_major_vec();
    log::debug!("row, col {}, {}",row_n,col_n);
    let mt = DenseMatrix::from_vec(row_n, col_n, &v);
    //let mt = DenseMatrix::from_2d_vec(&covs_val.vals_row_major());
    // default: alpha=0.0
    // confirmed to output the same answer as sklearn
    let lr = LogisticRegression::fit(&mt, &ys, Default::default()).unwrap();
    log::debug!("logreg");
    log::debug!("intercept {}", lr.intercept());
    log::debug!("coef {}", lr.coefficients());

    let mut wgts_cov: Vec<Wgt> = Vec::new();

    let intercept = lr.intercept().clone().to_row_vector()[0];
    if intercept.is_nan(){
        panic!("Intercept is NaN.");
    }
    let cov_id = CovId::new_const();
    let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));
    wgts_cov.push(wgt_const);

    let coefs = lr.coefficients().clone().to_row_vector();
    assert_eq!(coefs.len(), covs_val.covs_n());
    let cov_indexs = covs_val.cov_indexs().unwrap();
    for wgt_i in 0..covs_val.covs_n() {
        let coef = coefs[wgt_i];
        if coef.is_nan(){
            panic!("coef is NaN.");
        }
        let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        let wgt_cov = Wgt::construct_cov(cov_id, Coef::Linear(coef));
        wgts_cov.push(wgt_cov);
    }

    wgts_cov
}
 */

/*
// sklearn
pub fn logistic_regression_covs_sampleweights_sklearn(samples: &Samples, ps: &[f64]) -> Vec<Wgt> {
    let covs_val = samples.covs();
    if let None = covs_val {
        return vec![];
    }

    let covs_val = covs_val.unwrap();
    let v = covs_val.vals_row_major();
    let ys = samples.phe().inner_f64();

    let intercept_coef = logistic_regression_sampleweights_sklearn(&v, &ys, &ps);

    let intercept = intercept_coef[0];
    let coefs = &intercept_coef[1..];

    let mut wgts_cov: Vec<Wgt> = Vec::new();

    if intercept.is_nan() {
        panic!("Intercept is NaN.");
    }
    let cov_id = CovId::new_const();
    let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));
    wgts_cov.push(wgt_const);

    assert_eq!(coefs.len(), covs_val.covs_n());
    let cov_indexs = covs_val.cov_indexs().unwrap();
    for wgt_i in 0..covs_val.covs_n() {
        let coef = coefs[wgt_i];
        if coef.is_nan() {
            panic!("coef is NaN.");
        }
        let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        let wgt_cov = Wgt::construct_cov(cov_id, Coef::Linear(coef));
        wgts_cov.push(wgt_cov);
    }

    wgts_cov
}

fn logistic_regression_sampleweights_sklearn(x: &Vec<Vec<f64>>, y: &[f64], w: &[f64]) -> Vec<f64> {
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

/* // linfa
// [ref](https://github.com/zupzup/rust-ml-example/blob/main/src/main.rs)
pub fn logistic_regression_covs(samples: &Samples) -> Vec<Wgt> {
    let covs_val = samples.covs();
    // TODO: return None
    if let None = covs_val {
        return vec![];
    }

    let covs_val = covs_val.unwrap();
    //let v = covs_val.vals_row_major();
    let (v, row_n, col_n) = covs_val.vals_norm_row_major_vec();
    //let (v, row_n, col_n) = covs_val.vals_row_major_vec();
    //let (v, row_n, col_n) = covs_val.vals_column_major_vec();
    let ys = samples.phe().inner_bool();

    let cov_n = covs_val.covs_n();
    //let cov_n=v.len();
    log::debug!("cov_n {}", cov_n);
    assert_eq!(col_n, cov_n);

    let mt: Array2<f64> = Array::from_shape_vec((row_n, col_n), v).unwrap();
    log::debug!("mt {}, {}, {}", mt[[0, 0]], mt[[0, 1]], mt[[0, 2]]);
    log::debug!("mt {}, {}, {}", mt[[row_n-1, 0]], mt[[row_n-1, 1]], mt[[row_n-1, 2]]);
    let ys = Array1::from(ys);
    let d = Dataset::new(mt, ys);
    //let d = Dataset::new(mt, ys).map_targets(|x| if *x { "d" } else { "n" });

    // try 1.0
    println!("init 0.0");
    let init = Array1::from(vec![0.0f64; cov_n + 1]);
    //let init=Array1::from(vec![0.0f64;cov_n+1]);

    let max_iter = 10000;
    //let max_iter = 100;
    let tol = 1e-4;
    //let tol = 1e-8;
    let lr = LogisticRegression::default()
        .max_iterations(max_iter)
        .gradient_tolerance(tol)
        .alpha(0.0)
        .initial_params(init)
        .fit(&d)
        .expect("Could not train model.");

    println!("iter: {}", max_iter);
    println!("tol: {}", tol);
    println!("alpha: {}", 0.0);

    /*     // dense_matrix from 2d array
    // https://docs.rs/smartcore/latest/src/smartcore/linalg/naive/dense_matrix.rs.html#248
    //let (v, row_n, col_n) = covs_val.vals_row_major_vec();
    //log::debug!("row, col {}, {}",row_n,col_n);
    let mt = DenseMatrix::from_2d_vec(&v);
    //let mt = DenseMatrix::from_2d_vec(&covs_val.vals_row_major());
    // default: alpha=0.0
    // confirmed to output the same answer as sklearn
    let lr = LogisticRegression::fit(&mt, &ys, Default::default()).unwrap(); */
    log::debug!("logreg");
    log::debug!("intercept {:?}", lr.intercept());
    log::debug!("coef {:?}", lr.params());

    let mut wgts_cov: Vec<Wgt> = Vec::new();

    //let mut intercept=vec![0.0f64;1];
    //let mut intercept=Vec::with_capacity(1);
    //lr.intercept().copy_col_as_vec(0, &mut intercept);
    //let intercept = lr.intercept().clone().to_row_vector()[0];

    //let intercept=lr.intercept().col_iter().collect::<Vec<f64>>();
    //let intercept = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();

    let means = samples.covs().unwrap().means();
    let stds = samples.covs().unwrap().stds();

    let coefs_norm = lr.params().iter().map(|x| *x).collect::<Vec<f64>>();
    assert_eq!(coefs_norm.len(), covs_val.covs_n());

    // c=c' / std
    let coefs = coefs_norm
        .iter()
        .zip(stds.iter())
        .map(|(c, s)| c / s)
        .collect::<Vec<f64>>();

    //let intercept_norm = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();
    let intercept_norm = lr.intercept();
    log::debug!("intercept_norm {:?}", intercept_norm);
    if intercept_norm.is_nan() {
        panic!("Intercept is NaN.");
    }

    // i = i' - sum( c' * m / std) = i' - sum( c * m )
    let intercept = intercept_norm
        + coefs
            .iter()
            .zip(means.iter())
            .map(|(c, m)| -c * m)
            .sum::<f64>();


    //let intercept = lr.intercept();
    //log::debug!("intercept {:?}", intercept);
    //if intercept.is_nan() {
    //    panic!("Intercept is NaN.");
    //}
    let cov_id = CovId::new_const();
    let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));
    wgts_cov.push(wgt_const);

    //let coefs = lr.params().iter().map(|x| *x).collect::<Vec<f64>>();
    ////let coefs = lr.coefficients().clone().to_row_vector();
    //assert_eq!(coefs.len(), covs_val.covs_n());
    let cov_indexs = covs_val.cov_indexs().unwrap();
    for wgt_i in 0..covs_val.covs_n() {
        let coef = coefs[wgt_i];
        if coef.is_nan() {
            panic!("coef is NaN.");
        }
        let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
        let wgt_cov = Wgt::construct_cov(cov_id, Coef::Linear(coef));
        wgts_cov.push(wgt_cov);
    }

    wgts_cov
} */

#[cfg(test)]
mod tests {
    //use super::*;
    /*
    #[test]
    fn simple_example_1() {
        let log_reg = LogisticRegression::default();
        let x = array![[-1.0], [-0.01], [0.01], [1.0]];
        let y = array![0, 0, 1, 1];
        let dataset = Dataset::new(x, y);
        let res = log_reg.fit(&dataset).unwrap();
        assert_abs_diff_eq!(res.intercept(), 0.0);
        assert!(res.params().abs_diff_eq(&array![0.681], 1e-3));
        assert_eq!(
            &res.predict(dataset.records()),
            dataset.targets().as_single_targets()
        );
    } */
}
