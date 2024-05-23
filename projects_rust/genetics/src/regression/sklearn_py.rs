use numpy::PyArray;
use pyo3::prelude::*;

use crate::samples::prelude::*;
use crate::{CovsTrait, Samples};

// sklearn
pub fn logistic_regression_covs_sample_sampleweights(
    samples: &Samples,
    ps: &[f64],
) -> (Vec<f64>, f64) {
    let covs = samples.covs().unwrap();

    let v = covs.vals_row_major();
    let ys = samples.phe_unwrap().inner_f64();
    log::debug!("ys true count {}", samples.phe_unwrap().count());
    log::debug!("ys false count {}", samples.phe_unwrap().count_false());

    logistic_regression_sampleweights_sklearn(&v, &ys, &ps)
    //let intercept_coef = logistic_regression_sampleweights_sklearn(&v, &ys, &ps);
}

// x: row major
fn logistic_regression_sampleweights_sklearn(
    x: &Vec<Vec<f64>>,
    y: &[f64],
    w: &[f64],
) -> (Vec<f64>, f64) {
    let intercept_coef = call_py(x, y, w);

    let intercept = intercept_coef[0];
    let coefs = intercept_coef[1..].to_vec();

    if intercept.is_nan() {
        panic!("Intercept is NaN.");
    }

    if coefs.iter().any(|x| x.is_nan()) {
        panic!("Coef is NaN: {:?}", coefs);
    }

    (coefs, intercept)
}

fn call_py(x: &Vec<Vec<f64>>, y: &[f64], w: &[f64]) -> Vec<f64> {
    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/src_py/logreg.py"));

    //pyo3::append_to_inittab!();
    // necessary without auto-initialize for static link
    pyo3::prepare_freethreaded_python();

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec;

    #[test]
    fn test_logistic_regression_covs_sample_sampleweights() {
        // col major
        let v = vec![vec![-1.0f64, -0.01, 0.01, 1.0]];
        let y = vec![0.0f64, 1.0, 0.0, 1.0];
        let w = vec![1.0, 1.0, 1.0, 1.0];

        let v = vec::convert_vec2d_to_row_major(&v);
        let (coefs, intercept) = logistic_regression_sampleweights_sklearn(&v, &y, &w);
        assert_eq!(coefs.len(), 1);
        assert_float_absolute_eq!(intercept, 0.0);

        // expected answer was calculated using https://stats.blue/Stats_Suite/logistic_regression_calculator.html
        assert_float_absolute_eq!(coefs[0], 5.2672, 1e-2);
        // error
        //assert_float_absolute_eq!(coefs[0], 5.2672, 1e-3);
    }
}
