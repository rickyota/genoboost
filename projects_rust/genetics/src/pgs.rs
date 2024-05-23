// mysmartcore v0.3.1
//use mysmartcore::linalg::basic::matrix::DenseMatrix;
//use mysmartcore::linear::logistic_regression::LogisticRegression;
//use mysmartcore::linalg::basic::arrays::Array2;

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

/*
use numpy::PyArray;
use pyo3::prelude::*;
*/

use crate::regression;
use crate::vec;

// unstable
//fn score_to_prob(s: f64) -> f64 {
//    1.0f64 / (1.0 + (-s).exp())
//}

/* /// y=0/1
fn prob_y_ll(p: f64, y: f64) -> f64 {
    y*p.ln()+(1.0-y)*(1.0-p).ln()
} */

/// y=0/1
//fn prob_y_ll(p: f64, y: i32) -> f64 {
//    (y as f64) * p.ln() + (1.0 - y as f64) * (1.0 - p).ln()
//}

/// unstable
/// y=0/1
//fn prob_y_ll_unstable(p: f64, y: bool) -> f64 {
//    (y as i32 as f64) * p.ln() + (1.0 - y as i32 as f64) * (1.0 - p).ln()
//}

/// s -> log(p)
/// where p=1/(1+exp(-s))
/// numerically stable
///
// from linfa
/// A numerically stable version of the log of the logistic function.
///
/// Taken from scikit-learn
/// https://github.com/scikit-learn/scikit-learn/blob/0.23.1/sklearn/utils/_logistic_sigmoid.pyx
///
/// See the blog post describing this implementation:
/// http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression/
///fn log_logistic<F: linfa::Float>(x: F) -> F {
///    if x > F::zero() {
///        -(F::one() + (-x).exp()).ln()
///    } else {
///        x - (F::one() + x.exp()).ln()
///    }
///}
fn log_prob(s: f64) -> f64 {
    if s > 0.0 {
        -(1.0 + (-s).exp()).ln()
    } else {
        s - (1.0 + s.exp()).ln()
    }
}

/// s -> log(1-p)
/// numerically stable
fn log_1_prob(s: f64) -> f64 {
    -s + log_prob(s)
}

/// (s,y) -> y*ln(p)+(1-y)*ln(1-p)
fn prob_s_y_ll(s: f64, y: bool) -> f64 {
    if y {
        log_prob(s)
    } else {
        log_1_prob(s)
    }
}

//pub fn vals_norm(vals: &[f64]) -> Vec<f64> {
//    let vals_mean = vals.iter().sum::<f64>() / vals.len() as f64;
//    let vals_std = vals
//        .iter()
//        .map(|x| (x - vals_mean).powi(2))
//        .sum::<f64>()
//        .sqrt();
//    let vals_norm = vals
//        .iter()
//        .map(|x| (x - vals_mean) / vals_std)
//        .collect::<Vec<f64>>();
//
//    vals_norm
//}

// moved to vec.rs
//pub fn vals_norm(vals: &[f64], mean: f64, std: f64) -> Vec<f64> {
//    let vals_norm = vals.iter().map(|x| (x - mean) / std).collect::<Vec<f64>>();
//
//    vals_norm
//}

// mylinfa v0.7.0
// score len=n without padding
pub fn compute_loglikelihood(phe: &Vec<bool>, score: &[f64]) -> f64 {
    //let score_reg_iter: Vec<f64> = if vec::is_all_same_float(score) {
    let score_reg_iter: Vec<f64> = if vec::std_vector(score) < 1e-10 {
        log::debug!(
            "Skip regression for score since variance is small: {:?}",
            vec::std_vector(score)
        );
        //vec::is_all_same_float(score) {
        // regression without intercept
        // let v = vec![score.to_vec()];
        let coef =
            regression::logistic_regression_vec_no_intercept_1col(score.to_vec(), phe.clone());

        let score_reg_iter = score.iter().map(|s| s * coef).collect::<Vec<f64>>();
        score_reg_iter
    } else {
        let v = vec![score.to_vec()];
        let (coefs, intercept) = regression::logistic_regression_vec2(v, phe.clone());
        let coef = coefs[0];

        let score_reg_iter = score
            .iter()
            .map(|s| s * coef + intercept)
            .collect::<Vec<f64>>();
        score_reg_iter
    };

    //let v = vec![score.to_vec()];
    //let (coefs, intercept) = regression::logistic_regression_vec2(v, phe.clone());
    //let coef = coefs[0];

    //let score_reg_iter = score.iter().map(|s| s * coef + intercept);

    let llf: f64 = score_reg_iter
        .iter()
        .zip(phe.iter())
        .map(|(s, y)| prob_s_y_ll(*s, *y))
        .sum();

    llf
}

//// mysmartcore v0.3.1 or mysmartcore
//// score len=n without padding
//pub fn compute_loglikelihood(phe: &Vec<i32>, score: &[f64]) -> f64 {
//    //let v=vec![score.to_vec()];
//    let score_mean = score.iter().sum::<f64>() / score.len() as f64;
//    let score_std = score
//        .iter()
//        .map(|x| (x - score_mean).powi(2))
//        .sum::<f64>()
//        .sqrt();
//    let score_norm = vec::norm_vec(score, score_mean, score_std);
//    //let score_norm = vals_norm(score, score_mean, score_std);
//    let v = score_norm
//        .iter()
//        .map(|x| vec![*x])
//        .collect::<Vec<Vec<f64>>>();
//    //let v = score.iter().map(|x| vec![*x]).collect::<Vec<Vec<f64>>>();
//    let mt = DenseMatrix::from_2d_vec(&v);
//    //let mt = DenseMatrix::from_vec(score.len(),1, &score);
//
//    // confirmed to output the same answer as sklearn
//    let lr = LogisticRegression::fit(&mt, phe, Default::default()).unwrap();
//
//    let intercept_norm = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();
//    assert_eq!(intercept_norm.len(), 1);
//    //let mut intercept=Vec::with_capacity(1);
//    //lr.intercept().copy_col_as_vec(0, &mut intercept);
//    let intercept_norm = intercept_norm[0];
//    if intercept_norm.is_nan() {
//        panic!("Intercept for loglikelihood is NaN.");
//    }
//
//    //let intercept = lr.intercept().clone().to_row_vector()[0];
//    //assert_eq!(lr.coefficients().clone().to_row_vector().len(),1);
//
//    let coef_norm = lr.coefficients().iter().map(|x| *x).collect::<Vec<f64>>();
//    //let mut coef=Vec::with_capacity(1);
//    //lr.coefficients().copy_col_as_vec(0, &mut coef);
//    //let coef = lr.coefficients().clone().to_row_vector()[0];
//    assert_eq!(coef_norm.len(), 1);
//    let coef_norm = coef_norm[0];
//
//    // c=c' / std
//    let coef = coef_norm / score_std;
//    // i = i' - ( c' * m / std) = i' - ( c * m )
//    let intercept = intercept_norm - coef * score_mean;
//
//    let score_reg_iter = score.iter().map(|s| s * coef + intercept);
//
//    // concat and create numerically stable ver.
//    let prob_iter = score_reg_iter.map(|s| score_to_prob(s));
//    let llf: f64 = prob_iter
//        .zip(phe.iter())
//        .map(|(p, y)| prob_y_ll(p, *y))
//        .sum();
//
//    llf
//}

/* // mysmartcore v0.2.0
pub fn compute_loglikelihood(phe:&Vec<f64>, score:&[f64])->f64{
    let mt = DenseMatrix::from_vec(score.len(),1, &score);

    // default: alpha=0.0
    // confirmed to output the same answer as sklearn
    let lr = LogisticRegression::fit(&mt, phe, Default::default()).unwrap();

    let intercept = lr.intercept().clone().to_row_vector()[0];
    assert_eq!(lr.coefficients().clone().to_row_vector().len(),1);
    let coef = lr.coefficients().clone().to_row_vector()[0];

    let score_reg_iter=score.iter().map(|s| s*coef+intercept);

    let prob_iter=score_reg_iter.map(|s| score_to_prob(s));

    let llf: f64=prob_iter.zip(phe.iter()).map(|(p,y)| prob_y_ll(p,*y)).sum();

    llf
} */

pub fn nagelkerke_r2(phe: &Vec<bool>, score: &[f64], score_cov: &[f64]) -> f64 {
    assert_eq!(phe.len(), score.len());
    assert_eq!(phe.len(), score_cov.len());

    log::debug!("score_cov[0],[1] {:?}, {:?}", score_cov[0], score_cov[1]);

    let llf = compute_loglikelihood(phe, score);
    let llnull = compute_loglikelihood(phe, score_cov);

    let n = phe.len() as f64;
    let nagel = (1.0f64 - ((llnull - llf) * (2.0 / n)).exp()) / (1.0 - (llnull * (2.0 / n)).exp());

    nagel
}

/*
// mysmartcore v0.3.0
/// phe=0/1
pub fn nagelkerke_r2(phe: &Vec<i32>, score: &[f64], score_cov: &[f64]) -> f64 {
    assert_eq!(phe.len(), score.len());
    assert_eq!(phe.len(), score_cov.len());

    let llf = compute_loglikelihood(phe, score);
    let llnull = compute_loglikelihood(phe, score_cov);

    let n = phe.len() as f64;
    let nagel = (1.0f64 - ((llnull - llf) * (2.0 / n)).exp()) / (1.0 - (llnull * (2.0 / n)).exp());

    nagel
}
*/

/* /// smartcore v0.2.0
/// phe=0/1
pub fn nagelkerke_r2(phe: &Vec<f64>, score: &[f64],score_cov: &[f64]) -> f64 {

    assert_eq!(phe.len(),score.len());
    assert_eq!(phe.len(),score_cov.len());


    let llf=compute_loglikelihood(phe, score);
    let llnull=compute_loglikelihood(phe, score_cov);

    let n=phe.len() as f64;
    let nagel=(1.0f64-((llnull-llf)*(2.0/n)).exp())/(1.0-(llnull*(2.0/n)).exp());

    nagel

}
 */

/*
pub fn nagelkerke_r2(phe: &[f64], score: &[f64],score_cov: &[f64]) -> f64 {
    let py_app = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/src_py/accuracy.py"));

    pyo3::prepare_freethreaded_python();
    Python::with_gil(|py| -> f64 {
        PyModule::from_code(py, py_app, "accuracy", "accuracy").unwrap();
        let app: Py<PyAny> = PyModule::from_code(py, py_app, "", "")
            .unwrap()
            .getattr("nagelkerke_r2")
            .unwrap()
            .into();
        //let res = app.call0(py);
        // [ref](https://docs.rs/numpy/latest/numpy/array/struct.PyArray.html#method.from_vec2)
        let args = (
            PyArray::from_slice(py, phe),
            PyArray::from_slice(py, score),
            PyArray::from_slice(py, score_cov),
        );
        //let args=(PyArray::from_vec2(py, x).unwrap(),PyArray::from_vec(py,y),Pyarray::from_vec(py,w));
        //let args=(PyArray::from_vec2(py, &x).unwrap(), y.into_pyarray(py),w.into_pyarray(py));
        let res = app.call1(py, args);
        //let res = app.call1(py,(v ,));

        let ans: f64 = res.unwrap().as_ref(py).extract().unwrap();
        //let ans=res.unwrap().as_ref(py);
        //let ans: Array1<f64> =res.unwrap().as_ref(py).extract().unwrap();
        //let ans: PyArray1<f64> =res.unwrap().as_ref(py);

        ans
    })
}
 */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nagelkerke_r2() {
        // eff is nan
        //let phe=vec![0i32,0,1,1];
        let phe = vec![false, true, false, true];
        //let phe = vec![0i32, 1, 0, 1];
        //let phe=phe.iter().map(|x| *x as f64).collect::<Vec<f64>>();
        let score = vec![0.1, 0.2, 0.3, 0.9];
        let score_cov = vec![0.1, 0.2, 0.3, 0.9];
        let nagel = nagelkerke_r2(&phe, &score, &score_cov);

        assert_float_absolute_eq!(nagel, 0.0);
    }

    fn score_to_prob(s: f64) -> f64 {
        1.0f64 / (1.0 + (-s).exp())
    }

    fn prob_y_ll_unstable(p: f64, y: bool) -> f64 {
        (y as i32 as f64) * p.ln() + (1.0 - y as i32 as f64) * (1.0 - p).ln()
    }

    fn log_prob_unstable(s: f64) -> f64 {
        -(1.0 + (-s).exp()).ln()
    }

    #[test]
    fn test_prob_y_ll() {
        let s = 0.1f64;
        let y = true;

        let ll = prob_s_y_ll(s, y);

        let p = score_to_prob(s);
        let ll_unstable = prob_y_ll_unstable(p, y);
        assert_float_absolute_eq!(ll, ll_unstable);
    }

    #[test]
    fn test_log_prob() {
        let s = -1000000.0f64;
        // lp=-1000000.0f64
        let lp = log_prob(s);
        assert!(!lp.is_infinite());

        // lp_unstable=-inf
        let lp_unstable = log_prob_unstable(s);
        assert!(lp_unstable.is_infinite());
        assert!(lp_unstable.is_sign_negative());
    }
}
