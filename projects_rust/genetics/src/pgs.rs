




// mysmartcore v0.3.1
use mysmartcore::linalg::basic::matrix::DenseMatrix;
use mysmartcore::linear::logistic_regression::LogisticRegression;
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


fn score_to_prob(s: f64) -> f64 {
    1.0f64 / (1.0 + (-s).exp())
}


/* /// y=0/1
fn prob_y_ll(p: f64, y: f64) -> f64 {
    y*p.ln()+(1.0-y)*(1.0-p).ln()
} */


/// y=0/1
fn prob_y_ll(p: f64, y: i32) -> f64 {
    (y as f64)*p.ln()+(1.0-y as f64)*(1.0-p).ln()
}


// from linfa
/// A numerically stable version of the log of the logistic function.
///
/// Taken from scikit-learn
/// https://github.com/scikit-learn/scikit-learn/blob/0.23.1/sklearn/utils/_logistic_sigmoid.pyx
///
/// See the blog post describing this implementation:
/// http://fa.bianp.net/blog/2013/numerical-optimizers-for-logistic-regression/
//fn log_logistic<F: linfa::Float>(x: F) -> F {
//    if x > F::zero() {
//        -(F::one() + (-x).exp()).ln()
//    } else {
//        x - (F::one() + x.exp()).ln()
//    }
//}




/* // linfa
pub fn compute_loglikelihood(phe:&Vec<i32>, score:&[f64])->f64{
    // TMP
    0.0
} */



// mysmartcore v0.3.1 or mysmartcore
pub fn compute_loglikelihood(phe:&Vec<i32>, score:&[f64])->f64{
    //let v=vec![score.to_vec()];
   let v=score.iter().map(|x| vec![*x]).collect::<Vec<Vec<f64>>>();
    let mt = DenseMatrix::from_2d_vec( &v);
    //let mt = DenseMatrix::from_vec(score.len(),1, &score);

    // confirmed to output the same answer as sklearn
    let lr = LogisticRegression::fit(&mt, phe, Default::default()).unwrap();

    let intercept = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();
    assert_eq!(intercept.len(), 1);
    //let mut intercept=Vec::with_capacity(1);
    //lr.intercept().copy_col_as_vec(0, &mut intercept);
    let intercept=intercept[0];
    if intercept.is_nan(){
        panic!("Intercept is NaN.");
    }

    //let intercept = lr.intercept().clone().to_row_vector()[0];
    //assert_eq!(lr.coefficients().clone().to_row_vector().len(),1);

    let coef = lr.coefficients().iter().map(|x| *x).collect::<Vec<f64>>();
    //let mut coef=Vec::with_capacity(1);
    //lr.coefficients().copy_col_as_vec(0, &mut coef);
    //let coef = lr.coefficients().clone().to_row_vector()[0];
    assert_eq!(coef.len(), 1);
    let coef=coef[0];


    let score_reg_iter=score.iter().map(|s| s*coef+intercept);

    // concat and create numerically stable ver.
    let prob_iter=score_reg_iter.map(|s| score_to_prob(s));
    let llf: f64=prob_iter.zip(phe.iter()).map(|(p,y)| prob_y_ll(p,*y)).sum();

    llf
}


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


/// phe=0/1
pub fn nagelkerke_r2(phe: &Vec<i32>, score: &[f64],score_cov: &[f64]) -> f64 {

    assert_eq!(phe.len(),score.len());
    assert_eq!(phe.len(),score_cov.len());

    
    let llf=compute_loglikelihood(phe, score);
    let llnull=compute_loglikelihood(phe, score_cov);

    let n=phe.len() as f64;
    let nagel=(1.0f64-((llnull-llf)*(2.0/n)).exp())/(1.0-(llnull*(2.0/n)).exp());

    nagel

}



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
    fn test_nagelkerke_r2(){
        // eff is nan
        //let phe=vec![0i32,0,1,1];
        let phe=vec![0i32,1,0,1];
        //let phe=phe.iter().map(|x| *x as f64).collect::<Vec<f64>>();
        let score=vec![0.1,0.2,0.3,0.9];
        let score_cov=vec![0.1,0.2,0.3,0.9];
        let nagel=nagelkerke_r2(&phe, &score, &score_cov);

        assert_float_absolute_eq!(nagel,0.0);

    }
}
