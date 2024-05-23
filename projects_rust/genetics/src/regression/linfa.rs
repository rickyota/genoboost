// It seems that mylinfa can use f128.

// mylinfa v0.7.0
use mylinfa::prelude::*;
use mylinfa::traits::Fit;
use mylinfa_logistic::LogisticRegression;
use ndarray::prelude::*;

// linfa v0.7.0
//use linfa::prelude::*;
//use linfa::traits::Fit;
//use linfa_logistic::LogisticRegression;
//use ndarray::prelude::*;

use crate::samples::prelude::*;
use crate::{CovsTrait, Samples};

// use means stored in samples
pub fn logistic_regression_covs_samples(samples: &Samples) -> (Vec<f64>, f64) {
    let covs = samples.covs().unwrap();

    let (v, (row_n, col_n)) = covs.vals_norm_vec1d_row_major();
    let ys_bool = samples.phe_unwrap().inner_bool();
    log::debug!("ys true count {}", samples.phe_unwrap().count());
    log::debug!("ys false count {}", samples.phe_unwrap().count_false());

    let cov_n = covs.covs_n();
    log::debug!("cov_n {}", cov_n);
    assert_eq!(col_n, cov_n);

    let means = samples.covs().unwrap().means();
    let stds = samples.covs().unwrap().stds();
    log::debug!("means {:?}", means);
    log::debug!("stds {:?}", stds);

    logistic_regression_vec_linfa_normalize(v, col_n, row_n, ys_bool, &means, &stds)
}

// mylinfa v0.7.0
pub fn logistic_regression_vec_linfa_normalize(
    v: Vec<f64>, //row major; 1d
    col_n: usize,
    row_n: usize,
    ys: Vec<bool>,
    means: &Vec<f64>,
    stds: &Vec<f64>,
) -> (Vec<f64>, f64) {
    log::debug!("logreg");

    v.iter().for_each(|x| {
        if !x.is_finite() {
            panic!("Not finite value is included in vector: {:?}", x);
        }
    });

    stds.iter().for_each(|x| {
        if !x.is_finite() || *x == 0.0 {
            panic!(
                "Not finite value or 0.0 is included in stds: {:?}. Remove the column.",
                x
            );
        }
    });

    let mt: Array2<f64> = Array::from_shape_vec((row_n, col_n), v).unwrap();
    //let mt: Array2<f64> = Array::from_shape_vec((col_n, row_n), v).unwrap();
    //let mt: Array2<f64> = Array::from_shape_vec((row_n, col_n), v).unwrap();
    assert_eq!(means.len(), col_n);
    assert_eq!(stds.len(), col_n);

    //println!("mt {:?}", mt.shape());
    //println!("col {:?}", col_n);
    //println!("row {:?}", row_n);

    //log::debug!("mt {}, {}, {}", mt[[0, 0]], mt[[0, 1]], mt[[0, 2]]);
    //log::debug!(
    //    "mt {}, {}, {}",
    //    mt[[row_n - 1, 0]],
    //    mt[[row_n - 1, 1]],
    //    mt[[row_n - 1, 2]]
    //);

    let y = Array1::from(ys.clone());
    //let ys = Array1::from(ys_bool);
    let d = Dataset::new(mt, y);
    //let d = Dataset::new(mt, ys).map_targets(|x| if *x { "d" } else { "n" });

    let cov_n = means.len();

    assert_eq!(cov_n, col_n);

    // try 1.0
    //println!("init 0.0");
    let init = Array1::from(vec![0.0f64; cov_n + 1]);
    //let init=Array1::from(vec![0.0f64;cov_n+1]);

    // If the program stops here, stds(variance of mt) could be too small. Check the input.
    //let max_iter = 10000;
    // default
    let max_iter = 100;
    // default
    let tol = 1e-4;
    //let tol = 1e-8;
    let lr = LogisticRegression::default()
        .max_iterations(max_iter)
        .gradient_tolerance(tol)
        .alpha(0.0)
        .initial_params(init)
        .fit(&d)
        .expect("Could not train model.");

    //println!("iter: {}", max_iter);
    //println!("tol: {}", tol);
    //println!("alpha: {}", 0.0);

    /*     // dense_matrix from 2d array
    // https://docs.rs/smartcore/latest/src/smartcore/linalg/naive/dense_matrix.rs.html#248
    //let (v, row_n, col_n) = covs_val.vals_row_major_vec();
    //log::debug!("row, col {}, {}",row_n,col_n);
    let mt = DenseMatrix::from_2d_vec(&v);
    //let mt = DenseMatrix::from_2d_vec(&covs_val.vals_row_major());
    // default: alpha=0.0
    // confirmed to output the same answer as sklearn
    let lr = LogisticRegression::fit(&mt, &ys, Default::default()).unwrap(); */
    log::debug!("intercept {:?}", lr.intercept());
    log::debug!("coef {:?}", lr.params());
    log::debug!("labels {:?}", lr.labels());

    //let mut intercept=vec![0.0f64;1];
    //let mut intercept=Vec::with_capacity(1);
    //lr.intercept().copy_col_as_vec(0, &mut intercept);
    //let intercept = lr.intercept().clone().to_row_vector()[0];

    //let intercept=lr.intercept().col_iter().collect::<Vec<f64>>();
    //let intercept = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();

    let coefs_norm = lr.params().iter().map(|x| *x).collect::<Vec<f64>>();
    log::debug!("coef_norm {:?}", coefs_norm);
    if coefs_norm.iter().any(|x| x.is_nan()) {
        panic!("Coef is NaN: {:?}", coefs_norm);
    }

    let intercept_norm = lr.intercept();
    log::debug!("intercept_norm {:?}", intercept_norm);
    if intercept_norm.is_nan() {
        panic!("Intercept is NaN.");
    }

    //let x = linfa_logistic::label_classes(y).unwrap();
    //log::debug!("x {:?}", x);

    // label flip by linfa
    let is_label_flipped = is_label_flipped_linfa_labels(lr.labels());
    log::debug!("is_label_flipped {:?}", is_label_flipped);

    let coefs_norm = if is_label_flipped {
        coefs_norm.iter().map(|x| -*x).collect::<Vec<f64>>()
    } else {
        coefs_norm
    };

    let intercept_norm = if is_label_flipped {
        -intercept_norm
    } else {
        intercept_norm
    };

    // revert normalization
    // c=c' / std
    let coefs = coefs_norm
        .iter()
        .zip(stds.iter())
        .map(|(c, s)| c / s)
        .collect::<Vec<f64>>();

    // i = i' - sum( c' * m / std) = i' - sum( c * m )
    let intercept = intercept_norm
        + coefs
            .iter()
            .zip(means.iter())
            .map(|(c, m)| -c * m)
            .sum::<f64>();

    (coefs, intercept)
}

fn is_label_flipped_linfa_labels<F>(labels: &mylinfa_logistic::BinaryClassLabels<F, bool>) -> bool {
    //log::debug!("labels {:?}", labels.pos.class);
    //log::debug!("labels {:?}", labels.neg.class);
    // always +1.0
    //log::debug!("labels {:?}", lr.labels().pos.label);
    // always -1.0
    //log::debug!("labels {:?}", lr.labels().neg.label);

    // true if pos.class==false
    !labels.pos.class
}

// mylinfa v0.7.0
// mainly for v=score and all samples score are the same
// return coef
pub fn logistic_regression_vec_linfa_no_intercept_1col(
    v: Vec<f64>, //row major; 1d
    //col_n: usize,
    //row_n: usize,
    ys: Vec<bool>,
) -> f64 {
    log::debug!("logreg with no intercept");

    let col_n = 1;
    let row_n = v.len();

    let mt: Array2<f64> = Array::from_shape_vec((row_n, col_n), v).unwrap();
    let y = Array1::from(ys);
    //let y = Array1::from(ys.clone());
    let d = Dataset::new(mt, y);

    // try 1.0
    //println!("init 0.0");
    let init = Array1::from(vec![0.0f64; col_n]);

    // If the program stops here, stds(variance of mt) could be too small. Check the input.
    //let max_iter = 10000;
    // default
    let max_iter = 100;
    // default
    let tol = 1e-4;
    //let tol = 1e-8;
    let lr = LogisticRegression::default()
        .max_iterations(max_iter)
        .gradient_tolerance(tol)
        .alpha(0.0)
        .initial_params(init)
        .with_intercept(false)
        .fit(&d)
        .expect("Could not train model.");

    log::debug!("intercept {:?}", lr.intercept());
    log::debug!("coef {:?}", lr.params());
    log::debug!("labels {:?}", lr.labels());

    //////// not yet //////////

    //let mut intercept=vec![0.0f64;1];
    //let mut intercept=Vec::with_capacity(1);
    //lr.intercept().copy_col_as_vec(0, &mut intercept);
    //let intercept = lr.intercept().clone().to_row_vector()[0];

    //let intercept=lr.intercept().col_iter().collect::<Vec<f64>>();
    //let intercept = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();

    let coefs = lr.params().iter().map(|x| *x).collect::<Vec<f64>>();
    log::debug!("coef_norm {:?}", coefs);
    if coefs.iter().any(|x| x.is_nan()) {
        panic!("Coef is NaN: {:?}", coefs);
    }

    //let intercept_norm = lr.intercept();
    //log::debug!("intercept_norm {:?}", intercept_norm);
    //if intercept_norm.is_nan() {
    //    panic!("Intercept is NaN.");
    //}

    //let x = linfa_logistic::label_classes(y).unwrap();
    //log::debug!("x {:?}", x);

    // label flip by linfa
    let is_label_flipped = is_label_flipped_linfa_labels(lr.labels());
    log::debug!("is_label_flipped {:?}", is_label_flipped);

    let coefs = if is_label_flipped {
        coefs.iter().map(|x| -*x).collect::<Vec<f64>>()
    } else {
        coefs
    };

    //let intercept_norm = if is_label_flipped {
    //    -intercept_norm
    //} else {
    //    intercept_norm
    //};

    //(coefs, intercept)
    //(coefs, intercept)
    let coef = coefs[0];
    coef
}

// TODO: move to regression.rs and merge with logistic_regression_const_samples()
#[allow(dead_code)]
fn logistic_regression_vec_linfa_const_only_theory(ys: Vec<bool>) -> f64 {
    let pos_count = ys.iter().filter(|x| **x).count();
    let neg_count = ys.len() - pos_count;

    (pos_count as f64 / neg_count as f64).ln()
}

//fn is_label_flipped_linfa_labels<F, C: PartialOrd + std::fmt::Debug>(
//    labels: &mylinfa_logistic::BinaryClassLabels<F, C>,
//) -> bool {
//    log::debug!("labels {:?}", labels.pos.class);
//    log::debug!("labels {:?}", labels.neg.class);
//    // always +1.0
//    //log::debug!("labels {:?}", lr.labels().pos.label);
//    // always -1.0
//    //log::debug!("labels {:?}", lr.labels().neg.label);
//
//    // true if pos.class==false
//    !labels.pos.class
//}

/*
// linfa v0.7.0
fn logistic_regression_covs_vec_linfa(
    v: Vec<f64>, //1d
    row_n: usize,
    col_n: usize,
    ys: Vec<bool>,
    means: &Vec<f64>,
    stds: &Vec<f64>,
) -> (Vec<f64>, f64) {
    let mt: Array2<f64> = Array::from_shape_vec((row_n, col_n), v).unwrap();
    log::debug!("mt {}, {}, {}", mt[[0, 0]], mt[[0, 1]], mt[[0, 2]]);
    log::debug!(
        "mt {}, {}, {}",
        mt[[row_n - 1, 0]],
        mt[[row_n - 1, 1]],
        mt[[row_n - 1, 2]]
    );

    let y = Array1::from(ys.clone());
    //let ys = Array1::from(ys_bool);
    let d = Dataset::new(mt, y);
    //let d = Dataset::new(mt, ys).map_targets(|x| if *x { "d" } else { "n" });

    let cov_n = means.len();

    // try 1.0
    println!("init 0.0");
    let init = Array1::from(vec![0.0f64; cov_n + 1]);
    //let init=Array1::from(vec![0.0f64;cov_n+1]);

    //let max_iter = 10000;
    // default
    let max_iter = 100;
    // default
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
    log::debug!("labels {:?}", lr.labels());

    //let mut intercept=vec![0.0f64;1];
    //let mut intercept=Vec::with_capacity(1);
    //lr.intercept().copy_col_as_vec(0, &mut intercept);
    //let intercept = lr.intercept().clone().to_row_vector()[0];

    //let intercept=lr.intercept().col_iter().collect::<Vec<f64>>();
    //let intercept = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();

    let coefs_norm = lr.params().iter().map(|x| *x).collect::<Vec<f64>>();

    let intercept_norm = lr.intercept();
    log::debug!("intercept_norm {:?}", intercept_norm);
    if intercept_norm.is_nan() {
        panic!("Intercept is NaN.");
    }

    //let x = linfa_logistic::label_classes(y).unwrap();
    //log::debug!("x {:?}", x);

    /// NEVER use this!! use is_label_flipped_linfa_labels()
    /// linfa regard the label with more count as 1
    /// So, if y=true is less than y=false, flip the coefficient
    /// Should be consistent with [here](https://github.com/rust-ml/linfa/blob/4adece5fac27d22d15f218387fbcd7df990ae43d/algorithms/linfa-logistic/src/lib.rs#L270)
    fn is_label_flipped_linfa(ys: &[bool]) -> bool {
        let y0_n = ys.iter().filter(|x| !**x).count();
        let y1_n = ys.len() - y0_n;
        if y0_n == y1_n {
            panic!("y0_n==y1_n");
        }
        y0_n > y1_n
    }

    let coefs_norm = if is_label_flipped_linfa(&ys) {
        coefs_norm.iter().map(|x| -*x).collect::<Vec<f64>>()
    } else {
        coefs_norm
    };

    let intercept_norm = if is_label_flipped_linfa(&ys) {
        -intercept_norm
    } else {
        intercept_norm
    };

    // c=c' / std
    let coefs = coefs_norm
        .iter()
        .zip(stds.iter())
        .map(|(c, s)| c / s)
        .collect::<Vec<f64>>();

    // i = i' - sum( c' * m / std) = i' - sum( c * m )
    let intercept = intercept_norm
        + coefs
            .iter()
            .zip(means.iter())
            .map(|(c, m)| -c * m)
            .sum::<f64>();

    (coefs, intercept)
}
*/

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec;

    #[test]
    fn test_array2_shape() {
        // col=2, row=3
        let v = vec![1.0f64, 2.0, 3.0, 4.0, 5.0, 6.0];

        let col_n = 2;
        let row_n = 3;
        let mt: Array2<f64> = Array::from_shape_vec((row_n, col_n), v).unwrap();
        //let row_n = 2;
        //let col_n = 3;
        //let mt: Array2<f64> = Array::from_shape_vec((row_n, col_n), v).unwrap();
        assert_eq!(mt.shape(), &[row_n, col_n]);
        // [ [1, 2], [3, 4], [5, 6]]
        assert_eq!(mt[[0, 0]], 1.0f64);
        assert_eq!(mt[[0, 1]], 2.0f64);
        assert_eq!(mt[[1, 0]], 3.0f64);
        assert_eq!(mt[[1, 1]], 4.0f64);
        assert_eq!(mt[[2, 0]], 5.0f64);
        assert_eq!(mt[[2, 1]], 6.0f64);

        //assert_eq!(convert_vec2d_to_vec1_row_major(v), (v_ans, 2, 3));
    }

    #[test]
    fn test_array2_shape2() {
        // col=2, row=3
        let v = vec![1.0f64, 2.0, 3.0, 4.0, 5.0, 6.0];

        let col_n = 2;
        let row_n = 3;
        let mt: Array2<f64> = Array::from_shape_vec((col_n, row_n), v).unwrap();
        //let row_n = 2;
        //let col_n = 3;
        //let mt: Array2<f64> = Array::from_shape_vec((row_n, col_n), v).unwrap();
        assert_eq!(mt.shape(), &[col_n, row_n]);
        // [[1,2,3],[4,5,6]]
        assert_eq!(mt[[0, 0]], 1.0f64);
        assert_eq!(mt[[0, 1]], 2.0f64);
        assert_eq!(mt[[0, 2]], 3.0f64);
        assert_eq!(mt[[1, 0]], 4.0f64);
        assert_eq!(mt[[1, 1]], 5.0f64);
        assert_eq!(mt[[1, 2]], 6.0f64);
        //assert_eq!(convert_vec2d_to_vec1_row_major(v), (v_ans, 2, 3));
    }

    #[test]
    fn test_logistic_regression_covs_vec_linfa() {
        // col major
        let v = vec![vec![-1.0f64, -0.01, 0.01, 1.0]];
        let y = vec![false, true, false, true];
        let means = vec![0.0f64];
        let stds = vec::std_v2(&v);
        assert_eq!(stds.len(), 1);

        let v_norm = vec::convert_vec2d_to_row_major(&vec::norm_vec2d(v));
        let (v_norm_1d, (row_n, col_n)) = vec::convert_vec2d_to_vec1(v_norm);
        //let (v_norm_1d, row_n, col_n) = vec::convert_vec2d_to_vec1_row_major(v_norm);
        let (coefs, intercept) =
            logistic_regression_vec_linfa_normalize(v_norm_1d, col_n, row_n, y, &means, &stds);
        assert_eq!(coefs.len(), 1);
        assert_float_absolute_eq!(intercept, 0.0);

        // expected answer was calculated using https://stats.blue/Stats_Suite/logistic_regression_calculator.html
        assert_float_absolute_eq!(coefs[0], 5.2672, 1e-3);
    }

    #[test]
    #[should_panic]
    fn test_logistic_regression_covs_vec_linfa_2() {
        // col major
        let v = vec![vec![0.0f64, 0.0, 0.0, 0.0]];
        let y = vec![false, true, false, true];
        let means = vec![0.0f64];
        let stds = vec::std_v2(&v);
        assert_eq!(stds.len(), 1);

        let v_norm = vec::convert_vec2d_to_row_major(&vec::norm_vec2d(v));
        let (v_norm_1d, (row_n, col_n)) = vec::convert_vec2d_to_vec1(v_norm);
        //let (v_norm_1d, row_n, col_n) = vec::convert_vec2d_to_vec1_row_major(v_norm);
        // panic?
        //let (coefs, intercept) =
        let _ = logistic_regression_vec_linfa_normalize(v_norm_1d, col_n, row_n, y, &means, &stds);
        //assert_eq!(coefs.len(), 1);
        //asert_float_absolute_eq!(intercept, 0.0);
    }

    #[test]
    fn test_logistic_regression_vec_linfa_no_intercept_1col() {
        // col major
        let v = vec![5.0f64, 5.0, 5.0, 5.0, 5.0];
        let y = vec![false, true, false, true, false];

        let coef = logistic_regression_vec_linfa_no_intercept_1col(v, y.clone());

        let coef_theory = logistic_regression_vec_linfa_const_only_theory(y.clone());

        assert_float_absolute_eq!(5.0 * coef, coef_theory, 1e-3);
    }
}
