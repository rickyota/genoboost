// mysmartcore v0.3.1
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

//smartcore v3 or mysmartcore
//pub fn logistic_regression_covs(samples: &Samples) -> Vec<Wgt> {
//    if let Some(covs) = samples.covs() {
//        // load normed vals
//        let v = covs.vals_row_major_norm();
//        let ys = samples.phe_unwrap().inner_i32();
//
//        let cov_n = covs.covs_n();
//        //let cov_n=v.len();
//        log::debug!("cov_n {}", cov_n);
//        log::debug!("ys true count {}", samples.phe_unwrap().count());
//        log::debug!("ys false count {}", samples.phe_unwrap().count_false());
//        log::debug!("ys true count {}", ys.iter().filter(|x| **x).count());
//        log::debug!("ys false count {}", ys.iter().filter(|x| !**x).count());
//
//        log::debug!("norm v");
//        log::debug!("norm mt {}, {}", v[0][0], v[0][1]);
//
//        let means = samples.covs().unwrap().means();
//        let stds = samples.covs().unwrap().stds();
//        log::debug!("means {:?}", means);
//        log::debug!("stds {:?}", stds);
//
//        let (coefs, intercept) = logistic_regression_covs_vec(&v, &ys, &means, &stds);
//
//        assert_eq!(coefs.len(), covs.covs_n());
//
//        let mut wgts_cov: Vec<Wgt> = Vec::new();
//
//        //let mut intercept=vec![0.0f64;1];
//        //let mut intercept=Vec::with_capacity(1);
//        //lr.intercept().copy_col_as_vec(0, &mut intercept);
//        //let intercept = lr.intercept().clone().to_row_vector()[0];
//
//        //let intercept=lr.intercept().col_iter().collect::<Vec<f64>>();
//        //let intercept = lr.intercept().iter().map(|x| *x).collect::<Vec<f64>>();
//        //log::debug!("intercept {:?}", intercept);
//        //assert_eq!(intercept.len(), 1);
//        //let intercept = intercept[0];
//        if intercept.is_nan() {
//            panic!("Intercept is NaN.");
//        }
//        let cov_id = CovId::new_const();
//        let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));
//        wgts_cov.push(wgt_const);
//
//        let cov_indexs = covs.cov_indexs().unwrap();
//        for wgt_i in 0..covs.covs_n() {
//            let coef = coefs[wgt_i];
//            if coef.is_nan() {
//                panic!("coef is NaN.");
//            }
//            let cov_id = CovId::new_cov(cov_indexs[wgt_i].name().to_owned());
//            let wgt_cov = Wgt::construct_cov(cov_id, Coef::Linear(coef));
//            wgts_cov.push(wgt_cov);
//        }
//
//        wgts_cov
//    } else {
//        // regression on const
//
//        log::debug!("Regression on const value only.");
//
//        log::debug!("ys true count {}", samples.phe_unwrap().count());
//        log::debug!("ys false count {}", samples.phe_unwrap().count_false());
//
//        let y1_n = samples.phe_unwrap().count();
//        let y0_n = samples.phe_unwrap().count_false();
//
//        let intercept = logistic_regression_const(y0_n, y1_n);
//
//        let mut wgts_cov: Vec<Wgt> = Vec::new();
//        if intercept.is_nan() {
//            panic!("Intercept is NaN.");
//        }
//        let cov_id = CovId::new_const();
//        let wgt_const = Wgt::construct_cov(cov_id, Coef::Linear(intercept));
//        wgts_cov.push(wgt_const);
//
//        wgts_cov
//    }
//}

//smartcore v3 or mysmartcore
pub fn logistic_regression_covs_samples(samples: &Samples) -> (Vec<f64>, f64) {
    let covs = samples.covs().unwrap();

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

    let means = samples.covs().unwrap().means();
    let stds = samples.covs().unwrap().stds();
    log::debug!("means {:?}", means);
    log::debug!("stds {:?}", stds);

    logistic_regression_covs_vec(&v, &ys, &means, &stds)
}

fn logistic_regression_covs_vec(
    v: &Vec<Vec<f64>>,
    ys: &Vec<i32>,
    means: &Vec<f64>,
    stds: &Vec<f64>,
) -> (Vec<f64>, f64) {
    // dense_matrix from 2d array
    // https://docs.rs/smartcore/latest/src/smartcore/linalg/naive/dense_matrix.rs.html#248
    //let (v, row_n, col_n) = covs_val.vals_row_major_vec();
    //log::debug!("row, col {}, {}",row_n,col_n);
    let mt = DenseMatrix::from_2d_vec(&v);
    //let mt = DenseMatrix::from_2d_vec(&covs_val.vals_row_major());
    // default: alpha=0.0
    // confirmed to output the same result as sklearn
    let lr = LogisticRegression::fit(&mt, ys, Default::default()).unwrap();
    //let lr = LogisticRegression::fit(&mt, &ys, Default::default()).unwrap();
    log::debug!("logreg");
    log::debug!("intercept {:?}", lr.intercept());
    log::debug!("coef {:?}", lr.coefficients());

    //revert normalization
    //let means = samples.covs().unwrap().means();
    //let stds = samples.covs().unwrap().stds();

    let coefs_norm = lr.coefficients().iter().map(|x| *x).collect::<Vec<f64>>();
    //assert_eq!(coefs_norm.len(), covs.covs_n());

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

    (coefs, intercept)
}

fn logistic_regression_const(y0_n: usize, y1_n: usize) -> f64 {
    // regression on const
    let y1_n = y1_n as f64;
    let y0_n = y0_n as f64;

    // TODO: test
    let intercept = (y1_n / y0_n).ln();
    intercept
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec;

    #[test]
    fn test_logistic_regression_covs_vec() {
        // TOFIX: prepare larger data for stability
        let v = vec![vec![-1.0f64], vec![-0.01], vec![0.01], vec![1.0]];
        let y = vec![0i32, 1, 0, 1];
        let means = vec![0.0f64];
        let stds = vec::std_v2(&vec::convert_vec2_to_col_major(&v));
        let (coefs, intercept) = logistic_regression_covs_vec(&v, &y, &means, &stds);
        assert_eq!(coefs.len(), 1);
        assert_float_absolute_eq!(intercept, 0.0, 1e-2);
        // fail
        //assert_float_absolute_eq!(intercept, 0.0);

        // expected answer was calculated using https://stats.blue/Stats_Suite/logistic_regression_calculator.html
        // fail
        //assert_float_absolute_eq!(coefs[0], 5.2672, 1e-3);
    }

    #[test]
    fn test_logisticregression_const() {
        let y1_n = 100;
        let y0_n = 10;
        let intercept = logistic_regression_const(y0_n, y1_n);
        assert_float_absolute_eq!(intercept, (100.0f64 / 10.0).ln());
    }
}
