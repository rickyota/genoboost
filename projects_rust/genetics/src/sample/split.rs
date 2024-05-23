use rand::prelude::*;

use mysmartcore::linalg::basic::arrays::Array2;
// mysmartcore v0.3.1
// https://smartcorelib.org/user_guide/model_selection.html
use mysmartcore::linalg::basic::matrix::DenseMatrix;
use mysmartcore::model_selection::KFold;
use mysmartcore::model_selection::*;

/// output index is sorted
pub fn split_samples_random(
    n: usize,
    cvn: usize,
    // for cvn=1
    prop_train: Option<f64>,
    seed: Option<u64>,
) -> Vec<(Vec<usize>, Vec<usize>)> {
    // create cv samples
    //let sample_idx_cvs: Vec<(Vec<usize>, Vec<usize>)> = if cvn == 1 {

    if cvn == 0 {
        panic!("cvn must be >0");
    }

    if cvn == 1 {
        let mut rng: rand::rngs::StdRng = match seed {
            Some(seed) => StdRng::seed_from_u64(seed),
            None => StdRng::from_entropy(),
        };

        //if let Some(seed) = seed {
        //    StdRng::seed_from_u64(seed);
        //}
        let mut vec: Vec<usize> = (0..n).collect();
        vec.shuffle(&mut rng);
        //vec.shuffle(&mut rand::thread_rng());

        let prop_train = prop_train.expect("prop_train must be specified if cvn=1");
        if prop_train < 0.0 || prop_train > 1.0 {
            panic!("prop_train must be between 0.0 and 1.0");
        }
        //let prop_tr = 0.8;
        let n_tr = (prop_train * (n as f64)) as usize;
        let mut tr = vec[..n_tr].to_vec();
        let mut va = vec[n_tr..].to_vec();
        tr.sort();
        va.sort();
        vec![(tr, va)]
        //vec![(vec[..n_tr].to_vec(), vec[n_tr..].to_vec())]
    } else {
        // unnecessary panic
        //if prop_train.is_some() {
        //    panic!("prop_train must be None if cvn>1");
        //}

        let sample_matrix: DenseMatrix<f64> = DenseMatrix::zeros(n, 1);
        KFold::default()
            .with_n_splits(cvn)
            .with_shuffle(true)
            .with_seed(seed)
            .split(&sample_matrix)
            .into_iter()
            .collect()
    }
}

pub fn extract_sample_cvi(sample_idx_cvi: &[usize], samples_in: &[String]) -> Vec<String> {
    sample_idx_cvi
        .iter()
        .map(|idx| samples_in[*idx].clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vec;

    #[test]
    fn test_split_samples_random_cvn_1() {
        let n_in = 10;
        let cvn = 1;
        let seed = Some(0);
        let prop_train = Some(0.8);
        let sample_idx_cvs = split_samples_random(n_in, cvn, prop_train, seed);
        assert_eq!(sample_idx_cvs.len(), cvn);
        assert_eq!(sample_idx_cvs[0].0.len(), 8);
        assert_eq!(sample_idx_cvs[0].1.len(), 2);
        assert!(vec::is_sorted(&sample_idx_cvs[0].0));
        assert!(vec::is_sorted(&sample_idx_cvs[0].1));
    }

    #[test]
    fn test_split_samples_random() {
        let n_in = 10;
        let cvn = 2;
        let seed = Some(0);
        //let prop_train = Some(0.8);
        let mut sample_idx_cvs = split_samples_random(n_in, cvn, None, seed);

        assert_eq!(sample_idx_cvs.len(), cvn);
        for cvi in 0..cvn {
            assert_eq!(sample_idx_cvs[cvi].0.len(), 5);
            assert_eq!(sample_idx_cvs[cvi].1.len(), 5);
            assert!(vec::is_sorted(&sample_idx_cvs[cvi].0));
            assert!(vec::is_sorted(&sample_idx_cvs[cvi].1));
        }

        assert_eq!(sample_idx_cvs[0].0, sample_idx_cvs[1].1);
        // cvi=0 and 1 should be complement
        assert_eq!(sample_idx_cvs[0].0.sort(), sample_idx_cvs[1].1.sort());
        assert_eq!(sample_idx_cvs[0].1.sort(), sample_idx_cvs[0].1.sort());
    }

    #[test]
    fn test_split_samples_fixseed() {
        let n_in = 10;
        let cvn = 2;
        let seed = Some(0);
        let sample_idx_cvs = split_samples_random(n_in, cvn, None, seed);
        let sample_idx_cvs_2 = split_samples_random(n_in, cvn, None, seed);

        assert_eq!(sample_idx_cvs, sample_idx_cvs_2);
    }
}
