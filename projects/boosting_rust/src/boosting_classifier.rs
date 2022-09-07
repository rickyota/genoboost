/// May not use but this is an option
//  use crate::boosting::sample_weight;
use genetics_rust::{alloc, genot_bi_new, vec, B8};

pub trait Boosting {
    //fn fit(&self, &Dataset);
    fn fit<D: DatasetTrait>(&self, fout: &str, dataset: D);
    fn score<D: DatasetTrait>(&self, fout: &str, dataset: D, iterations: &[usize]);
}

pub struct AdaboostClassifier {
    iter: usize,
}

impl AdaboostClassifier {
    fn new(iteration: usize) -> Self {
        AdaboostClassifier { iteration }
    }
}

pub struct BaseClassifier {
    iter: usize,
}

pub struct ModelfreeClassifier {
    iter: usize,
    // default: 1.0
    learning_rate: f64, //learning_rate: Option<f64>
                        //base: BaseClassifier,
}

impl ModelfreeClassifier {
    fn new(iteration: usize, learning_rate: f64) -> Self {
        ModelfreeClassifier {
            iteration,
            learning_rate,
            //base: BaseClassifier { iter },
        }
    }
}

impl Boosting for ModelfreeClassifier {
    fn fit<D: DatasetTrait>(&self, fout: &str, dataset: D) {
        println!("fit");
        let n = dataset.samples_n();
        let m = dataset.snvs_n();
        let mut scores = vec![0.0f64; n];

        let mut ws = vec![0.0f64; n];

        // len != n only for ps
        let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        ps.resize(n + 32, 0.0);

        //let len_n = genot_bi_new::len_n(n);
        //let len_n = mistake::get_len_n(n);

        //let mut pred: Vec<B8> = Vec::with_capacity(len_n);

        //let mut losss_gt: Vec<f64> = Vec::with_capacity(m);

        //sample_weight::renew_ps(&mut ps, &ws, n);

        println!("fit done");
    }
    fn score<D: DatasetTrait>(&self, fout: &str, dataset: D, iterations: &[usize]) {}
}

#[cfg(test)]
mod tests {
    use super::*;
}
