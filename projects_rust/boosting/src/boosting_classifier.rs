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
        log::debug!("fit");
        let n = dataset.samples_n();
        let m = dataset.snvs_n();
        let mut scores = vec![0.0f64; n];

        let mut ws = vec![0.0f64; n];

        // len != n only for ps
        let mut ps_pad: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
        ps_pad.resize(n + 32, 0.0);

        log::debug!("fit done");
    }
    fn score<D: DatasetTrait>(&self, fout: &str, dataset: D, iterations: &[usize]) {}
}

#[cfg(test)]
mod tests {
    use super::*;
}
