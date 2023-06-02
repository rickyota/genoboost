// Dirty implementation for trait

use crate::mut_ref_in_struct::Predict;

mod impl_sub;

// Do not implement this Trait and use AsRef instead
pub trait PredictTrait {
    fn n(&self) -> usize;
    fn val(&self) -> &[u8];
}

// Corresponds to Vec, String
pub struct PredictVec {
    predict_v: Vec<u8>,
    n: usize,
}

// splitting impl is allowed
impl PredictVec {
    fn n(&self) -> usize {
        self.n
    }
}

impl PredictVec {
    fn val_len(&self) -> usize {
        self.predict_v.len()
    }
}

/*
// splitting impl trait is not allowed
impl PredictTrait for PredictVec {
    fn n(&self) -> usize {
        self.n
    }
}

impl PredictTrait for PredictVec {
    fn val(&self) -> &[u8] {
        &self.predict_v
    }
}
*/

pub fn test() {}
