//! Compressed Vector
//!
//!
//!
//!
// Maybe, cvec and cbits are unnecessary with cvecref etc.

use super::CMatrix;
use super::B8;
use super::{BaseCMatrix, BaseCMatrixMut, BaseCVec, BaseCVecMut};

/// Compressed Vector
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CVec(CMatrix);

impl BaseCMatrix for CVec {
    #[inline]
    fn inner(&self) -> &[B8] {
        self.0.inner()
    }
    #[inline]
    fn col_n(&self) -> usize {
        self.0.col_n()
    }
    #[inline]
    fn row_n(&self) -> usize {
        self.0.row_n()
    }
    #[inline]
    fn digit_n(&self) -> usize {
        self.0.digit_n()
    }
}

impl BaseCMatrixMut for CVec {
    fn inner_mut(&mut self) -> &mut [B8] {
        self.0.inner_mut()
    }
}

impl BaseCVec for CVec {}
impl BaseCVecMut for CVec {}

impl CVec {
    //pub fn new(digit: usize, vec: Vec<B8>) -> CVec {
    pub fn new(digit_n: usize, vec: &[B8]) -> CVec {
        CVec(CMatrix::new(1, vec.len(), digit_n, vec))
    }
    pub fn new_empty(col_n: usize, digit: usize) -> CVec {
        CVec(CMatrix::new_zeros(1, col_n, digit))
    }

    pub fn inner_mut_vector(&mut self) -> &mut Vec<B8> {
        self.0.inner_mut_vector()
    }

    /*
    pub fn new_vec_bool(vec: Vec<bool>) -> CVec {
        let col_n = vec.len();
        CVec(CMatrix::new_vec_bool(1, col_n, vec))
    }
     */
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let vec = vec![3, 2, 1, 0];
        CVec::new(2, &vec);
    }

    #[test]
    fn test_get() {
        let vec = vec![3, 2, 1, 0];
        let v = CVec::new(2, &vec);
        assert_eq!(v.get_val_unchecked_v(0), 3);
        assert_eq!(v.get_val_unchecked_v(1), 2);
        assert_eq!(v.get_val_unchecked_v(2), 1);
        assert_eq!(v.get_val_unchecked_v(3), 0);
    }
}
