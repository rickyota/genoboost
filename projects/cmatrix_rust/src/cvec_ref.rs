//! Compressed Vector Reference
//! One row of CMatrix
//!

use super::CMatrixRef;
use super::B8;
use super::{BaseCMatrix, BaseCVec};

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct CVecRef<'a>(CMatrixRef<'a>);

impl<'a> BaseCMatrix for CVecRef<'a> {
    #[inline]
    fn inner(&self) -> &[B8] {
        &self.0.inner()
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

impl<'a> BaseCVec for CVecRef<'a> {}

impl<'a> CVecRef<'a> {
    pub fn new(inner: &'a [B8], col_n: usize, digit_n: usize) -> Self {
        let m = CVecRef(CMatrixRef::new(inner, 1, col_n, digit_n));
        m.check_new();
        m
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
