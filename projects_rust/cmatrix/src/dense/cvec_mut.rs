//! Compressed Vector Mutable Reference
//!
//!
//!

use super::CMatrixMut;
use super::B8;
use super::{BaseCMatrix, BaseCMatrixMut, BaseCVec, BaseCVecMut};

#[derive(Debug, PartialEq, Eq)]
pub struct CVecMut<'a>(CMatrixMut<'a>);

impl<'a> BaseCMatrix for CVecMut<'a> {
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

impl<'a> BaseCMatrixMut for CVecMut<'a> {
    fn inner_mut(&mut self) -> &mut [B8] {
        self.0.inner_mut()
    }
}

impl<'a> BaseCVec for CVecMut<'a> {}
impl<'a> BaseCVecMut for CVecMut<'a> {}

impl<'a> CVecMut<'a> {
    pub fn new(inner: &'a mut [B8], col_n: usize, digit_n: usize) -> Self {
        let m = CVecMut(CMatrixMut::new(inner, 1, col_n, digit_n));
        m.check_new();
        m
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
