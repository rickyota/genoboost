//! Compressed Bits Reference
//!
//!
//!

use super::CMatrixMut;
use super::B8;
use super::{BaseCBits, BaseCBitsMut, BaseCMatrix, BaseCMatrixMut, BaseCVec, BaseCVecMut};

/// Compressed Vector
#[derive(Debug, PartialEq, Eq)]
pub struct CBitsMut<'a>(CMatrixMut<'a>);

impl<'a> BaseCMatrix for CBitsMut<'a> {
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

impl<'a> BaseCMatrixMut for CBitsMut<'a> {
    #[inline]
    fn inner_mut(&mut self) -> &mut [B8] {
        self.0.inner_mut()
    }
}

impl<'a> BaseCVec for CBitsMut<'a> {}
impl<'a> BaseCVecMut for CBitsMut<'a> {}
impl<'a> BaseCBits for CBitsMut<'a> {}
impl<'a> BaseCBitsMut for CBitsMut<'a> {}

impl<'a> CBitsMut<'a> {
    pub fn new(inner: &'a mut [B8], col_n: usize) -> Self {
        let m = CBitsMut(CMatrixMut::new(inner, 1, col_n, 1));
        m.check_new();
        m
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
