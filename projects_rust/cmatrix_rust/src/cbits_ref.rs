//! Compressed Bits Reference
//!
//!
//!

use super::CMatrixRef;
use super::B8;
use super::{BaseCBits, BaseCMatrix, BaseCVec};

/// Compressed Bits
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct CBitsRef<'a>(CMatrixRef<'a>);

impl<'a> BaseCMatrix for CBitsRef<'a> {
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

impl<'a> BaseCVec for CBitsRef<'a> {}
impl<'a> BaseCBits for CBitsRef<'a> {}

impl<'a> CBitsRef<'a> {
    pub fn new(inner: &'a [B8], col_n: usize) -> Self {
        let m = CBitsRef(CMatrixRef::new(inner, 1, col_n, 1));
        m.check_new();
        m
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
