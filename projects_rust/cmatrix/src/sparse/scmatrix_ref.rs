//! Compressed Matrix Reference
//!

use super::BaseCMatrix;
use super::B8;

/// Compressed Matrix
/// Mainly used for genotype and phenotype.
///
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct SCMatrixRef<'a> {
    inner: &'a [B8],
    row_n: usize,
    col_n: usize,
    digit_n: usize,
}

impl<'a> BaseCMatrix for CMatrixRef<'a> {
    #[inline]
    fn inner(&self) -> &[B8] {
        &self.inner
    }
    #[inline]
    fn col_n(&self) -> usize {
        self.col_n
    }
    #[inline]
    fn row_n(&self) -> usize {
        self.row_n
    }
    #[inline]
    fn digit_n(&self) -> usize {
        self.digit_n
    }
}

impl<'a> CMatrixRef<'a> {
    pub(super) fn new(inner: &'a [B8], row_n: usize, col_n: usize, digit_n: usize) -> Self {
        let m = Self {
            inner,
            row_n,
            col_n,
            digit_n,
        };
        m.check_new();
        m
    }
}

/*
// no time to use
// return scalar
pub struct ScalIter<'a> {
    pub(super) cmatrix: CMatrixRef<'a>,
    pub(super) index_iter: ScalIndexIter,
}

// iterating over scalar
impl<'a> Iterator for ScalIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(ar) = self.index_iter.next() {
            return Some(self.cmatrix.get_val_ar(ar));
        }
        None
    }
}

pub struct BoolIter<'a> {
    cmatrix: CMatrixRef<'a>,
    index_iter: BoolIndexIter,
}

// iterating over bool
impl<'a> Iterator for BoolIter<'a> {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(ar) = self.index_iter.next() {
            return Some(self.cmatrix.get_ar(ar));
        }
        None
    }
}
 */

#[cfg(test)]
mod tests {
    //use super::*;

    /*     #[test]
    fn test_scal_iter() {
        let vec = vec![0, 1, 2, 1, 2, 3];
        let mut m = CMatrix::new_vec(2, 3, 2, vec);

        let mut iter = m.iter_scal();

        assert_eq!(iter.next().unwrap(), 0);
        assert_eq!(iter.next().unwrap(), 1);
        assert_eq!(iter.next().unwrap(), 2);
        assert_eq!(iter.next().unwrap(), 1);
        assert_eq!(iter.next().unwrap(), 2);
        assert_eq!(iter.next().unwrap(), 3);
        assert_eq!(iter.next(), None);
    } */
}
