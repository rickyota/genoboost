//! Compressed Bits
//! One row, one digit of CMatrix
//!
//!
//!
//!

use super::CMatrix;
use super::B8;
use super::{BaseCBits, BaseCBitsMut, BaseCMatrix, BaseCMatrixMut, BaseCVec, BaseCVecMut};

/// Compressed Vector
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CBits(CMatrix);

impl BaseCMatrix for CBits {
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

impl BaseCMatrixMut for CBits {
    #[inline]
    fn inner_mut(&mut self) -> &mut [B8] {
        self.0.inner_mut()
    }

    //#[inline]
    //fn inner_mut_vector(&mut self) -> &mut Vec<B8> {
    //    self.0.inner_mut_vector()
    //}
}

impl BaseCVec for CBits {}
impl BaseCVecMut for CBits {}

impl BaseCBits for CBits {}
impl BaseCBitsMut for CBits {}

impl CBits {
    pub fn new(vec: &[bool]) -> CBits {
        let col_n = vec.len();
        CBits(CMatrix::new_vec_bool(1, col_n, vec))
    }
    pub fn new_empty(col_n: usize) -> CBits {
        CBits(CMatrix::new_zeros(1, col_n, 1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let vec = vec![true, false];
        CBits::new(&vec);
    }

    #[test]
    fn test_get() {
        let vec = vec![true, false];
        let v = CBits::new(&vec);
        assert_eq!(v.get_unchecked(0, 0, 0), true);
        assert_eq!(v.get_unchecked(0, 1, 0), false);
    }

    #[test]
    fn test_count_true() {
        let y = vec![
            true, false, true, false, true, false, true, false, true, false,
        ];
        let v = CBits::new(&y);
        assert_eq!(v.count(), 5);
    }

    // test count_simd()
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[test]
    fn test_count_simd() {
        let y = vec![
            true, false, true, false, true, false, true, false, true, false,
        ];
        let v = CBits::new(&y);
        assert_eq!(v.count(), 5);
    }

    #[test]
    fn test_count_nosimd() {
        let y = vec![
            true, false, true, false, true, false, true, false, true, false,
        ];
        let v = CBits::new(&y);
        assert_eq!(v.count(), 5);
    }
}
