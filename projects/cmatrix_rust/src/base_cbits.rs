//! Base trait for CVec family
//!
//! To avoid name collision to base_cmatrix, use _v
//!

use super::iterator::BoolIter;
use super::CBitsRef;
use super::{BaseCVec, BaseCVecMut};

/// for CBits, CBitsRef, CBitsMut
pub trait BaseCBits: BaseCVec {
    #[inline]
    fn get_unchecked_b(&self, col_i: usize) -> bool {
        self.get_unchecked(0, col_i, 0)
    }

    #[inline]
    fn get_b(&self, col_i: usize) -> bool {
        self.get(0, col_i, 0)
    }

    #[inline]
    fn inner_vals_col_digit_b(&self) -> Vec<bool> {
        self.inner_vals_col_digit(0, 0)
    }

    fn count(&self) -> usize {
        self.iter().filter(|b| *b).count()
    }

    fn count_false(&self) -> usize {
        self.col_n() - self.count()
    }

    fn to_cbits_ref_b(&self) -> CBitsRef {
        CBitsRef::new(self.inner(), self.col_n())
    }

    fn iter(&self) -> BoolIter {
        BoolIter::new(self.to_cbits_ref_b())
    }
}

/// for CBits, CBitsMut
pub trait BaseCBitsMut: BaseCVecMut + BaseCBits {
    fn set_bool_unchecked_b(&mut self, b: bool, col_i: usize) {
        self.set_bool_unchecked(b, 0, col_i, 0);
    }

    fn set_bool_b(&mut self, b: bool, col_i: usize) {
        self.set_bool(b, 0, col_i, 0);
    }
}
