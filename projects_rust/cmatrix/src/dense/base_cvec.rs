//! Base trait for CVec family
//!
//! To avoid name collision to base_cmatrix, use _v
//!

use super::base_cmatrix::{BaseCMatrix, BaseCMatrixMut};
use super::iterator::ScalIter;
use super::B8;
use super::{CVecMut, CVecRef};

// only implement fn directly called by outsize
// other fn use fn in BaseCMatrix
// row_n=1 is assured since BaseCVec is bound for CVec and CBits
/// for CVec, CVecRef, CVecMut
pub trait BaseCVec: BaseCMatrix {
    #[inline]
    fn get_unchecked_v(&self, col_i: usize, digit_i: usize) -> bool {
        self.get_unchecked(0, col_i, digit_i)
    }

    #[inline]
    fn get_v(&self, col_i: usize, digit_i: usize) -> bool {
        self.get(0, col_i, digit_i)
    }

    #[inline]
    fn get_val_unchecked_v(&self, col_i: usize) -> u8 {
        self.get_val_unchecked(0, col_i)
    }

    #[inline]
    fn get_val_v(&self, col_i: usize) -> u8 {
        self.get_val(0, col_i)
    }

    fn inner_col_digit_v(&self, digit_i: usize) -> &[B8] {
        self.inner_col_digit(0, digit_i)
    }

    fn inner_vals_col_digit_v(&self, digit_i: usize) -> Vec<bool> {
        self.inner_vals_col_digit(0, digit_i)
    }

    fn inner_vals_col_v(&self) -> Vec<u8> {
        self.inner_vals_col(0)
    }

    fn as_cvec_ref_v(&self) -> CVecRef {
        CVecRef::new(self.inner(), self.col_n(), self.digit_n())
    }

    fn iter_scal(&self) -> ScalIter {
        ScalIter::new(self.as_cvec_ref_v())
    }
}

/// for CVec, CVecMut
pub trait BaseCVecMut: BaseCMatrixMut + BaseCVec {
    #[inline]
    fn set_v(&mut self, val: u8, col_i: usize) {
        self.set(val, 0, col_i)
    }
    #[inline]
    fn set_unchecked_v(&mut self, val: u8, col_i: usize) {
        self.set_unchecked(val, 0, col_i)
    }

    #[inline]
    fn set_by_bool_unchecked_v(&mut self, b0: bool, b1: bool, col_i: usize) {
        self.set_by_bool_unchecked(b0, b1, 0, col_i)
    }

    #[inline]
    fn set_by_bool_init_unchecked_v(&mut self, b0: bool, b1: bool, col_i: usize) {
        self.set_by_bool_init_unchecked(b0, b1, 0, col_i)
    }

    fn set_col_unchecked_v(&mut self, vec: &[u8]) {
        // TODO: set_vec_unchecked(vec, row_i=0)
        self.set_vec_unchecked(vec);
    }

    fn as_cvec_mut_v(&mut self) -> CVecMut {
        let (col_n, digit_n) = (self.col_n(), self.digit_n());
        CVecMut::new(self.inner_mut(), col_n, digit_n)
    }


    /*
    // TODO: check val<4
    /// assume row=1
    fn set_col_unchecked_v(&mut self, vec: &[u8]) {
        if self.row_n() != 1 {
            panic!("");
        }
        if self.digit_n() != 2 {
            panic!("");
        }
        vec.iter()
            .enumerate()
            .for_each(|(ci, x)| self.set_v(*x, ci));
        /*
        for ci in 0..self.col_n() {
            self.set(vec[ci], ci);
            /*
            let count = vec[ci];
            self.set_bool((count >> 1) != (count & 1), ci, 0);
            self.set_bool((count & 2) == 2, ci, 1);
            */
        }
         */
    }
     */

    /*     fn set_col_bool_v(&mut self, vec: &[bool]) {
        if self.row_n() != 1 {
            panic!("");
        }
        if self.digit_n() != 1 {
            panic!("");
        }
        for ci in 0..self.col_n() {
            let count = vec[ci];
            self.set_bool_unchecked(count, 0, ci, 0);
            //self.set_val((count >> 1) != (count & 1), ci, 0);
            //self.set_val((count & 2) == 2, ci, 1);
        }
    } */
}
