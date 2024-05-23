//! Base trait for CVec family
//!
//! To avoid name collision to base_cmatrix, use _v
//!

use super::calc;
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

    #[allow(unreachable_code)]
    fn count(&self) -> usize {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            return count_simd(self.inner(), self.col_n());
        }
        //self.count_nosimd()
        count_nosimd(self.inner(), self.col_n())

        //self.iter().filter(|b| *b).count()
        // how to do this?
        // -> use self.count_nosimd()
        //return count(self);
    }

    // TODO use cmatrix::calc::count_nosimd() is faster
    //fn count_nosimd(&self) -> usize {
    //    self.iter().filter(|b| *b).count()
    //}

    fn count_false(&self) -> usize {
        self.col_n() - self.count()
    }

    /*     // TODO: integrate to count()
    fn count_simd(&self) -> usize {
        count_simd(self.inner(), self.col_n())
        /*         fn pop(x: u32) -> usize {
            unsafe { core::arch::x86_64::_popcnt64(x as i64) as usize }
        }

        let n = self.col_n();
        let ys = self.inner();

        let mut count = 0;

        // padding should be false
        for ni in 0..(n / 32 + 1) {
            let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
            count = count + pop(ys_b32);
        }

        count */
    } */

    fn as_cbits_ref_b(&self) -> CBitsRef {
        CBitsRef::new(self.inner(), self.col_n())
    }

    fn iter(&self) -> BoolIter {
        BoolIter::new(self.as_cbits_ref_b())
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn count_simd(ys: &[u8], n: usize) -> usize {
    let mut count = 0;

    // padding should be false
    for ni in 0..(n / 32 + 1) {
        let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        count = count + calc::popcnt_simd(ys_b32);
    }

    count
}

fn count_nosimd(ys: &[u8], n: usize) -> usize {
    let mut count = 0;

    // padding should be false
    for ni in 0..(n / 32 + 1) {
        let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        count = count + calc::popcnt_nosimd(ys_b32);
    }

    count
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
