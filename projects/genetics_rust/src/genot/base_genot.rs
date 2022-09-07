//! module for predictions-related func.
//! len: length of vector
//! size: sample size or snv size
//! predictions shape is m * 2 * n and different from 2 * m * n of C++ ver.
//! This does not require m when getting index of part of predictions
//!
//!
//! TODO: check using _unchecked or not.
//! ex. loading plink is safe to use _unchecked.
//!
//!

use super::genot_iterator::{GenotCountIter, GenotIter, GenotIterMut};
use super::{GenotMut, GenotSnvMut, GenotSnvRef};
use crate::samples::prelude::*;
use cmatrix::prelude::*;

/*
// 8 x 1 bit
pub type B8 = u8;
// 4 x 2 bit
pub type B8_2 = u8;
// si=0 : dominant
// si=1 : recessive
pub const THRESHOLD_SNV: [f64; 2] = [0.5, 1.5];

// count to pred
// 0 => (0,0)
// 1 => (1,0)
// 2 => (1,1)
// 3 => (0,1)
const COUNT_PRED0_AR: [bool; 4] = [false, true, true, false];
const COUNT_PRED1_AR: [bool; 4] = [false, false, true, true];

// plink bed code to pred
// 00(2); 1,1
// 01(NA); 0,1
// 10(1); 1,0
// 11(0); 0,0
const BED_PRED0_AR: [bool; 4] = [true, false, true, false];
const BED_PRED1_AR: [bool; 4] = [true, true, false, false];
 */

// TODO: write down where to implement Inner
pub trait BaseGenot {
    type Inner: BaseCMatrix;
    fn genot_inner(&self) -> &Self::Inner;

    fn m(&self) -> usize {
        self.genot_inner().row_n()
    }
    fn n(&self) -> usize {
        self.genot_inner().col_n()
    }
    fn len_n(&self) -> usize {
        self.genot_inner().col_len()
    }
    fn len(&self) -> usize {
        self.genot_inner().len()
    }
    fn get_unchecked(&self, mi: usize, ni: usize) -> u8 {
        self.genot_inner().get_val_unchecked(mi, ni)
    }
    fn get(&self, mi: usize, ni: usize) -> u8 {
        self.genot_inner().get_val(mi, ni)
    }
    fn predict_snv_s(&self, mi: usize, si: usize) -> &[B8] {
        self.genot_inner().inner_col_digit(mi, si)
    }

    fn vals_snv(&self, mi: usize) -> Vec<u8> {
        self.genot_inner().inner_vals_col(mi)
    }
    fn vals_snv_s(&self, mi: usize, si: usize) -> Vec<bool> {
        self.genot_inner().inner_vals_col_digit(mi, si)
    }
    fn vals_snv_s_u8(&self, mi: usize, si: usize) -> Vec<u8> {
        self.genot_inner().inner_vals_col_digit_u8(mi, si)
    }

    fn iter_snv(&self) -> GenotIter {
        GenotIter::new(self.genot_inner().iter_row())
    }

    fn to_genot_snv(&self, mi: usize) -> GenotSnvRef {
        GenotSnvRef::new(self.genot_inner().to_cvec_ref(mi))
    }
}

pub trait BaseGenotMut: BaseGenot
where
    Self::Inner: BaseCMatrixMut,
{
    //type Inner: BaseCMatrixMut;
    fn genot_inner_mut(&mut self) -> &mut Self::Inner;
    fn iter_snv_mut(&mut self) -> GenotIterMut {
        GenotIterMut::new(self.genot_inner_mut().iter_row_mut())
    }
    fn as_genot_mut(&mut self) -> GenotMut {
        GenotMut::new(self.genot_inner_mut().to_cmatrix_mut())
    }
    fn as_genot_snv_mut(&mut self, mi: usize) -> GenotSnvMut {
        GenotSnvMut::new(self.genot_inner_mut().to_cvec_mut(mi))
    }
    fn as_genot_snvs_mut(&mut self, m_begin: usize, m_end: usize) -> GenotMut {
        GenotMut::new(self.genot_inner_mut().rows_mut(m_begin, m_end))
    }
}

pub trait BaseGenotSnv: BaseGenot
where
    Self::Inner: BaseCVec,
{
    //type Inner: BaseCVec;
    //fn genot_inner_snv(&self) -> &Self::Inner;

    // Which to call?
    // use CVec:: or GenotTwin::
    // -> should be CVec since cvec is inside of struct
    // also, boundry is BaseCVec
    // also, no need to think about row_i
    fn vals(&self) -> Vec<u8> {
        //self.counts_snv(0)
        self.genot_inner().inner_vals_col_v()
    }
    fn get_val(&self, ni: usize) -> u8 {
        self.genot_inner().get_val_unchecked_v(ni)
    }
    fn get_val_unchecked(&self, ni: usize) -> u8 {
        self.genot_inner().get_val_unchecked_v(ni)
    }
    fn predict_s(&self, si: usize) -> &[B8] {
        self.genot_inner().inner_col_digit_v(si)
    }
    fn iter(&self) -> GenotCountIter {
        GenotCountIter::new(self.genot_inner().iter_scal())
    }

    // counts are
    // (d2, n2, d1, n1, d0, n0, dm, nm)
    fn stat_contingency_table(
        &self,
        phe: &Phe,
    ) -> (usize, usize, usize, usize, usize, usize, usize, usize) {
        let mut d2_ = 0;
        let mut n2_ = 0;
        let mut d1_ = 0;
        let mut n1_ = 0;
        let mut d0_ = 0;
        let mut n0_ = 0;
        let mut dm_ = 0;
        let mut nm_ = 0;

        self.iter().zip(phe.iter()).for_each(|(count, y)| {
            if y {
                match count {
                    2 => d2_ += 1,
                    1 => d1_ += 1,
                    0 => d0_ += 1,
                    3 => dm_ += 1,
                    _ => panic!("Wrong count."),
                }
            } else {
                match count {
                    2 => n2_ += 1,
                    1 => n1_ += 1,
                    0 => n0_ += 1,
                    3 => nm_ += 1,
                    _ => panic!("Wrong count."),
                }
            }
        });
        (d2_, n2_, d1_, n1_, d0_, n0_, dm_, nm_)
    }
}

pub trait BaseGenotSnvMut: BaseGenotMut + BaseGenotSnv
where
    Self::Inner: BaseCVecMut,
{
    //type Inner: BaseCVecMut;
    //fn genot_inner_mut_snv(&mut self) -> &mut Self::Inner;
    #[inline]
    fn set(&mut self, val: u8, ni: usize) {
        self.genot_inner_mut().set_v(val, ni);
    }
    #[inline]
    fn set_unchecked(&mut self, val: u8, ni: usize) {
        self.genot_inner_mut().set_unchecked_v(val, ni);
    }

    // TODO: direct from bed code
    /// for plink bed
    /// 00(2); 1,1
    /// 10(1); 1,0
    /// 11(0); 0,0
    /// 01(NA); 0,1
    #[inline]
    fn set_bed_code_unchecked(&mut self, bcode: u8, ni: usize) {
        let b0 = (bcode & 1) == 0;
        let b1 = (bcode & 2) == 0;
        self.genot_inner_mut().set_by_bool_unchecked_v(b0, b1, ni);
    }

    fn as_genot_snv_mut_snv(&mut self) -> GenotSnvMut {
        GenotSnvMut::new(self.genot_inner_mut().to_cvec_mut_v())
    }
}
