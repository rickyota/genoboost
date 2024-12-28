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
use cmatrix;
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

    fn vals_snvs(&self) -> Vec<Vec<u8>> {
        self.genot_inner().inner_vals_u8()
        //self.iter_snv().map(|genot_snv| genot_snv.vals()).collect()
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
        GenotSnvRef::new(self.genot_inner().as_cvec_ref(mi))
    }
}

// BaseGenotRef trait is unnecessary
pub trait BaseGenotMut: BaseGenot
where
    Self::Inner: BaseCMatrixMut,
{
    fn genot_inner_mut(&mut self) -> &mut Self::Inner;
    fn iter_snv_mut(&mut self) -> GenotIterMut {
        GenotIterMut::new(self.genot_inner_mut().iter_row_mut())
    }
    fn as_genot_mut(&mut self) -> GenotMut {
        GenotMut::new(self.genot_inner_mut().as_cmatrix_mut())
    }
    fn as_genot_snv_mut(&mut self, mi: usize) -> GenotSnvMut {
        GenotSnvMut::new(self.genot_inner_mut().as_cvec_mut(mi))
    }
    fn as_genot_snvs_mut(&mut self, m_begin: usize, m_end: usize) -> GenotMut {
        GenotMut::new(self.genot_inner_mut().rows_mut(m_begin, m_end))
    }
    fn split_genot(&mut self, m: usize) -> (GenotMut, GenotMut) {
        let (cm1, cm2) = self.genot_inner_mut().split_rows(m);
        (GenotMut::new(cm1), GenotMut::new(cm2))
    }
}

pub trait BaseGenotSnv: BaseGenot
where
    Self::Inner: BaseCVec,
{
    //type Inner: BaseCVec;

    // not necessary: use inner: BaseCMatrix not BaseCVec
    //fn genot_inner_snv(&self) -> &Self::Inner;

    // Which to call?
    // use CVec:: or GenotTwin::
    // -> should be CVec since cvec is inside of struct
    // also, boundary is BaseCVec
    // also, no need to think about row_i
    fn vals(&self) -> Vec<u8> {
        self.genot_inner().inner_vals_col_v()
    }
    fn get_val(&self, ni: usize) -> u8 {
        self.genot_inner().get_val_unchecked_v(ni)
    }
    fn get_val_unchecked(&self, ni: usize) -> u8 {
        self.genot_inner().get_val_unchecked_v(ni)
    }
    // TODO: rename -> inheritance?
    fn predict_s(&self, si: usize) -> &[B8] {
        self.genot_inner().inner_col_digit_v(si)
    }
    fn iter(&self) -> GenotCountIter {
        GenotCountIter::new(self.genot_inner().iter_scal())
    }
    fn as_genot_snv(&self) -> GenotSnvRef {
        GenotSnvRef::new(self.genot_inner().as_cvec_ref_v())
    }

    fn count_missing(&self) -> usize {
        self.iter().filter(|&x| x == 3).count()
    }

    /// BE CAREFUL for order; different from contingency table
    fn count_table(&self) -> (usize, usize, usize, usize) {
        let n = self.n();
        let pred_s0m = self.predict_s(0);
        let pred_s1m = self.predict_s(1);

        let mut sums = (0usize, 0, 0);

        fn sum_byte_count(p0: u32, p1: u32) -> (usize, usize, usize) {
            // do not use d0/n0 since could be mixed with padding
            // TODO: calculate sum here? but need to deal with padding
            let c2 = p0 & p1;
            let c1 = p0 & (!p1);
            let cm = (!p0) & p1;

            (
                cmatrix::popcnt(c2),
                cmatrix::popcnt(c1),
                cmatrix::popcnt(cm),
            )
            //(pop(d2), pop(n2), pop(d1), pop(n1), pop(dm), pop(nm))
        }

        fn add_tuple3(
            sums: (usize, usize, usize),
            sums_: (usize, usize, usize),
        ) -> (usize, usize, usize) {
            (sums.0 + sums_.0, sums.1 + sums_.1, sums.2 + sums_.2)
        }

        for ni in 0..(n / 32 + 1) {
            let pred_s0_b32 =
                u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
            let pred_s1_b32 =
                u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
            //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());

            let sums_32 = sum_byte_count(pred_s0_b32, pred_s1_b32);

            sums = add_tuple3(sums, sums_32);
        }

        //println!("in table afr for: {} sec",  start_time.elapsed().as_micros());

        let c2 = sums.0;
        let c1 = sums.1;
        let cm = sums.2;

        let call = self.n();

        let c0 = call - (c2 + c1 + cm);

        (c0, c1, c2, cm)
    }

    /// Frequency of a2
    fn maf(&self) -> f64 {
        //use std::time::Instant;
        //let start_time = Instant::now();

        //println!("in table afr init: {} sec",  start_time.elapsed().as_micros());

        //let n = self.n();
        //let pred_s0m = self.predict_s(0);
        //let pred_s1m = self.predict_s(1);

        //let mut sums = (0usize, 0, 0);

        ////let mut sums: (usize, usize, usize, usize, usize, usize) = (0usize, 0, 0, 0, 0, 0);
        ////let mut sums=vec![0usize; 6];

        //fn sum_byte_count(p0: u32, p1: u32) -> (usize, usize, usize) {
        //    // do not use d0/n0 since could be mixed with padding
        //    // TODO: calculate sum here? but need to deal with padding
        //    let c2 = p0 & p1;
        //    let c1 = p0 & (!p1);
        //    let cm = (!p0) & p1;

        //    (
        //        cmatrix::popcnt(c2),
        //        cmatrix::popcnt(c1),
        //        cmatrix::popcnt(cm),
        //    )
        //    //(pop(d2), pop(n2), pop(d1), pop(n1), pop(dm), pop(nm))
        //}

        //fn add_tuple3(
        //    sums: (usize, usize, usize),
        //    sums_: (usize, usize, usize),
        //) -> (usize, usize, usize) {
        //    (sums.0 + sums_.0, sums.1 + sums_.1, sums.2 + sums_.2)
        //}

        //for ni in 0..(n / 32 + 1) {
        //    let pred_s0_b32 =
        //        u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //    let pred_s1_b32 =
        //        u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //    //let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());

        //    let sums_32 = sum_byte_count(pred_s0_b32, pred_s1_b32);

        //    sums = add_tuple3(sums, sums_32);
        //}

        ////println!("in table afr for: {} sec",  start_time.elapsed().as_micros());

        //let c2 = sums.0;
        //let c1 = sums.1;
        //let cm = sums.2;

        let (_c0, c1, c2, cm) = self.count_table();

        let call = self.n();

        //let c0 = call - (c2 + c1 + cm);

        let cnomissing = call - cm;

        //println!("in table afr last: {} sec",  start_time.elapsed().as_micros());

        let maf = ((c2 * 2 + c1) as f64) / ((cnomissing * 2) as f64);
        // TODO: Sometimes fails in rap
        //assert!((0.0 <= maf) && (maf <= 1.0));
        maf
    }

    /// slow
    /// counts are
    /// (d2, n2, d1, n1, d0, n0, dm, nm)
    fn stat_contingency_table_nosimd(
        &self,
        phe: &Phe,
    ) -> (usize, usize, usize, usize, usize, usize, usize, usize) {
        assert_eq!(self.n(), phe.n());

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

    /// tuple is faster than vec; 80%
    // TODO: make len 64; also make padding 64
    //#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    //#[target_feature(enable = "avx2")]
    fn stat_contingency_table(
        &self,
        phe: &Phe,
    ) -> (usize, usize, usize, usize, usize, usize, usize, usize) {
        //use std::time::Instant;
        //let start_time = Instant::now();

        assert_eq!(self.n(), phe.n());

        //println!("in table afr init: {} sec",  start_time.elapsed().as_micros());

        //#[cfg(target_arch = "x86")]
        //use std::arch::x86::*;
        //#[cfg(target_arch = "x86_64")]
        //use std::arch::x86_64::*;
        //use std::convert::TryInto;

        //use core::arch::x86_64::*;

        let n = phe.n();
        let ys = phe.inner();
        let pred_s0m = self.predict_s(0);
        let pred_s1m = self.predict_s(1);

        let mut sums: (usize, usize, usize, usize, usize, usize) = (0usize, 0, 0, 0, 0, 0);
        //let mut sums=vec![0usize; 6];

        // TODO: counting both and substract case could be faster
        fn sum_byte(y: u32, p0: u32, p1: u32) -> (usize, usize, usize, usize, usize, usize) {
            // do not use d0/n0 since could be mixed with padding
            let d2 = y & p0 & p1;
            let n2 = (!y) & p0 & p1;
            let d1 = y & p0 & (!p1);
            let n1 = (!y) & p0 & (!p1);
            let dm = y & (!p0) & p1;
            let nm = (!y) & (!p0) & p1;
            (
                cmatrix::popcnt(d2),
                cmatrix::popcnt(n2),
                cmatrix::popcnt(d1),
                cmatrix::popcnt(n1),
                cmatrix::popcnt(dm),
                cmatrix::popcnt(nm),
            )
        }

        fn add_tuple6(
            sums: (usize, usize, usize, usize, usize, usize),
            sums_: (usize, usize, usize, usize, usize, usize),
        ) -> (usize, usize, usize, usize, usize, usize) {
            (
                sums.0 + sums_.0,
                sums.1 + sums_.1,
                sums.2 + sums_.2,
                sums.3 + sums_.3,
                sums.4 + sums_.4,
                sums.5 + sums_.5,
            )
        }

        for ni in 0..(n / 32 + 1) {
            let pred_s0_b32 =
                u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
            let pred_s1_b32 =
                u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
            let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());

            let sums_32 = sum_byte(ys_b32, pred_s0_b32, pred_s1_b32);

            sums = add_tuple6(sums, sums_32);
            //sums = sums
            //    .iter()
            //    .zip(sums_32.iter())
            //    .map(|(s_, x_)| *s_ + *x_)
            //    .collect();
        }

        //println!("in table afr for: {} sec",  start_time.elapsed().as_micros());

        /*
        fn sum_byte(y: u8, p0: u8, p1: u8) -> Vec<usize> {
            // do not use d0/n0 since could be mixed with padding
            let d2 = y & p0 & p1;
            let n2 = (!y) & p0 & p1;
            let d1 = y & p0 & (!p1);
            let n1 = (!y) & p0 & (!p1);
            let dm = y & (!p0) & p1;
            let nm = (!y) & (!p0) & p1;
            fn pop(x: u8) -> usize {
                unsafe { core::arch::x86_64::_popcnt64(x as i64) as usize }
            }
            vec![pop(d2), pop(n2), pop(d1), pop(n1), pop(dm), pop(nm)]
        }

        // TODO: use tuple should be faster
        let sums = ys
            .iter()
            .zip(pred_s0m.iter())
            .zip(pred_s1m.iter())
            .map(|((y, p0), p1)| sum_byte(*y,*p0,*p1))
            .fold(vec![0usize; 6], |sum, x| {
                sum.iter().zip(x.iter()).map(|(s_, x_)| *s_ + *x_).collect()
            });
             */

        let d2 = sums.0;
        let n2 = sums.1;
        let d1 = sums.2;
        let n1 = sums.3;
        let dm = sums.4;
        let nm = sums.5;

        let dall = phe.count();
        //let dall = phe.count_simd();
        //let dall = phe.count();
        let nall = n - dall;
        let d0 = dall - (d2 + d1 + dm);
        let n0 = nall - (n2 + n1 + nm);
        //println!("in table afr last: {} sec",  start_time.elapsed().as_micros());

        (d2, n2, d1, n1, d0, n0, dm, nm)
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

    /// self must be all false
    ///
    /// for plink bed
    /// code(count); b0, b1
    /// 00(2); 1,1
    /// 10(1); 1,0
    /// 11(0); 0,0
    /// 01(NA); 0,1
    #[inline]
    fn set_bed_code_init_unchecked(&mut self, bcode: u8, ni: usize) {
        // skip (b0, b1)=(false, false)
        // large proportion should be 11
        if bcode != 3 {
            let b0 = (bcode & 1) == 0;
            let b1 = (bcode & 2) == 0;
            self.genot_inner_mut()
                .set_by_bool_init_unchecked_v(b0, b1, ni);
        }
    }

    #[inline]
    fn set_bed_code_unchecked(&mut self, bcode: u8, ni: usize) {
        let b0 = (bcode & 1) == 0;
        let b1 = (bcode & 2) == 0;
        self.genot_inner_mut().set_by_bool_unchecked_v(b0, b1, ni);
    }

    fn as_genot_snv_mut_snv(&mut self) -> GenotSnvMut {
        GenotSnvMut::new(self.genot_inner_mut().as_cvec_mut_v())
    }

    fn fill_missing_mode(&mut self) {
        // count 0,1,2
        let mut counts_allele = vec![0usize; 4];

        let n = self.n();
        for ni in 0..n {
            counts_allele[self.get_val_unchecked(ni) as usize] += 1;
        }

        let mut mode: usize = 4;
        let mut mode_counts = 0;
        for i in 0..=2 {
            if counts_allele[i] > mode_counts {
                mode_counts = counts_allele[i];
                mode = i;
            }
        }
        let mode = mode as u8;
        assert_ne!(mode, 4);

        for ni in 0..n {
            if self.get_val_unchecked(ni) == 3 {
                self.set_unchecked(mode, ni);
            }
        }
        // check all are non-missing?
        // -> performance...
    }

    fn fill_missing_mode_maf(&mut self, maf: f64) {
        let mode: u8 = super::maf_to_mode(maf);
        //let mode: u8 = if maf < 1.0 / 3.0f64 {
        //    0
        //} else if maf > 2.0 / 3.0f64 {
        //    2
        //} else {
        //    1
        //};

        let n = self.n();

        for ni in 0..n {
            if self.get_val_unchecked(ni) == 3 {
                self.set_unchecked(mode, ni);
            }
        }
        // check all are non-missing?
        // -> performance...
    }

    // {0:2, 1:1, 2:0, 3:3}
    const REV_AR: [u8; 4] = [2, 1, 0, 3];
    fn reverse_allele(&mut self) {
        let n = self.n();
        for ni in 0..n {
            let val = self.get_val_unchecked(ni);
            self.set_unchecked(Self::REV_AR[val as usize], ni);
        }
    }

    /// for group aggregation
    #[inline]
    fn or_binary(&mut self, other: &GenotSnvRef) {
        // *self cannot be 2*
        // self has 0/1/3
        // other can be 0/1/2/3
        // self  \   others
        //      0 1 2 3
        // -----------
        // 0 | 0 1 1 3
        // 1 |  1 1 1 3
        // 3 | 3 3 3 3
        //

        for ni in 0..self.n() {
            let val = self.get_val_unchecked(ni);
            if val == 3 {
                self.set_unchecked(3, ni);
            } else {
                let val_other = other.get_val_unchecked(ni);
                if val_other == 3 {
                    self.set_unchecked(3, ni);
                } else {
                    // val = 0/1
                    // val_other = 0/1/2
                    self.set_unchecked(val | ((val_other != 0) as u8), ni);
                }
            }
        }

        //self.genot_inner_mut()
        //    .or_binary_v(&g.genot_inner().as_cvec_ref_v());
    }
}
