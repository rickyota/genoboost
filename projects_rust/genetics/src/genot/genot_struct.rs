//! module for predictions-related func.
//! len: length of vector
//! size: sample size or snv size
//! predictions shape is m * 2 * n and different from 2 * m * n of C++ ver.
//! This does not require m when getting index of part of predictions
// TODO: many of functions here should move to genetics_rust::{plink,genotype}
// rename to GenotTwin->Genot

use super::base_genot::{BaseGenot, BaseGenotMut, BaseGenotSnv, BaseGenotSnvMut};
use cmatrix::prelude::*;

// 8 x 1 bit
pub type B8 = u8;
// 4 x 2 bit
pub type B8_2 = u8;
// si=0 : dominant
// si=1 : recessive
pub const THRESHOLD_SNV: [f64; 2] = [0.5, 1.5];

/*
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

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Genot(CMatrix);

#[derive(Debug, PartialEq, Eq)]
pub struct GenotMut<'a>(CMatrixMut<'a>);

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct GenotSnv(CVec);
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct GenotSnvRef<'a>(CVecRef<'a>);

#[derive(Debug, PartialEq, Eq)]
pub struct GenotSnvMut<'a>(CVecMut<'a>);

/*
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct GenotSnvModel(CBits);
 */

impl BaseGenot for Genot {
    type Inner = CMatrix;
    fn genot_inner(&self) -> &Self::Inner {
        &self.0
    }
}

impl BaseGenotMut for Genot {
    //type Inner = CMatrix;
    fn genot_inner_mut(&mut self) -> &mut Self::Inner {
        &mut self.0
    }
}

impl Genot {
    #[inline]
    pub fn len_arg(m: usize, n: usize) -> usize {
        CMatrix::len_arg(m, n, 2)
    }

    pub fn byte(m: usize, n: usize) -> usize {
        // len = 1 byte
        Genot::len_arg(m, n)
    }

    pub fn new(m: usize, n: usize, vec: &[u8]) -> Genot {
        Genot(CMatrix::new(m, n, 2, vec))
    }
    pub fn new_zeros(m: usize, n: usize) -> Genot {
        Genot(CMatrix::new_zeros(m, n, 2))
    }

    pub fn extract_snvs(self, use_snvs: &[bool]) -> Genot {
        let Genot(m) = self;
        let m_use = m.extract_snvs(use_snvs);
        Genot(m_use)
    }
}

impl<'a> BaseGenot for GenotMut<'a> {
    type Inner = CMatrixMut<'a>;
    fn genot_inner(&self) -> &Self::Inner {
        &self.0
    }
}

impl<'a> BaseGenotMut for GenotMut<'a> {
    fn genot_inner_mut(&mut self) -> &mut Self::Inner {
        &mut self.0
    }
}

impl<'a> GenotMut<'a> {
    pub fn new(cmat_mut: CMatrixMut) -> GenotMut {
        GenotMut(cmat_mut)
    }
}

impl BaseGenot for GenotSnv {
    type Inner = CVec;
    fn genot_inner(&self) -> &Self::Inner {
        &self.0
    }
}

impl BaseGenotMut for GenotSnv {
    fn genot_inner_mut(&mut self) -> &mut Self::Inner {
        &mut self.0
    }
}
impl BaseGenotSnv for GenotSnv {}
impl BaseGenotSnvMut for GenotSnv {}

impl GenotSnv {
    pub fn new(vec: &[u8]) -> GenotSnv {
        GenotSnv(CVec::new(2, vec))
    }
    pub fn new_empty(n: usize) -> GenotSnv {
        GenotSnv(CVec::new_empty(n, 2))
    }
}

impl<'a> BaseGenot for GenotSnvRef<'a> {
    type Inner = CVecRef<'a>;
    fn genot_inner(&self) -> &Self::Inner {
        &self.0
    }
}
impl<'a> BaseGenotSnv for GenotSnvRef<'a> {}

impl<'a> GenotSnvRef<'a> {
    pub fn new(cvec_ref: CVecRef) -> GenotSnvRef {
        GenotSnvRef(cvec_ref)
    }
    pub fn from_cvecref(cvec: CVecRef) -> GenotSnvRef {
        GenotSnvRef(cvec)
    }
}

impl<'a> BaseGenot for GenotSnvMut<'a> {
    type Inner = CVecMut<'a>;
    fn genot_inner(&self) -> &Self::Inner {
        &self.0
    }
}

impl<'a> BaseGenotMut for GenotSnvMut<'a> {
    fn genot_inner_mut(&mut self) -> &mut Self::Inner {
        &mut self.0
    }
}

impl<'a> BaseGenotSnv for GenotSnvMut<'a> {}
impl<'a> BaseGenotSnvMut for GenotSnvMut<'a> {}

impl<'a> GenotSnvMut<'a> {
    pub fn new(cvec: CVecMut) -> GenotSnvMut {
        GenotSnvMut(cvec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::samples::phe::Phe;

    #[test]
    fn test_new() {
        let vec = vec![1, 2, 3, 0, 1, 2];
        Genot::new(2, 3, &vec);
    }

    #[test]
    fn test_stat_contingency_table() {
        let vec = vec![2, 2, 1, 3, 3, 2, 0, 0, 3];
        let g = GenotSnv::new(&vec);

        let y = vec![true, true, true, true, true, false, false, false, false];
        let phe = Phe::new(&y);

        let t8 = g.stat_contingency_table_nosimd(&phe);

        assert_eq!(t8, (2usize, 1, 1, 0, 0, 2, 2, 1));
    }

    #[test]
    fn test_stat_contingency_table_simd() {
        let vec = vec![2, 2, 1, 3, 3, 2, 0, 0, 3];
        let g = GenotSnv::new(&vec);

        let y = vec![true, true, true, true, true, false, false, false, false];
        let phe = Phe::new(&y);

        let t8 = g.stat_contingency_table(&phe);
        let t8_nosimd = g.stat_contingency_table_nosimd(&phe);

        assert_eq!(t8, t8_nosimd);
    }

    fn is_eq_f64(v: f64, w: f64, e: f64) -> bool {
        (v - w).abs() < e
    }

    #[test]
    fn test_maf_simd() {
        let vec = vec![2, 2, 1, 3, 3, 2, 0, 0, 3];
        let g = GenotSnv::new(&vec);

        let maf = g.maf();

        // (2*3+1)/(6*2)
        assert!(is_eq_f64(maf, 7.0f64 / 12.0f64, 1e-7))

        //wrong
        // (2*3+1)/6
        //assert!(is_eq_f64(maf, 7.0f64 / 6.0f64, 1e-7))
    }

    #[test]
    fn test_set_bed_code_init_unchecked() {
        let mut g = Genot::new_zeros(2, 3);
        let mut gref = g.as_genot_snv_mut(0);

        gref.set_bed_code_init_unchecked(3, 1);
        assert_eq!(gref.get_val(1), 0);

        gref.set_bed_code_init_unchecked(2, 1);
        assert_eq!(gref.get_val(1), 1);

        gref.set_bed_code_init_unchecked(0, 1);
        assert_eq!(gref.get_val(1), 2);

        // cannot overwrite by _init
        gref.set_bed_code_init_unchecked(1, 1);
        assert_ne!(gref.get_val(1), 3);
    }

    #[test]
    fn test_set_bed_code_unchecked() {
        let mut g = Genot::new_zeros(2, 3);
        let mut gref = g.as_genot_snv_mut(0);

        gref.set_bed_code_unchecked(3, 1);
        assert_eq!(gref.get_val(1), 0);

        gref.set_bed_code_unchecked(2, 1);
        assert_eq!(gref.get_val(1), 1);

        gref.set_bed_code_unchecked(0, 1);
        assert_eq!(gref.get_val(1), 2);

        gref.set_bed_code_unchecked(1, 1);
        assert_eq!(gref.get_val(1), 3);
    }

    #[test]
    fn test_fill_missing_mode() {
        let vec = vec![0, 1, 1, 2, 1, 3];
        let vec_exp = vec![0u8, 1, 1, 2, 1, 1];
        let mut g = GenotSnv::new(&vec);

        g.as_genot_snv_mut_snv().fill_missing_mode();

        assert_eq!(&g.vals(), &vec_exp);
    }

    #[test]
    fn test_fill_missing_mode2() {
        // when missing is the mode
        let vec = vec![0, 1, 1, 3, 3, 3];
        let vec_exp = vec![0, 1, 1, 1, 1, 1];
        let mut g = GenotSnv::new(&vec);

        g.as_genot_snv_mut_snv().fill_missing_mode();

        assert_eq!(&g.vals(), &vec_exp);
    }

    //TODO: fill_missing_mode_maf()

    #[test]
    fn test_reverse_allele() {
        let vec = vec![0, 1, 1, 2, 1, 3];
        let vec_exp = vec![2u8, 1, 1, 0, 1, 3];
        let mut g = GenotSnv::new(&vec);

        g.as_genot_snv_mut_snv().reverse_allele();

        assert_eq!(&g.vals(), &vec_exp);
    }
}
