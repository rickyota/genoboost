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
    pub fn new(m: usize, n: usize, vec: Vec<u8>) -> Genot {
        Genot(CMatrix::new(m, n, 2, vec))
    }
    pub fn new_empty(m: usize, n: usize) -> Genot {
        Genot(CMatrix::new_empty(m, n, 2))
    }

    pub fn extract_snvs(self, use_snvs: &[bool]) -> Genot {
        let Genot(m) = self;
        let m_use = m.extract_snvs(use_snvs);
        Genot(m_use)
    }
}

// GenotTwinMut
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

// GenotTwinSnv
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
    pub fn new(vec: Vec<u8>) -> GenotSnv {
        GenotSnv(CVec::new(2, vec))
    }
    pub fn new_empty(n: usize) -> GenotSnv {
        GenotSnv(CVec::new_empty(n, 2))
    }
}

// GenotTwinSnvRef
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

// GenotTwinSnvMut
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

    #[test]
    fn test_new() {
        let vec = vec![1, 2, 3, 0, 1, 2];
        Genot::new(2, 3, vec);
    }
}
