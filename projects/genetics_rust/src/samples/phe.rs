use super::base_phe::{BasePhe, BasePheMut};
use cmatrix::prelude::*;

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Phe(CBits);

// Phe
impl BasePhe for Phe {
    type Inner = CBits;
    fn phe_inner(&self) -> &Self::Inner {
        &self.0
    }
}

impl BasePheMut for Phe {
    fn phe_inner_mut(&mut self) -> &mut Self::Inner {
        &mut self.0
    }
}

impl Phe {
    pub fn new(vec: Vec<bool>) -> Phe {
        Phe(CBits::new(vec))
    }
    pub fn new_empty(n: usize) -> Phe {
        Phe(CBits::new_empty(n))
    }
}

/* // PheRef
impl<'a> BasePhe for PheRef<'a> {
    type Inner = CBitsRef<'a>;
    fn phe_inner(&self) -> &Self::Inner {
        &self.0
    }
} */

/* // PheMut
impl<'a> BasePhe for PheMut<'a> {
    type Inner = CBitsMut<'a>;
    fn phe_inner(&self) -> &Self::Inner {
        &self.0
    }
}

impl<'a> BasePheMut for PheMut<'a> {
    fn phe_inner_mut(&mut self) -> &mut Self::Inner {
        &mut self.0
    }
}
 */
