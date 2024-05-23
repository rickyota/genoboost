//use super::phe::PheIter;
use cmatrix::prelude::*;

pub trait BasePhe {
    type Inner: BaseCBits;
    fn phe_inner(&self) -> &Self::Inner;

    #[inline]
    fn n(&self) -> usize {
        self.phe_inner().col_n()
    }
    // byte len
    #[inline]
    fn len_n(&self) -> usize {
        self.phe_inner().col_len()
    }
    #[inline]
    fn inner(&self) -> &[B8] {
        self.phe_inner().inner()
    }
    #[inline]
    fn get(&self, ni: usize) -> bool {
        self.phe_inner().get_b(ni)
    }
    #[inline]
    fn get_unchecked(&self, ni: usize) -> bool {
        self.phe_inner().get_unchecked_b(ni)
    }
    fn count(&self) -> usize {
        self.phe_inner().count()
    }
    /*     fn count_simd(&self) -> usize {
        self.phe_inner().count_simd()
    } */
    fn count_false(&self) -> usize {
        self.phe_inner().count_false()
    }
    fn count_label(&self, is_case: bool) -> usize {
        if is_case {
            self.count()
        } else {
            self.count_false()
        }
    }
    fn iter(&self) -> PheIter {
        PheIter::new(self.phe_inner().iter())
    }

    /// len=n without padding
    fn inner_f64(&self) -> Vec<f64> {
        self.iter().map(|x| x as i32 as f64).collect()
    }
    fn inner_i32(&self) -> Vec<i32> {
        self.iter().map(|x| x as i32).collect()
    }
    fn inner_bool(&self) -> Vec<bool> {
        self.iter().collect()
    }
}

pub trait BasePheMut: BasePhe
where
    Self::Inner: BaseCBitsMut,
{
    fn phe_inner_mut(&mut self) -> &mut Self::Inner;

    #[inline]
    fn set_unchecked(&mut self, b: bool, ni: usize) {
        self.phe_inner_mut().set_bool_unchecked_b(b, ni);
    }

    #[inline]
    fn set(&mut self, b: bool, ni: usize) {
        self.phe_inner_mut().set_bool_b(b, ni);
    }
}

pub struct PheIter<'a>(BoolIter<'a>);

impl<'a> PheIter<'a> {
    pub fn new(iter: BoolIter) -> PheIter {
        PheIter(iter)
    }
}

impl<'a> Iterator for PheIter<'a> {
    type Item = bool;
    fn next(&mut self) -> Option<Self::Item> {
        return self.0.next();
    }
}
