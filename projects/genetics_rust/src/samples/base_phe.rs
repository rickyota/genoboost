use super::iterator::PheIter;
use cmatrix::prelude::*;

pub trait BasePhe {
    type Inner: BaseCBits;
    fn phe_inner(&self) -> &Self::Inner;

    #[inline]
    fn n(&self) -> usize {
        self.phe_inner().col_n()
    }
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
    fn inner_f64(&self) -> Vec<f64> {
        self.iter().map(|x| x as i32 as f64).collect()
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
