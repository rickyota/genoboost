//! Compressed Matrix Mutable Reference
//!

use super::B8;
use super::{BaseCMatrix, BaseCMatrixMut};

/// Compressed Matrix Mut
///
#[derive(Debug, PartialEq, Eq)]
pub struct CMatrixMut<'a> {
    inner: &'a mut [B8],
    row_n: usize,
    col_n: usize,
    digit_n: usize,
}

impl<'a> BaseCMatrix for CMatrixMut<'a> {
    #[inline]
    fn inner(&self) -> &[B8] {
        &self.inner
    }
    #[inline]
    fn col_n(&self) -> usize {
        self.col_n
    }
    #[inline]
    fn row_n(&self) -> usize {
        self.row_n
    }
    #[inline]
    fn digit_n(&self) -> usize {
        self.digit_n
    }
}

impl<'a> BaseCMatrixMut for CMatrixMut<'a> {
    #[inline]
    fn inner_mut(&mut self) -> &mut [B8] {
        &mut self.inner
    }
}

/*
// TODO: move to iterator.rs
// TODO:  Item=CVec<'a>
impl<'a> Iterator for RowIterMut<'a> {
    type Item = CMatrixMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.row_i < self.cmatrix.row_n() {
            let row_i = self.row_i;
            self.row_i += 1;
            // lifetime error
            //return Some(self.cmatrix.rows_iter_mut(mi, mi + 1));

            let row_begin = row_i;
            let row_end = row_i + 1;
            let index_begin = self.cmatrix.index_row(row_begin);
            let index_end = self.cmatrix.index_row(row_end);
            let col_n = self.cmatrix.col_n();
            let digit = self.cmatrix.digit_n();

            unsafe {
                let s = &mut self.cmatrix.inner_mut()[index_begin..index_end];
                let ptr = s.as_mut_ptr();
                let len = index_end - index_begin;
                let p = std::slice::from_raw_parts_mut(ptr, len);

                return Some(CMatrixMut {
                    inner: p,
                    row_n: row_end - row_begin,
                    col_n,
                    digit_n: digit,
                });
            }
        }
        None
    }
}
*/

/*
// Is this possible?
// mut scalar
// -> no, shouldn't
pub struct ScalIterMut<'a> {
    cmatrix: CMatrixMut<'a>,
    index_iter: ScalIndexIter,
}
*/

impl<'a> CMatrixMut<'a> {
    pub fn new(inner: &'a mut [B8], row_n: usize, col_n: usize, digit_n: usize) -> Self {
        let m = Self {
            inner,
            row_n,
            col_n,
            digit_n,
        };
        m.check_new();
        m
    }
}

#[cfg(test)]
mod tests {
    use super::super::BaseCVec;
    use super::super::CMatrix;
    use super::*;

    #[test]
    fn test_row_iter_mut() {
        let vec = vec![0, 1, 2, 1, 2, 3];
        let mut m = CMatrix::new(2, 3, 2, &vec);

        let mut iter = m.iter_row_mut();
        assert_eq!(iter.next().unwrap().inner_vals_col_v(), vec![0, 1, 2]);
        assert_eq!(iter.next().unwrap().inner_vals_col_v(), vec![1, 2, 3]);
    }

    #[test]
    fn test_rows_mut_twice() {
        let vec = vec![0, 1, 2, 1, 2, 3];
        let m_ans = CMatrix::new(2, 3, 2, &vec);
        let mut m = CMatrix::new(2, 3, 2, &vec);
        let mut mref = m.as_cmatrix_mut();
        let mut mref2 = mref.rows_mut(0, 1);
        mref2.set_unchecked(3, 0, 0);
        assert_ne!(m, m_ans);
    }

    #[test]
    fn test_split_genot() {
        let vec = vec![0, 1, 2, 1, 2, 3, 1, 1, 1];
        let mut m0_ans = CMatrix::new(2, 3, 2, &[0, 1, 2, 1, 2, 3]);
        let mut m1_ans = CMatrix::new(1, 3, 2, &[1, 1, 1]);

        let mut m = CMatrix::new(3, 3, 2, &vec);
        let mut mref = m.as_cmatrix_mut();
        let (mref0, mref1) = mref.split_rows(1);
        assert_ne!(mref0, m0_ans.as_cmatrix_mut());
        assert_ne!(mref1, m1_ans.as_cmatrix_mut());
    }
}
