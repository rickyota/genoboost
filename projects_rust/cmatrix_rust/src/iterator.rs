//! Iterator over index in CMatrix and over rows, values
//!
//!

use super::{BaseCBits, BaseCMatrix, BaseCMatrixMut, BaseCVec};
use super::{CBitsRef, CMatrixMut, CMatrixRef, CVecMut, CVecRef};

pub struct BoolIter<'a> {
    cbits: CBitsRef<'a>,
    col_i: usize,
}

impl<'a> Iterator for BoolIter<'a> {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        if self.col_i < self.cbits.col_n() {
            let col_i = self.col_i;
            self.col_i += 1;
            return Some(self.cbits.get_b(col_i));
        }
        None
    }
}

impl<'a> BoolIter<'a> {
    pub fn new(cvec: CBitsRef<'a>) -> BoolIter {
        BoolIter {
            cbits: cvec,
            col_i: 0,
        }
    }
}

pub struct ScalIter<'a> {
    cbits: CVecRef<'a>,
    col_i: usize,
}

impl<'a> ScalIter<'a> {
    pub fn new(cvec: CVecRef<'a>) -> ScalIter {
        ScalIter {
            cbits: cvec,
            col_i: 0,
        }
    }
}

// return vec
// Iterator is also unsafe in Vec
// https://doc.rust-lang.org/src/core/slice/iter.rs.html#86
impl<'a> Iterator for ScalIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.col_i < self.cbits.col_n() {
            let col_i = self.col_i;
            self.col_i += 1;
            return Some(self.cbits.get_val_unchecked_v(col_i));
        }
        None
    }
}

pub struct RowIter<'a> {
    cmatrix: CMatrixRef<'a>,
    row_i: usize,
}

impl<'a> RowIter<'a> {
    pub fn new(cmatrix: CMatrixRef<'a>) -> RowIter {
        RowIter { cmatrix, row_i: 0 }
    }
}

// return vec
// Iterator is also unsafe in Vec
// https://doc.rust-lang.org/src/core/slice/iter.rs.html#86
impl<'a> Iterator for RowIter<'a> {
    type Item = CVecRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.row_i < self.cmatrix.row_n() {
            let row_i = self.row_i;
            self.row_i += 1;

            let s = self.cmatrix.inner_row(row_i);
            //let s = self.cmatrix.inner_rows(row_i, row_i + 1);
            let len = s.len();
            let ptr = s.as_ptr();

            unsafe {
                let p = std::slice::from_raw_parts(ptr, len);
                return Some(CVecRef::new(
                    p,
                    self.cmatrix.col_n(),
                    self.cmatrix.digit_n(),
                ));
            }

            //return Some(self.cmatrix.to_cvec_ref(row_i));
        }
        None
    }
}

pub struct RowIterMut<'a> {
    cmatrix: CMatrixMut<'a>,
    row_i: usize,
}

impl<'a> RowIterMut<'a> {
    pub fn new(cmatrix: CMatrixMut<'a>) -> RowIterMut {
        RowIterMut { cmatrix, row_i: 0 }
    }
}

// return vec
impl<'a> Iterator for RowIterMut<'a> {
    type Item = CVecMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.row_i < self.cmatrix.row_n() {
            let row_i = self.row_i;
            self.row_i += 1;
            // lifetime error
            //return Some(self.cmatrix.rows_iter_mut(mi, mi + 1));

            let (col_n, digit_n) = (self.cmatrix.col_n(), self.cmatrix.digit_n());

            let s = self.cmatrix.inner_row_mut(row_i);
            //let s = self.cmatrix.inner_rows_mut(row_i, row_i + 1);
            let ptr = s.as_mut_ptr();
            let len = s.len();
            unsafe {
                let p = std::slice::from_raw_parts_mut(ptr, len);

                return Some(CVecMut::new(p, col_n, digit_n));
            }
        }
        None
    }
}

/*
// no opportunity to use
// reimplement when to use
pub struct ScalIndexIter {
    row_i: usize,
    col_i: usize,
    row_n: usize,
    col_n: usize,
}

impl ScalIndexIter {
    pub fn new(row_n: usize, col_n: usize) -> ScalIndexIter {
        ScalIndexIter {
            row_i: 0,
            col_i: 0,
            row_n,
            col_n,
        }
    }
    pub fn new_any(row_i: usize, col_i: usize, row_n: usize, col_n: usize) -> ScalIndexIter {
        ScalIndexIter {
            row_i,
            col_i,
            row_n,
            col_n,
        }
    }
}

/// iterating index(row_i, col_i)
/// (0,0)->(0,1)->...(0,col-1)->(1,0)->...(row-1,col-1)
impl Iterator for ScalIndexIter {
    type Item = [usize; 2];

    fn next(&mut self) -> Option<Self::Item> {
        // checking mi only is fine
        if self.row_i < self.row_n {
            let row_i = self.row_i;
            let col_i = self.col_i;
            if self.col_i + 1 == self.col_n {
                self.row_i += 1;
                self.col_i = 0;
            } else {
                self.col_i += 1;
            }
            return Some([row_i, col_i]);
        }
        None
    }
}

pub struct BoolIndexIter {
    row_i: usize,
    col_i: usize,
    digit_i: usize,
    row_n: usize,
    col_n: usize,
    digit_n: usize,
}

impl BoolIndexIter {
    pub fn new(row_n: usize, col_n: usize, digit_n: usize) -> BoolIndexIter {
        BoolIndexIter {
            row_i: 0,
            col_i: 0,
            digit_i: 0,
            row_n,
            col_n,
            digit_n,
        }
    }
    pub fn new_any(
        row_i: usize,
        col_i: usize,
        digit_i: usize,
        row_n: usize,
        col_n: usize,
        digit_n: usize,
    ) -> BoolIndexIter {
        BoolIndexIter {
            row_i,
            col_i,
            digit_i,
            row_n,
            col_n,
            digit_n,
        }
    }
}

/// iterating index(row_i, digit_i, col_i)
/// (0,0,0)->(0,0,1)->...(0,0,col-1)->(0,1,0)->...(row-1,digit-1,col-1)
impl Iterator for BoolIndexIter {
    type Item = [usize; 3];

    fn next(&mut self) -> Option<Self::Item> {
        // checking mi only is fine
        if self.row_i < self.row_n {
            let row_i = self.row_i;
            let digit_i = self.digit_i;
            let col_i = self.col_i;
            // move next
            if (self.col_i + 1 == self.col_n) & (self.digit_i + 1 == self.digit_n) {
                self.row_i += 1;
                self.digit_i = 0;
                self.col_i = 0;
            } else if self.col_i + 1 == self.col_n {
                self.digit_i += 1;
                self.col_i = 0;
            } else {
                self.col_i += 1;
            }
            return Some([row_i, col_i, digit_i]);
        }
        None
    }
}
  */

/*
pub struct RowIterMut<'a> {
    cmatrix: CMatrixMut<'a>,
    row_i: usize,
}

// FIXME: incomprehensible error
// cannot call func in BaseCMatrix e.g., self.cmatrix.row_n()
impl<'a> RowIterMut<'a> {
    pub fn new(cmatrix: CMatrixMut<'a>) -> RowIterMut {
        RowIterMut { cmatrix, row_i: 0 }
    }

    pub fn new_any(cmatrix: CMatrixMut<'a>, row_i: usize) -> RowIterMut {
        RowIterMut { cmatrix, row_i }
    }
}

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

#[cfg(test)]
mod tests {
    use rayon::iter::ParallelBridge;
    use rayon::iter::ParallelIterator;

    use super::super::{BaseCVec, BaseCVecMut};
    use super::super::{CBits, CMatrix, CVec};
    use super::*;

    #[test]
    fn test_iter() {
        let vec = vec![true, false];
        let m = CBits::new(&vec);

        let mut iter = m.iter();
        assert_eq!(iter.next().unwrap(), true);
        assert_eq!(iter.next().unwrap(), false);
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_iter_scal() {
        let vec = vec![3, 2, 1, 0];
        let m = CVec::new(2, &vec);

        let mut iter = m.iter_scal();
        assert_eq!(iter.next().unwrap(), 3);
        assert_eq!(iter.next().unwrap(), 2);
        assert_eq!(iter.next().unwrap(), 1);
        assert_eq!(iter.next().unwrap(), 0);
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_row_iter() {
        let vec = vec![0, 1, 2, 1, 2, 3];
        let m = CMatrix::new(2, 3, 2, &vec);

        let mut iter = m.iter_row();
        assert_eq!(iter.next().unwrap().inner_vals_col_v(), vec![0, 1, 2]);
        assert_eq!(iter.next().unwrap().inner_vals_col_v(), vec![1, 2, 3]);
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_row_iter_mut() {
        let vec = vec![0, 1, 2, 1, 2, 3];
        let m_ans = CMatrix::new(2, 3, 2, &vec);
        let mut m = CMatrix::new(2, 3, 2, &vec);

        let mut iter = m.iter_row_mut();
        assert_eq!(iter.next().unwrap().inner_vals_col_v(), vec![0, 1, 2]);
        assert_eq!(iter.next().unwrap().inner_vals_col_v(), vec![1, 2, 3]);
        assert_eq!(iter.next(), None);

        let mut iter = m.iter_row_mut();
        iter.next().unwrap().set_v(3, 0);

        assert_ne!(m, m_ans);
    }

    #[test]
    fn test_row_iter_mut_par() {
        let vec = vec![0, 1, 2, 1, 2, 3];
        let mut m = CMatrix::new(2, 3, 2, &vec);

        m.iter_row_mut()
            .par_bridge()
            .for_each(|mut row| row.set_v(3, 0));

        assert_eq!(m.inner_vals_col(0), vec![3, 1, 2]);
        assert_eq!(m.inner_vals_col(1), vec![3, 2, 3]);
    }

    /*
    #[test]
    fn test_scal_index_iter() {
        let mut iter = ScalIndexIter {
            row_i: 0,
            col_i: 0,
            row_n: 2,
            col_n: 3,
        };

        assert_eq!(iter.next().unwrap(), [0, 0]);
        assert_eq!(iter.next().unwrap(), [0, 1]);
        assert_eq!(iter.next().unwrap(), [0, 2]);
        assert_eq!(iter.next().unwrap(), [1, 0]);
        assert_eq!(iter.next().unwrap(), [1, 1]);
        assert_eq!(iter.next().unwrap(), [1, 2]);
        assert_eq!(iter.next(), None);
    }
     */
}
