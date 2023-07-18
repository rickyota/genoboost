//! Compressed Matrix
//!
//! Matrix, Vector, Scalar, Bool
//!

use super::B8;
use super::{BaseCMatrix, BaseCMatrixMut};
//use crate::alloc;

/// Compressed Matrix
/// Mainly used for genotype and phenotype.
///
/// ax2 is the same as matrix.
/// Row values are compressed to ax0 and ax1 axes.
/// ax0 is a vec of bits, where one bit corresponds to one row.
/// By collecting ax1 bits, a value is expressed.
///
/// ex.  row=2, col=3, val=[[0,1,2],[1,2,3]] (in 2bit)
/// ax0=2, ax1=2(fixed), ax2=3
/// ax1=0 shows dominant model, and ax1=1 shows recessive model
/// inner = [ 0b11000000, 0, 0, 0, 0, 0b01100000, 0, 0, 0, 0, 0b01100000, 0, 0, 0, 0, 0b00100000, 0, 0, 0, 0 ]
///
///
///
/// index is same as np.array https://docs.rs/ndarray/latest/ndarray/struct.ArrayBase.html#indexing-and-dimension
/// matrix=[4,3]
/// [[ [0, 0], [0, 1], [0, 2] ],  // row 0
///   [ [1, 0], [1, 1], [1, 2] ],  // row 1
///   [ [2, 0], [2, 1], [2, 2] ],  // row 2
///   [ [3, 0], [3, 1], [3, 2] ]]  // row 3
///           \       \       \
///   column 0  \     column 2
///            column 1
///
/// ref: https://athemathmo.github.io/rulinalg/doc/rulinalg/matrix/struct.Matrix.html
///
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct CMatrix {
    inner: Vec<B8>,
    // row: m
    row_n: usize,
    // column: n
    col_n: usize,
    // digit of bits; s
    digit_n: usize,
}

/*
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct CMatrixRef<'a> {
    inner: &'a [B8],
    row_n: usize,
    col_n: usize,
    digit: usize,
}

#[derive(Debug, PartialEq, Eq)]
pub struct CMatrixMut<'a> {
    inner: &'a mut [B8],
    row_n: usize,
    col_n: usize,
    digit: usize,
}
 */

impl BaseCMatrix for CMatrix {
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

impl<'a> BaseCMatrixMut for CMatrix {
    fn inner_mut(&mut self) -> &mut [B8] {
        &mut self.inner
    }
}

/*
// impossible
impl<'a> AsRef<CMatrixRef<'a>> for CMatrix {
    #[inline]
    fn as_ref(&self) -> &CMatrixRef<'a> {
        &self.to_cmatrix_ref()
    }
}
*/

impl CMatrix {
    /// padding=32
    /// TODO: -> 64
    #[inline]
    fn col_len_arg(col_n: usize) -> usize {
        col_n / 8 + 5
    }

    #[inline]
    pub fn len_arg(row_n: usize, col_n: usize, digit_n: usize) -> usize {
        let col_len = CMatrix::col_len_arg(col_n);
        col_len * row_n * digit_n
    }

    /// called by outside
    #[inline]
    fn row_index_arg(row: usize, col_n: usize, digit_n: usize) -> (usize, usize) {
        let col_len = CMatrix::col_len_arg(col_n);
        (row * digit_n * col_len, (row + 1) * digit_n * col_len)
    }

    /// Vec is not aligned
    pub fn new_zeros(row_n: usize, col_n: usize, digit_n: usize) -> CMatrix {
        CMatrix::check_new(row_n, col_n, digit_n);

        let len = CMatrix::len_arg(row_n, col_n, digit_n);

        //log::debug!("create vec {}", len);
        let vec: Vec<B8> = vec![0u8; len];
        //let mut vec: Vec<B8> = Vec::with_capacity(len);
        //let mut vec: Vec<B8> = alloc::with_capacity_align_u8(len);
        //vec.resize(len, 0x00);

        let m = CMatrix {
            inner: vec,
            row_n,
            col_n,
            digit_n,
        };

        m.check_new();
        m

        //CMatrix::new_compressed_vec(row_n, col_n, digit, vec)
    }

    // TODO: make vec &[u8]
    /// Assume max(vec) <=3 but not checked.
    /// digit=2 only. use new_vec_bool() for digit=1
    /// input is vec with row*col
    /// val -> bit in func depending on digit
    //pub fn new(row: usize, col: usize, digit: usize, vec: Vec<u8>) -> CMatrix {
    pub fn new(row: usize, col: usize, digit: usize, vec: &[u8]) -> CMatrix {
        CMatrix::check_new(row, col, digit);
        CMatrix::check_new_vec(row, col, digit, vec);
        let mut m = CMatrix::new_zeros(row, col, digit);
        m.set_vec_unchecked(vec);
        m.check_new();
        m
    }

    // for CVec
    // all of digit>=1 is 0x00 (undefined)
    // -> panic
    pub fn new_vec_bool(row: usize, col: usize, vec: &[bool]) -> CMatrix {
        let digit = 1;
        CMatrix::check_new(row, col, digit);
        CMatrix::check_new_vec_bool(row, col, digit, &vec);
        let mut m = CMatrix::new_zeros(row, col, digit);
        m.set_vec_bool(&vec);
        m.check_new();
        m
    }

    /*
    fn val_u8_to_bit(&mut self, vec: &[u8]) {
        let col_n = self.col_n();
        for row_i in 0..self.row_n() {
            let mut row = self.rows_mut(row_i, row_i + 1);
            row.set_col_unchecked(&vec[row_i * col_n..(row_i + 1) * col_n]);
        }
    }

    fn val_bool_to_bit(&mut self, vec: &[bool]) {
        let col_n = self.col_n();
        for row_i in 0..self.row_n() {
            let mut row = self.rows_mut(row_i, row_i + 1);
            row.set_col_bool(&vec[row_i * col_n..(row_i + 1) * col_n]);
        }
    }
     */

    fn check(&self) {
        if self.inner().len() != self.len() {
            panic!("inner length of CMatrix is wrong.")
        }
    }

    fn check_new(row_n: usize, col_n: usize, digit: usize) {
        if row_n == 0 {
            panic!("row >=1 is required.");
        }
        if col_n == 0 {
            panic!("col >=1 is required.");
        }
        if digit == 0 {
            panic!("digit >=1 is required.");
        }
        if digit >= 3 {
            unimplemented!("digit <= 2 is required.");
        }
    }

    fn check_new_vec(row_n: usize, col_n: usize, digit: usize, vec: &[u8]) {
        if digit == 1 {
            panic!("Use CMatrix::new_vec_bool() instead.");
        } else if digit > 2 {
            panic!("Unimplemented for digit>2.");
        }
        // check len
        if row_n * col_n != vec.len() {
            panic!("Length of vec is wrong.");
        }
        // max is 3
    }
    fn check_new_vec_bool(row_n: usize, col_n: usize, digit: usize, vec: &[bool]) {
        if digit != 1 {
            panic!("Wrong digit");
        }

        // check len
        if row_n * col_n != vec.len() {
            panic!("Length of vec is wrong.");
        }
        // digit =1
    }

    pub fn extract_snvs(self, use_snvs: &[bool]) -> CMatrix {
        let CMatrix {
            mut inner,
            row_n,
            col_n,
            digit_n,
        } = self;

        assert_eq!(row_n, use_snvs.len());

        // overwrite inner
        // ok since overwriting from top
        // use retain?
        // ri: row_index of converted row
        let mut ri_to = 0;
        // ri: row_index of  row before
        for (ri, b) in use_snvs.iter().enumerate() {
            if *b {
                let (index_b, index_e) = CMatrix::row_index_arg(ri, col_n, digit_n);
                let (index_to_b, index_to_e) = CMatrix::row_index_arg(ri_to, col_n, digit_n);
                // replace like inner[index_to_b..index_to_e] = inner[index_b..index_e];
                // TODO: how to avoid clone?
                let a = inner[index_b..index_e].to_owned();
                //let a = &inner[index_b..index_e];
                inner.splice(index_to_b..index_to_e, a);
                ri_to += 1;
            }
        }

        let row_n_use = use_snvs.iter().filter(|b| **b).count();
        assert_eq!(ri_to, row_n_use);

        let len_new = CMatrix::len_arg(row_n_use, col_n, digit_n);
        // make size and capacity fit
        inner.truncate(len_new);
        inner.shrink_to_fit();

        let m = CMatrix {
            inner,
            row_n: row_n_use,
            col_n,
            digit_n,
        };
        m.check_new();
        m.check();
        m
    }
}

/*
// These index should not be implemented
// since index should return **actual** reference to the element.
//
// https://stackoverflow.com/questions/33770989/implementing-the-index-operator-for-matrices-with-multiple-parameters
// using tuple is an option
//impl Index<(usize, usize)> for CMatrix {
// IndexMut cannot be implemented since `let x=v[1]` should cannot be done here...
impl<'a> Index<[usize; 2]> for CMatrixMut<'a> {
    type Output = u8;
    #[inline]
    fn index(&self, idx: [usize; 2]) -> &Self::Output {
        let [row_i, ni] = idx;
        let digit0 = self.inner_val_bool(row_i, ni, 0);
        let digit1 = self.inner_val_bool(row_i, ni, 1);

        if (!digit0) & (!digit1) {
            &0
        } else if digit0 & (!digit1) {
            &1
        } else if digit0 & digit1 {
            &2
        } else {
            &3
        }
    }
}
impl<'a> Index<[usize; 3]> for CMatrixMut<'a> {
    type Output = bool;
    #[inline]
    fn index(&self, idx: [usize; 3]) -> &Self::Output {
        let [row_i, ni, digit] = idx;
        let b = self.inner_val_bool(row_i, ni, digit);
        //&b.clone()
        //&b
        match b {
            true => &true,
            false => &false,
        }
    }
}

// this might be allowed?
impl<'a> Index<std::ops::Range<usize>> for CMatrixRef<'a> {
    type Output = CMatrixRef<'a>;
    fn index(&self, idx: std::ops::Range<usize>) -> &Self::Output {
        let first = idx.next().unwrap();
        let last = idx.last().unwrap();
        &self.rows(first, last + 1)
    }
}

*/

#[cfg(test)]
mod tests {
    use super::*;

    fn setup0() -> CMatrix {
        CMatrix::new_zeros(2, 3, 2)
    }

    fn setup1() -> CMatrix {
        let vec = vec![1, 2, 3, 0, 1, 2];
        CMatrix::new(2, 3, 2, &vec)
    }

    fn setup2() -> CMatrix {
        let vec = vec![1, 2, 3, 0, 1, 2];
        CMatrix::new(3, 2, 2, &vec)
    }

    #[test]
    fn test_new_len() {
        let row = 5;
        let col = 15;
        let digit = 2;
        let len_ans = 60;
        assert_eq!(len_ans, CMatrix::new_zeros(row, col, digit).len());
    }

    /*     #[test]
    fn test_new_large() {
        let row = 10000;
        let col = 100000;
        let digit = 2;
        //let len_ans = 60;
        //assert_eq!(len_ans, CMatrix::new_empty(row, col, digit).len());
        CMatrix::new_empty(row, col, digit);
    } */

    #[test]
    fn test_new_vec() {
        let m = setup1();
        assert_eq!(m.get_val(0, 0), 1);
        assert_eq!(m.get_val(0, 1), 2);
        assert_eq!(m.get_val(0, 2), 3);
        assert_eq!(m.get_val(1, 0), 0);
        assert_eq!(m.get_val(1, 1), 1);
        assert_eq!(m.get_val(1, 2), 2);
    }

    #[test]
    #[should_panic]
    fn test_new_vec_wrong_len() {
        let vec = vec![1, 2, 3, 0, 1];
        CMatrix::new(2, 3, 2, &vec);
    }

    /*     #[test]
    #[should_panic]
    fn test_new_vec_large() {
        let vec = vec![4, 2, 3, 0, 1, 2];
        CMatrix::new(2, 3, 2, vec);
    } */

    #[test]
    fn test_new_bool() {
        let vec = vec![true, false, false, true];
        let m = CMatrix::new_vec_bool(2, 2, &vec);
        assert_eq!(m.get(0, 0, 0), true);
        assert_eq!(m.get(0, 1, 0), false);
        assert_eq!(m.get(1, 0, 0), false);
        assert_eq!(m.get(1, 1, 0), true);
    }

    #[test]
    fn test_as_cmatrix_mut() {
        let mut m = setup1();
        let mref = m.as_cmatrix_mut();
        assert_eq!(mref.row_n(), 2);
    }

    #[test]
    fn test_inner() {
        // row=2, col=3, vec = [1, 2, 3, 0, 1, 2];
        // row0, dom = 0b11000000
        // row0, rec = 0b01100000
        // row1, dom = 0b01100000
        // row1, rec = 0b00100000

        let m = setup1();
        assert_eq!(
            m.inner(),
            [
                0b11000000, 0, 0, 0, 0, 0b01100000, 0, 0, 0, 0, 0b01100000, 0, 0, 0, 0, 0b00100000,
                0, 0, 0, 0
            ]
        );
    }

    #[test]
    fn test_get() {
        let m = setup1();
        //let vec = vec![1, 2, 3, 0, 1, 2];
        assert_eq!(m.get(0, 0, 0), true);
        assert_eq!(m.get(0, 0, 1), false);
        assert_eq!(m.get(0, 1, 0), true);
        assert_eq!(m.get(0, 1, 1), true);
        assert_eq!(m.get(0, 2, 0), false);
        assert_eq!(m.get(0, 2, 1), true);
        assert_eq!(m.get(1, 0, 0), false);
        assert_eq!(m.get(1, 0, 1), false);
    }

    #[test]
    #[should_panic]
    fn test_get_val_outside_col() {
        let m = setup1();
        m.get_val(0, 4);
    }

    #[test]
    #[should_panic]
    fn test_get_val_outside_row() {
        let m = setup1();
        m.get_val(3, 0);
    }

    #[test]
    #[should_panic]
    fn test_get_outside_digit() {
        let m = setup1();
        m.get(0, 0, 2);
    }

    #[test]
    fn test_inner_col_digit() {
        let m = setup1();

        //let vec = [[1, 2, 3], [0, 1, 2]];
        assert_eq!(m.inner_col_digit(0, 0), &[0b110_00000u8, 0, 0, 0, 0]);
    }

    #[test]
    fn test_inner_vals_col_digit() {
        let m = setup1();

        //let vec = vec![1, 2, 3, 0, 1, 2];
        assert_eq!(m.inner_vals_col_digit(0, 0), vec![true, true, false]);
    }

    #[test]
    fn test_inner_vals_col_digit_u8() {
        let m = setup1();

        //let vec = vec![1, 2, 3, 0, 1, 2];
        assert_eq!(m.inner_vals_col_digit_u8(0, 0), vec![1, 1, 0]);
    }

    #[test]
    fn test_inner_vals_col() {
        let m = setup1();

        //let vec = vec![1, 2, 3, 0, 1, 2];
        assert_eq!(m.inner_vals_col(0), vec![1, 2, 3]);
        assert_eq!(m.inner_vals_col(1), vec![0, 1, 2]);
    }

    #[test]
    fn test_rows() {
        let m = setup2();

        //vec = vec![[1, 2], [3, 0], [1, 2]];
        let rs = m.rows(1, 3);
        assert_eq!(rs.inner_vals_col(0), vec![3, 0]);
        assert_eq!(rs.inner_vals_col(1), vec![1, 2]);
    }

    /// empty rows
    #[test]
    fn test_rows_empty() {
        let m = setup2();

        let rs = m.rows(1, 1);
        assert_eq!(rs.row_n(), 0);
        assert_eq!(rs.inner().len(), 0);
    }

    #[test]
    fn test_as_vec_ref() {
        let m = setup2();

        //vec = vec![[1, 2], [3, 0], [1, 2]];
        let rs = m.as_cvec_ref(1);
        assert_eq!(rs.inner_vals_col(0), vec![3, 0]);
    }

    #[test]
    fn test_set() {
        let mut m = setup1();
        m.set(3, 1, 2);
        assert_eq!(m.get_val(1, 2), 3);
    }

    #[test]
    fn test_set_by_bool_unchecked() {
        let mut m = setup1();

        m.set_by_bool_unchecked(false, false, 1, 2);
        assert_eq!(m.get_val(1, 2), 0);

        m.set_by_bool_unchecked(true, false, 1, 2);
        assert_eq!(m.get_val(1, 2), 1);

        m.set_by_bool_unchecked(true, true, 1, 2);
        assert_eq!(m.get_val(1, 2), 2);

        m.set_by_bool_unchecked(false, true, 1, 2);
        assert_eq!(m.get_val(1, 2), 3);
    }

    #[test]
    fn test_set_by_bool_init_unchecked() {
        let mut m = setup0();

        m.set_by_bool_init_unchecked(false, false, 1, 2);
        assert_eq!(m.get_val(1, 2), 0);

        m.set_by_bool_init_unchecked(true, false, 1, 2);
        assert_eq!(m.get_val(1, 2), 1);

        m.set_by_bool_init_unchecked(true, true, 1, 2);
        assert_eq!(m.get_val(1, 2), 2);

        // cannot overwrite by _init
        m.set_by_bool_init_unchecked(false, true, 1, 2);
        assert_ne!(m.get_val(1, 2), 3);

        // but this will do
        m.set_by_bool_unchecked(false, true, 1, 2);
        assert_eq!(m.get_val(1, 2), 3);
    }
}
