//! Base trait for CMatrix family
//!
//!
// TODO: row_i or row?

// use super::cmatrix_mut::RowIterMut;
//use super::cmatrix_ref::ScalIter;
// use super::iterator::{BoolIndexIter, ScalIndexIter};
use super::iterator::{RowIter, RowIterMut};
use super::{CMatrixMut, CMatrixRef};
use super::{CVecMut, CVecRef, B8};

/// Compressed Matrix
/// Mainly used for genotype and phenotype.
///
/// When accessing index, using iterator is efficient since no index boundary check is run.
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
///
/// Mapping score to digits
/// Anything func is fine but should be consistent to digit0 in digit=1.
///  -> Is this necessary since digit1 never exist.
/// Now,
/// For digit=1 (digit1=0)
/// digit0 <-> scalar
/// 0         <-> 0
/// 1          <-> 1
///
/// For digit=2
/// (digit0, 1) <-> scalar
/// (0       , 0) <-> 0   // consistent to digit=1
/// (1       , 0)  <-> 1    // consistent to digit=1
/// (1       , 1)  <-> 2
/// (0       , 1)  <-> 3
///
///
/// for CMatrix, CMatrixRef, CMatrixMut
/// use unchecked if guaranteed
pub trait BaseCMatrix {
    fn inner(&self) -> &[B8];

    fn col_n(&self) -> usize;

    fn row_n(&self) -> usize;

    fn digit_n(&self) -> usize;

    #[inline]
    fn col_len(&self) -> usize {
        self.col_n() / 8 + 5
    }

    #[inline]
    fn len(&self) -> usize {
        self.col_len() * self.row_n() * self.digit_n()
    }

    /// row_i <= row_n() since = in row_end
    #[inline]
    fn index_row_unchecked(&self, row_i: usize) -> usize {
        row_i * self.digit_n() * self.col_len()
        //row_i * 2 * self.col_len()
    }
    /// needs checker digit_i<self.digit()
    #[inline]
    fn index_row_digit_unchecked(&self, row_i: usize, digit_i: usize) -> usize {
        //row_i * 2 * self.col_len() + digit_i * self.col_len()
        row_i * self.digit_n() * self.col_len() + digit_i * self.col_len()
    }

    #[inline]
    fn index_byte_unchecked(&self, row_i: usize, col_i: usize, digit_i: usize) -> usize {
        self.index_row_digit_unchecked(row_i, digit_i) + (col_i >> 3)
        //row_i * self.digit() * self.col_len() + digit_i * self.col_len() + (col_i >> 3)
        //row_i * 2 * self.col_len() + digit_i * self.col_len() + (col_i >> 3)
    }

    #[inline]
    fn get_byte_unchecked(&self, row_i: usize, col_i: usize, digit_i: usize) -> B8 {
        self.inner()[self.index_byte_unchecked(row_i, col_i, digit_i)]
    }

    /*
    /// index = [row_i, col_i, digit_i]
    #[inline]
    fn get_byte_ar_unchecked(&self, index: [usize; 3]) -> B8 {
        self.get_byte_unchecked(index[0], index[1], index[2])
    }
     */

    #[inline]
    fn get_unchecked(&self, row_i: usize, col_i: usize, digit_i: usize) -> bool {
        (self.get_byte_unchecked(row_i, col_i, digit_i) & (0x80 >> (col_i & 7))) != 0x00
        //(self.inner_col_digit(index[0], index[2])[index[1] >> 3] & (0x80 >> (index[1] & 7))) != 0x00
    }

    #[inline]
    fn get(&self, row_i: usize, col_i: usize, digit_i: usize) -> bool {
        self.check_arg(row_i, col_i, digit_i);

        self.get_unchecked(row_i, col_i, digit_i)
    }

    /*
    /// index = [row_i, col_i, digit_i]
    #[inline]
    fn get_ar_unchecked(&self, index: [usize; 3]) -> bool {
        self.get_unchecked(index[0], index[1], index[2])
    }
     */

    // digit=1 is wrong here
    // panic if digit>2?
    // -> should be but performance...
    // or creating treat for genot and do not contain get_val() there.
    // or in index_byte? -> no, index_byte works in any situation
    #[inline]
    fn get_val_unchecked(&self, row_i: usize, col_i: usize) -> u8 {
        let digit0 = self.get_unchecked(row_i, col_i, 0);
        let digit1 = self.get_unchecked(row_i, col_i, 1);

        if (!digit0) && (!digit1) {
            0
        } else if digit0 && (!digit1) {
            1
        } else if digit0 && digit1 {
            2
        } else {
            3
        }
    }

    #[inline]
    fn get_val(&self, row_i: usize, col_i: usize) -> u8 {
        self.check_arg_digit2();
        self.check_arg(row_i, col_i, 0);

        self.get_val_unchecked(row_i, col_i)
    }

    /*
    /// ar = [row_i, col_i]
    #[inline]
    fn get_val_ar(&self, ar: [usize; 2]) -> u8 {
        self.get_val(ar[0], ar[1])
    }
     */

    fn inner_vals_u8(&self) -> Vec<Vec<u8>> {
        (0..self.row_n())
            .map(|ri| self.inner_vals_col(ri))
            .collect()
        // cannot call fn in CVecRef
        //self.iter_row().map(|r| r.inner_vals_colv_v()).collect()
        //self.iter_snv().map(|genot_snv| genot_snv.vals()).collect()
    }

    fn inner_vals_col_digit(&self, row_i: usize, digit_i: usize) -> Vec<bool> {
        self.check_arg(row_i, 0, digit_i);

        (0..self.col_n())
            .map(|ci| self.get_unchecked(row_i, ci, digit_i))
            .collect::<Vec<bool>>()
        /*
        // what is into_iter() for?
        (0..self.col_n())
            .into_iter()
            .map(|ci| self.get_unchecked(row_i, ci, digit_i))
            .collect::<Vec<bool>>()
             */
    }

    /// assume s=2
    fn inner_vals_col_digit_u8(&self, row_i: usize, digit_i: usize) -> Vec<u8> {
        self.check_arg_digit2();
        self.check_arg(row_i, 0, digit_i);

        (0..self.col_n())
            .into_iter()
            .map(|ci| self.get_unchecked(row_i, ci, digit_i) as u8)
            .collect::<Vec<u8>>()
    }

    fn inner_vals_col(&self, row_i: usize) -> Vec<u8> {
        self.check_arg_digit2();
        self.check_arg(row_i, 0, 0);

        (0..self.col_n())
            .into_iter()
            .map(|ci| self.get_val_unchecked(row_i, ci))
            .collect::<Vec<u8>>()
    }

    #[inline]
    fn inner_col_digit(&self, row_i: usize, digit_i: usize) -> &[B8] {
        self.check_arg(row_i, 0, digit_i);
        let len = self.col_len();
        let index_begin = self.index_row_digit_unchecked(row_i, digit_i);
        &self.inner()[index_begin..index_begin + len]
    }

    #[inline]
    fn inner_row(&self, row: usize) -> &[B8] {
        self.inner_rows(row, row + 1)
    }

    #[inline]
    fn inner_rows(&self, row_begin: usize, row_end: usize) -> &[B8] {
        self.check_arg_row_range(row_begin, row_end);
        let index_begin = self.index_row_unchecked(row_begin);
        let index_end = self.index_row_unchecked(row_end);
        &self.inner()[index_begin..index_end]
    }

    fn rows(&self, row_begin: usize, row_end: usize) -> CMatrixRef {
        CMatrixRef::new(
            self.inner_rows(row_begin, row_end),
            row_end - row_begin,
            self.col_n(),
            self.digit_n(),
        )
        /*
        let index_begin = self.index_row(row_begin);
        let index_end = self.index_row(row_end);
        let col_n = self.col_n();
        let digit = self.digit_n();
        CMatrixRef::new(
            &self.inner()[index_begin..index_end],
            row_end - row_begin,
            col_n,
            digit,
        )
         */
    }

    fn as_cmatrix_ref(&self) -> CMatrixRef {
        CMatrixRef::new(self.inner(), self.row_n(), self.col_n(), self.digit_n())
    }

    fn as_cvec_ref(&self, row_i: usize) -> CVecRef {
        CVecRef::new(
            self.inner_rows(row_i, row_i + 1),
            self.col_n(),
            self.digit_n(),
        )
    }

    fn iter_row(&self) -> RowIter {
        RowIter::new(self.as_cmatrix_ref())
    }

    /*
    fn iter_scal(&mut self) -> ScalIter {
        ScalIter {
            cmatrix: self.to_cmatrix_ref(),
            index_iter: ScalIndexIter::new(self.row_n(), self.col_n()),
        }
    }
     */

    fn check_new(&self) {
        // len
        if self.inner().len() != self.len() {
            panic!("wrong len");
        }
    }

    /// If some arguments are not applicable, put 0.
    fn check_arg(&self, row_i: usize, col_i: usize, digit_i: usize) {
        if row_i >= self.row_n() {
            panic!(
                "row_i {} should be smaller than row_n {}",
                row_i,
                self.row_n()
            );
        }
        if col_i >= self.col_n() {
            panic!(
                "col_i {} should be smaller than col_n {}",
                col_i,
                self.col_n()
            );
        }
        if digit_i >= self.digit_n() {
            panic!(
                "digit_i {} should be smaller than digit_n {}",
                digit_i,
                self.digit_n()
            );
        }
    }

    /// check row index range
    fn check_arg_row_range(&self, row_begin: usize, row_end: usize) {
        if row_end > self.row_n() {
            panic!(
                "row_end {} should be smaller than or the same as row_n {}",
                row_end,
                self.row_n()
            );
        }
        // allow begin=row_n, then get empty vec
        //if row_begin >= self.row_n() {
        if row_begin > self.row_n() {
            panic!(
                "row_begin {} should be smaller than row_n {}",
                row_begin,
                self.row_n()
            );
        }
        // allow begin=end, then get empty vec
        //if row_begin >= row_end {
        if row_begin > row_end {
            panic!(
                "row_begin {} should be smaller than row_end {}",
                row_begin, row_end,
            );
        }
    }

    /// Call when digit_n=2 only is allowed.
    fn check_arg_digit2(&self) {
        if self.digit_n() != 2 {
            panic!(
                "This function is allowed to use when digit_n()==2 but {}",
                self.digit_n()
            );
        }
    }
}

/// for CMatrix, CMatrixMut
pub trait BaseCMatrixMut: BaseCMatrix {
    fn inner_mut(&mut self) -> &mut [B8];

    #[inline]
    fn inner_row_mut(&mut self, row_i: usize) -> &mut [B8] {
        self.inner_rows_mut(row_i, row_i + 1)
    }

    #[inline]
    fn inner_rows_mut(&mut self, row_begin: usize, row_end: usize) -> &mut [B8] {
        self.check_arg_row_range(row_begin, row_end);
        let index_begin = self.index_row_unchecked(row_begin);
        let index_end = self.index_row_unchecked(row_end);
        &mut self.inner_mut()[index_begin..index_end]
    }

    #[inline]
    fn inner_col_digit_mut(&mut self, row_i: usize, digit_i: usize) -> &mut [B8] {
        self.check_arg(row_i, 0, digit_i);
        let len = self.col_len();
        let index_begin = self.index_row_digit_unchecked(row_i, digit_i);
        &mut self.inner_mut()[index_begin..index_begin + len]
    }

    /*     // assume row_n=1
    #[inline]
    fn inner_col_digit_row0_mut(&mut self, digit_i: usize) -> &mut [B8] {
        let len = self.col_len();
        &mut self.inner_mut()[digit_i * len..(digit_i + 1) * len]
    } */

    /*     #[inline]
    fn set_bool_row0_unchecked(&mut self, b: bool, col_i: usize, digit_i: usize) {
        match b {
            true => {
                self.inner_col_digit_row0_mut(digit_i)[col_i >> 3] |= 0x80 >> (col_i & 7);
            }
            false => {
                self.inner_col_digit_row0_mut(digit_i)[col_i >> 3] &= !(0x80 >> (col_i & 7));
            }
        }
    } */

    #[inline]
    fn set_bool_unchecked(&mut self, b: bool, row_i: usize, col_i: usize, digit_i: usize) {
        match b {
            true => {
                self.inner_col_digit_mut(row_i, digit_i)[col_i >> 3] |= 0x80 >> (col_i & 7);
            }
            false => {
                self.inner_col_digit_mut(row_i, digit_i)[col_i >> 3] &= !(0x80 >> (col_i & 7));
            }
        }
    }

    #[inline]
    fn set_bool(&mut self, b: bool, row_i: usize, col_i: usize, digit_i: usize) {
        self.check_arg(row_i, col_i, digit_i);

        self.set_bool_unchecked(b, row_i, col_i, digit_i);
    }

    /*     #[inline]
    fn set_row0_unchecked(&mut self, val: u8, col_i: usize) {
        self.set_bool_row0_unchecked((val >> 1) != (val & 1), col_i, 0);
        self.set_bool_row0_unchecked((val & 2) == 2, col_i, 1);
    } */

    /// now val -> (bool, bool) is 2->(1,1), 1->(1,0), 0->(0,0), m->(0,1)
    // TODO: genotype encoding should be 2=(1,1), 1=(1,0), 0=(0,1), m=(0,0) -> since (0,0) could be mixed with padding
    #[inline]
    fn set_unchecked(&mut self, val: u8, row_i: usize, col_i: usize) {
        self.set_bool_unchecked((val >> 1) != (val & 1), row_i, col_i, 0);
        self.set_bool_unchecked((val & 2) == 2, row_i, col_i, 1);
    }

    /*     #[inline]
    fn set_row0_unchecked(&mut self, val: u8, col_i: usize) {
        self.set_unchecked(val, 0, col_i)
    } */

    #[inline]
    fn set(&mut self, val: u8, row_i: usize, col_i: usize) {
        self.check_arg(row_i, col_i, 0);
        self.check_arg_set(val);
        self.set_unchecked(val, row_i, col_i)
    }

    /// If some arguments are not applicable, put 0.
    fn check_arg_set(&self, val: u8) {
        if val >= 4 {
            panic!("value {} should be smaller than 4", val);
        }
    }

    #[inline]
    fn set_by_bool_unchecked(&mut self, b0: bool, b1: bool, row_i: usize, col_i: usize) {
        self.set_bool_unchecked(b0, row_i, col_i, 0);
        self.set_bool_unchecked(b1, row_i, col_i, 1);
    }

    /// set bool assuming original is false
    /// do nth if b is false, should be fast
    #[inline]
    fn set_by_bool_init_unchecked(&mut self, b0: bool, b1: bool, row_i: usize, col_i: usize) {
        if b0 == true {
            self.set_bool_unchecked(b0, row_i, col_i, 0);
        }
        if b1 == true {
            self.set_bool_unchecked(b1, row_i, col_i, 1);
        }
    }

    /*     // TODO: check val<4
    /// assume row=1
    fn set_col_unchecked(&mut self, vec: &[u8]) {
        if self.row_n() != 1 {
            panic!("");
        }
        if self.digit_n() != 2 {
            panic!("");
        }
        vec.iter()
            .enumerate()
            .for_each(|(ci, x)| self.set_v(*x, ci));
        /*
        for ci in 0..self.col_n() {
            self.set(vec[ci], ci);
            /*
            let count = vec[ci];
            self.set_bool((count >> 1) != (count & 1), ci, 0);
            self.set_bool((count & 2) == 2, ci, 1);
            */
        }
         */
    } */

    // TODO: check val<4
    // use set() since rowi and coli is guaranteed
    fn set_vec_unchecked(&mut self, vec: &[u8]) {
        if self.digit_n() != 2 {
            panic!("");
        }
        let col_n = self.col_n();
        for row_i in 0..self.row_n() {
            vec[row_i * col_n..(row_i + 1) * col_n]
                .iter()
                .enumerate()
                .for_each(|(ci, x)| self.set(*x, row_i, ci));
            /*
            let mut row = self.rows_mut(row_i, row_i + 1);
            row.set_col_unchecked_v(&vec[row_i * col_n..(row_i + 1) * col_n]);
             */
        }
    }

    /*     // TODO: check val<4
    // DONE improve: not beautiful to use set_col() in BaseCVec
    fn set_vec(&mut self, vec: &[u8]) {
        if self.digit_n() != 2 {
            panic!("");
        }
        let col_n = self.col_n();
        for row_i in 0..self.row_n() {
            vec[row_i * col_n..(row_i + 1) * col_n]
                .iter()
                .enumerate()
                .for_each(|(ci, x)| self.set(*x, row_i, ci));
            /*
            let mut row = self.rows_mut(row_i, row_i + 1);
            row.set_col_unchecked_v(&vec[row_i * col_n..(row_i + 1) * col_n]);
             */
        }
    } */

    /*     fn set_col_bool(&mut self, vec: &[bool]) {
        if self.row_n() != 1 {
            panic!("");
        }
        if self.digit_n() != 1 {
            panic!("");
        }
        for ci in 0..self.col_n() {
            let count = vec[ci];
            self.set_bool_row0_unchecked(count, ci, 0);
            //self.set_val((count >> 1) != (count & 1), ci, 0);
            //self.set_val((count & 2) == 2, ci, 1);
        }
    } */

    fn rows_mut(&mut self, row_begin: usize, row_end: usize) -> CMatrixMut {
        /*
        let index_begin = self.index_row_unchecked(row_begin);
        let index_end = self.index_row_unchecked(row_end);
         */
        //let col_n = self.col_n();
        //let digit = self.digit_n();
        let (col_n, digit_n) = (self.col_n(), self.digit_n());
        CMatrixMut::new(
            self.inner_rows_mut(row_begin, row_end),
            row_end - row_begin,
            col_n,
            digit_n,
        )
    }

    // DONE not beautiful to use set_col() in BaseCVec
    fn set_vec_bool(&mut self, vec: &[bool]) {
        if self.digit_n() != 1 {
            panic!("");
        }
        let col_n = self.col_n();
        for row_i in 0..self.row_n() {
            vec[row_i * col_n..(row_i + 1) * col_n]
                .iter()
                .enumerate()
                .for_each(|(ci, x)| self.set_bool_unchecked(*x, row_i, ci, 0));

            //let mut row = self.rows_mut(row_i, row_i + 1);
            //row.set_col_bool(&vec[row_i * col_n..(row_i + 1) * col_n]);
        }
    }

    fn as_cmatrix_mut(&mut self) -> CMatrixMut {
        let (row_n, col_n, digit_n) = (self.row_n(), self.col_n(), self.digit_n());
        CMatrixMut::new(self.inner_mut(), row_n, col_n, digit_n)
    }

    fn as_cvec_mut(&mut self, row_i: usize) -> CVecMut {
        let (col_n, digit_n) = (self.col_n(), self.digit_n());
        CVecMut::new(self.inner_row_mut(row_i), col_n, digit_n)
        //CVecMut::new(self.inner_rows_mut(row_i, row_i + 1), col_n, digit_n)
    }

    fn iter_row_mut(&mut self) -> RowIterMut {
        RowIterMut::new(self.as_cmatrix_mut())
    }

    /*
    // somehow error in next()
    fn rows_iter_mut(&mut self, row_begin: usize, row_end: usize) -> CMatrixMut {
        let index_begin = self.index_row(row_begin);
        let index_end = self.index_row(row_end);
        let col_n = self.col_n();
        let digit = self.digit();

        unsafe {
            let s = &mut self.inner_mut()[index_begin..index_end];
            let ptr = s.as_mut_ptr();
            let p = std::slice::from_raw_parts_mut(ptr, 2);

            CMatrixMut {
                inner: p,
                row_n: row_end - row_begin,
                col_n,
                digit,
            }
        }
    }
     */
}
