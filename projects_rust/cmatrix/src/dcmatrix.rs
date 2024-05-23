//! Diff Compressed Matrix
//!
//! In addition to Compressed Matrix, offer
//! 1. storing only diff if there are similar vector in other rows
//! 2. storing sample index only for rare vector

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct DCMatrix {
    cmatrix: CMatrix,
    flag_diff: CBits,
    // row_i -> ref_row_i
    reference_diff: HashMap<usize, usize>,
    // row_i -> (ref_row_i, ind_diff_s0 index)
    //reference_diff: HashMap<usize,(usize,usize)>,
    // sample of diff from ori
    // use HashMap for snv id for random access
    // use Vec for sample id for sequencial access
    // slow using HashMap??
    ind_diff_s0: HashMap<usize, Vec<u32>>,
    ind_diff_s1: HashMap<usize, Vec<u32>>,
    // put all in one vec
    //ind_diff_s0: Vec<u32>,
    //ind_diff_s1: Vec<u32>,
    // to implement in calc_loss
    //reference_empty: Vec<u8>,
    flag_rare: CBits,
    ind_rare_s0: HashMap<Vec<u32>>,
    ind_rare_s1: HashMap<Vec<u32>>,
    //ind_rare_s0: Vec<u32>,
    //ind_rare_s1: Vec<u32>,
}

impl DCMatrix {
    // impossible to return &[B8]
    fn inner_col_digit_vec(&self, row_i: usize, digit_i: usize) -> Vec<u8> {
        //self.cmatrix.inner_col_digit(0, digit_i)
        if flag_diff.get(row_i) {
            // use diff
            let ref_row_i = reference_diff.get(row_i).unwrap();
            let ref_vec = self.cmatrix.inner_col_digit(ref_row_i, digit_i).clone();
            let ind_diff_s0_row = self.ind_diff_s0.get(row_i).unwrap();
            let ind_diff_s1_row = self.ind_diff_s1.get(row_i).unwrap();

            ref_vec.set_s0(ind_diff_s0_row);
            ref_vec.set_s1(ind_diff_s1_row);
            ref_vec
        } else if flag_rare.get(row_i) {
            // use rare
            // TODO: or you do not need this?
            // implement calc_loss for ind_rare
            // -> if you implement this, implementing flag_diff is easy as well.
            // -> No, for rare variant, iterate over non-zero samples is better.
        } else {
            // clone cmatrix
        }
    }

    // use  &[88] if not diff
    fn inner_col_digit_vec_diff(&self, row_i: usize, digit_i: usize) -> Option<Vec<u8>> {
        //self.cmatrix.inner_col_digit(0, digit_i)
        if flag_diff.get(row_i) || flag_rare.get(row_i) {
            Some(self.inner_col_digit_vec(row_i, digit_i))
        } else {
            None
        }
    }
}
