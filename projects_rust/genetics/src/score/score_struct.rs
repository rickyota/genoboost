use std::collections::HashMap;
use std::path::Path;

use crate::alloc;
use crate::genot;
use crate::genot_calc::Table3By3Ar;
use crate::textfile;

#[derive(Default, Debug)]
pub struct SampleScore {
    n: usize,
    // unncessary; use scores_pad.len()
    //len_n: usize,
    scores_pad: Vec<f64>,
}

//// padding <-> supress
//fn supress_len(v: &[f64], n: usize) -> &[f64] {
//    v.map(|x| &x[..n])
//}
//fn supress_len_mut(v: Option<&mut [f64]>, n: usize) -> Option<&mut [f64]> {
//    //v.map(|x| &x[..self.n])
//    //let ws_mut=self.ws.as_deref_mut();
//    v.map(|x| &mut x[..n])
//}

impl SampleScore {
    pub fn new(n: usize) -> Self {
        let len_n = n + 32;
        let mut scores_pad = alloc::with_capacity_align::<f64>(n + 32);
        scores_pad.resize(len_n, 0.0f64);
        Self { n, scores_pad }
    }

    pub fn new_vec(v: Vec<f64>) -> Self {
        let n = v.len();
        let len_n = n + 32;
        let scores_pad = alloc::vec_align_f64(v, len_n);

        Self { n, scores_pad }
    }

    pub fn new_empty() -> Self {
        Self {
            n: 0,
            scores_pad: alloc::with_capacity_align::<f64>(0),
        }
    }

    pub fn new_from(
        file: &Path,
        // TODO
        //file_buf: &[u8],
        samples_name: &[String],
        samples_name_val: Option<&[String]>,
    ) -> (Self, Self) {
        let score_buf = if file.extension().unwrap() == "zst" {
            log::debug!("zst: {:?}", file);
            textfile::read_file_to_end(file, Some("zst")).unwrap_or_else(|_| {
                panic!("Cannot read file: {:?}", file);
            })
        } else {
            textfile::read_file_to_end(file, None).unwrap_or_else(|_| {
                panic!("Cannot read file: {:?}", file);
            })
        };

        Self::new_from_vec(&score_buf, samples_name, samples_name_val)

        //let cols = [0usize, 1];
        //let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(&score_buf, &cols, true);

        //let mut name_score = HashMap::with_capacity(vss[0].len());
        //for vi in 0..vss[0].len() {
        //    name_score.insert(vss[0][vi].clone(), vss[1][vi].parse::<f64>().unwrap());
        //}

        //// for training samples
        //let scores = samples_name
        //    .iter()
        //    .map(|x| name_score.get(x).unwrap().clone())
        //    .collect::<Vec<f64>>();

        //let scores = Self::new_vec(scores);

        //// for validation samples
        //let scores_val = if let Some(samples_name_val) = samples_name_val {
        //    let scores_val = samples_name_val
        //        .iter()
        //        .map(|x| name_score.get(x).unwrap().clone())
        //        .collect::<Vec<f64>>();

        //    Self::new_vec(scores_val)
        //} else {
        //    Self::new(0)
        //};

        //(scores, scores_val)
    }

    pub fn new_from_vec(
        //file: &Path,
        score_buf: &[u8],
        samples_name: &[String],
        samples_name_val: Option<&[String]>,
    ) -> (Self, Self) {
        //let score_buf = textfile::read_file_to_end(file, None).unwrap_or_else(|_| {
        //    panic!("Cannot read file: {:?}", file);
        //});

        let cols = [0usize, 1];
        let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(&score_buf, &cols, true);

        let mut name_score = HashMap::with_capacity(vss[0].len());
        for vi in 0..vss[0].len() {
            name_score.insert(vss[0][vi].clone(), vss[1][vi].parse::<f64>().unwrap());
        }

        // for training samples
        let scores = samples_name
            .iter()
            .map(|x| name_score.get(x).unwrap().clone())
            .collect::<Vec<f64>>();

        let scores = Self::new_vec(scores);

        // for validation samples
        let scores_val = if let Some(samples_name_val) = samples_name_val {
            let scores_val = samples_name_val
                .iter()
                .map(|x| name_score.get(x).unwrap().clone())
                .collect::<Vec<f64>>();

            Self::new_vec(scores_val)
        } else {
            Self::new(0)
        };

        (scores, scores_val)
    }

    pub fn n(&self) -> usize {
        self.n
    }

    // len_n >= n
    pub fn len_n(&self) -> usize {
        self.scores_pad().len()
        //self.len_n
    }

    pub fn scores_pad_mut(&mut self) -> &mut [f64] {
        &mut self.scores_pad
    }

    pub fn scores_mut(&mut self) -> &mut [f64] {
        let n = self.n();
        &mut self.scores_pad_mut()[..n]
    }

    pub fn scores_pad(&self) -> &[f64] {
        &self.scores_pad
    }

    pub fn scores(&self) -> &[f64] {
        &self.scores_pad()[..self.n()]
    }

    // make n..score.len() zero
    pub fn clear_pad(&mut self) {
        let n = self.n();
        self.scores_pad_mut()[n..]
            .iter_mut()
            .for_each(|x| *x = 0.0f64);
    }

    pub fn check_no_nan(&self) {
        self.scores().iter().enumerate().for_each(|(i, x)| {
            if x.is_nan() {
                panic!("score contains nan at index: {}", i);
            }
        });
    }

    pub fn clone_align(&self) -> Self {
        let scores_pad = self.scores_pad();
        let mut scores_pad_new = alloc::with_capacity_align::<f64>(scores_pad.len());
        scores_pad_new.resize(scores_pad.len(), 0.0f64);
        scores_pad_new.copy_from_slice(scores_pad);
        Self {
            n: self.n,
            scores_pad: scores_pad_new,
        }
    }

    // ignore for padding
    pub fn add_scores(&mut self, other: &Self) {
        self.scores_mut()
            .iter_mut()
            .zip(other.scores().iter())
            .for_each(|(x, y)| *x += y);
    }
}

// TODO: merge with Table3By3Ar? -> troublesome for Eq and 0.0f64?
// TODO: use trait
// ((x00, x01, x02),(x10,x11,x12),(x20,x21,x22)
#[derive(Debug, Clone, Copy)]
pub struct Sum3by3Ar([f64; 9]);

impl Sum3by3Ar {
    pub fn from_table(table: Table3By3Ar) -> Self {
        let mut x = [0.0f64; 9];
        let table_ar = table.inner();
        for i in 0..9 {
            x[i] = table_ar[i] as f64;
        }
        Sum3by3Ar(x)
    }

    pub fn sum(&self) -> f64 {
        self.0.iter().sum()
    }

    pub fn interaction_genotype(maf_1: f64, maf_2: f64) -> Self {
        if maf_1.is_nan() || maf_2.is_nan() {
            panic!("maf is nan.");
        }

        let x_1 = standardized_genotype_ar(maf_1);
        let x_2 = standardized_genotype_ar(maf_2);

        Self::vec_multi_matrix_ar(x_1, x_2)
    }

    pub fn add_scalar(&self, a: f64) -> Self {
        let mut res = [0.0f64; 9];
        for i in 0..9 {
            res[i] = self.0[i] + a;
        }
        Sum3by3Ar(res)
    }

    pub fn substract(&self, other: &Self) -> Self {
        let mut res = [0.0f64; 9];
        for i in 0..9 {
            res[i] = self.0[i] - other.0[i];
        }
        Sum3by3Ar(res)
    }

    pub fn multiply_scalar(&self, a: f64) -> Self {
        let mut res = [0.0f64; 9];
        for i in 0..9 {
            res[i] = self.0[i] * a;
        }
        Sum3by3Ar(res)
    }

    pub fn multiply(&self, other: &Self) -> Self {
        let mut res = [0.0f64; 9];
        for i in 0..9 {
            res[i] = self.0[i] * other.0[i];
        }
        Sum3by3Ar(res)
    }

    // self * other * other
    pub fn multiply_pow(&self, other: &Self) -> Self {
        let mut res = [0.0f64; 9];
        for i in 0..9 {
            res[i] = self.0[i] * other.0[i] * other.0[i];
        }
        Sum3by3Ar(res)
    }

    pub fn linearconst(&self, c: f64, a: f64) -> Self {
        let mut res = [0.0f64; 9];
        for i in 0..9 {
            res[i] = c + a * self.0[i];
        }
        Sum3by3Ar(res)
    }

    pub fn vec_multi_matrix_ar(x: [f64; 3], y: [f64; 3]) -> Self {
        // using "for" should be slower?
        Self([
            x[0] * y[0],
            x[0] * y[1],
            x[0] * y[2],
            x[1] * y[0],
            x[1] * y[1],
            x[1] * y[2],
            x[2] * y[0],
            x[2] * y[1],
            x[2] * y[2],
        ])
    }
}

pub fn standardized_genotype_ar(maf: f64) -> [f64; 3] {
    [-2.0 * maf, 1.0 - 2.0 * maf, 2.0 - 2.0 * maf]
}

// TODO: use Sum3by3([f64;9])
// For g1, g2
// ((x00, x01, x02),(x10,x11,x12),(x20,x21,x22))
#[derive(Debug, Clone, Copy)]
pub struct Sum3by3(((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)));
//type Sum3by3 = ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64));

impl Sum3by3 {
    pub fn new(x: ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64))) -> Self {
        Sum3by3(x)
    }

    pub fn to_tuple(self) -> ((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)) {
        let Sum3by3(x) = self;
        x
    }

    pub fn linearconst(self, c: f64, a: f64) -> Sum3by3 {
        let Sum3by3((x0, x1, x2)) = self;

        fn linearconst_3(x: (f64, f64, f64), c: f64, a: f64) -> (f64, f64, f64) {
            let (x0, x1, x2) = x;
            (c + a * x0, c + a * x1, c + a * x2)
        }

        Sum3by3((
            linearconst_3(x0, c, a),
            linearconst_3(x1, c, a),
            linearconst_3(x2, c, a),
        ))
    }

    pub fn sum(x: Self) -> f64 {
        let Sum3by3((x0, x1, x2)) = x;
        fn sum_3(a: (f64, f64, f64)) -> f64 {
            a.0 + a.1 + a.2
        }
        sum_3(x0) + sum_3(x1) + sum_3(x2)
    }

    pub fn sum_multi(x: Self, y: Self) -> f64 {
        let Sum3by3((x0, x1, x2)) = x;
        let Sum3by3((y0, y1, y2)) = y;
        fn sum_mult_3(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
            a.0 * b.0 + a.1 * b.1 + a.2 * b.2
        }
        sum_mult_3(x0, y0) + sum_mult_3(x1, y1) + sum_mult_3(x2, y2)
    }

    /// sum(x*y^2)
    pub fn sum_multi_pow(x: Self, y: Self) -> f64 {
        let Sum3by3((x0, x1, x2)) = x;
        let Sum3by3((y0, y1, y2)) = y;
        fn sum_mult_3(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
            a.0 * b.0 * b.0 + a.1 * b.1 * b.1 + a.2 * b.2 * b.2
        }
        sum_mult_3(x0, y0) + sum_mult_3(x1, y1) + sum_mult_3(x2, y2)
    }

    // TODO: clean
    pub fn check_no_nan(&self) {
        let Sum3by3((x0, x1, x2)) = *self;
        let x = vec![x0, x1, x2];
        x.iter().for_each(|(x_0, x_1, x_2)| {
            if x_0.is_nan() || x_1.is_nan() || x_2.is_nan() {
                panic!("score contains nan.")
            }
        });
    }
}

fn vec_multi_matrix(x: (f64, f64, f64), y: (f64, f64, f64)) -> Sum3by3 {
    let (x0, x1, x2) = x;
    let (y0, y1, y2) = y;
    Sum3by3::new((
        (x0 * y0, x0 * y1, x0 * y2),
        (x1 * y0, x1 * y1, x1 * y2),
        (x2 * y0, x2 * y1, x2 * y2),
    ))
}

pub fn standardized_genotype(maf: f64) -> (f64, f64, f64) {
    (-2.0 * maf, 1.0 - 2.0 * maf, 2.0 - 2.0 * maf)
}

pub fn interaction_genotype(maf_1: f64, maf_2: f64) -> Sum3by3 {
    if maf_1.is_nan() || maf_2.is_nan() {
        panic!("maf is nan.");
    }

    let x_1 = standardized_genotype(maf_1);
    let x_2 = standardized_genotype(maf_2);

    vec_multi_matrix(x_1, x_2)
}

// TODO: test
// For g1, g2
// ((x00, x01, x02, x03),(x10,x11,x12,x13),(x20,x21,x22,x23),(x30,x31,x32,x33))
#[derive(Debug, Clone, Copy)]
pub struct Sum4by4([f64; 16]);
//pub struct Sum4by4(((f64, f64, f64), (f64, f64, f64), (f64, f64, f64)));

impl Sum4by4 {
    pub fn val(&self, g1: usize, g2: usize) -> f64 {
        self.0[g1 * 4 + g2]
    }

    pub fn to_4x4(self) -> [[f64; 4]; 4] {
        let x = self.0;
        //let Sum4by4(x) = self;
        [
            [x[0], x[1], x[2], x[3]],
            [x[4], x[5], x[6], x[7]],
            [x[8], x[9], x[10], x[11]],
            [x[12], x[13], x[14], x[15]],
        ]
    }

    pub fn linearconst(&self, c: f64, a: f64) -> Sum4by4 {
        let x = self.0;
        let mut res = [0.0f64; 16];
        for i in 0..16 {
            res[i] = c + a * x[i];
        }
        Sum4by4(res)
    }
}

fn mode_to_score(mode: u8, scores: (f64, f64, f64)) -> f64 {
    match mode {
        0 => scores.0,
        1 => scores.1,
        2 => scores.2,
        _ => panic!("Wrong genotype."),
    }
}

pub fn standardized_genotype_missing(
    maf: f64,
    missing_to_mode: bool,
    missing_to_mean: bool,
) -> [f64; 4] {
    let (sg0, sg1, sg2) = (-2.0 * maf, 1.0 - 2.0 * maf, 2.0 - 2.0 * maf);

    // TODO: fn
    let sgm = if missing_to_mode {
        let mode = genot::maf_to_mode(maf);
        mode_to_score(mode, (sg0, sg1, sg2))
    } else if missing_to_mean {
        0.0
        //const_ti + 2.0 * maf.unwrap() * alpha_ti
    } else {
        f64::NAN
    };

    [sg0, sg1, sg2, sgm]
}

fn vec_multi_matrix_missing(x: [f64; 4], y: [f64; 4]) -> Sum4by4 {
    let mut res = [0.0f64; 16];
    for i in 0..4 {
        for j in 0..4 {
            res[i * 4 + j] = x[i] * y[j];
        }
    }
    Sum4by4(res)
}

// mv in Self
pub fn interaction_genotype_missing(
    maf_1: f64,
    maf_2: f64,
    missing_to_mode: bool,
    missing_to_mean: bool,
) -> Sum4by4 {
    if maf_1.is_nan() || maf_2.is_nan() {
        panic!("maf is nan.");
    }

    let x_1 = standardized_genotype_missing(maf_1, missing_to_mode, missing_to_mean);
    let x_2 = standardized_genotype_missing(maf_2, missing_to_mode, missing_to_mean);

    vec_multi_matrix_missing(x_1, x_2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clear_pad() {
        let mut samples = SampleScore::new(10);
        let scores_pad = samples.scores_pad_mut();
        scores_pad.iter_mut().for_each(|x| *x = 1.0f64);
        samples.clear_pad();
        let len_n = samples.len_n();

        let scores_pad = samples.scores_pad_mut();
        for i in 0..10 {
            assert_eq!(scores_pad[i], 1.0f64);
        }
        assert_eq!(len_n, 42);
        for i in 10..len_n {
            assert_eq!(scores_pad[i], 0.0f64);
        }
    }
}
