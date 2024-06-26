use std::collections::HashMap;

use crate::sample;
use crate::textfile;
use crate::vec;
use crate::CovKind;
use crate::GenotFile;
//use crate::{CovKind, Cov};
//use super::prelude::*;

// usually accompanied with Samples
// should implement these in Samples?
pub trait CovsTrait {
    fn covs(&self) -> Option<&Covs>;
    fn vals(&self) -> Option<&Vec<Vec<f64>>> {
        self.covs().unwrap().vals()
    }
    fn cov_indexs(&self) -> Option<&Vec<CovId>> {
        self.covs().unwrap().cov_indexs()
    }
    fn covs_n(&self) -> usize {
        match self.covs() {
            None => 0,
            Some(x) => x.cov_indexs().unwrap().len(),
        }
    }
}

// not Vec<Cov> since vals should be matrix
// TODO: Covs must not be None since "const" should be in
// -> or "const " should not be here?? since "const" is not covariates
#[derive(Clone)]
pub struct Covs {
    cov_indexs: Vec<CovId>,
    // cov_n x sample_n
    // col-major
    vals: Vec<Vec<f64>>,
}

impl CovsTrait for Covs {
    // or return &Covs
    // Here, you should return Covs since if you have Covs::covs(), covs should exist.
    // but from Samples, you should implement as Option<&Covs>, since might not exist
    fn covs(&self) -> Option<&Covs> {
        Some(self)
    }
    fn vals(&self) -> Option<&Vec<Vec<f64>>> {
        Some(&self.vals)
    }
    fn cov_indexs(&self) -> Option<&Vec<CovId>> {
        Some(&self.cov_indexs)
    }
}

impl Covs {
    //pub fn new_data_from_vars(vars: Vec<Cov>) -> Covs {
    //    let cov_indexs = (&vars)
    //        .iter()
    //        .map(|x| CovId::new(x.name().to_owned(), x.kind()))
    //        .collect();
    //    let vals = vars.into_iter().map(|x| x.vals_consume()).collect();
    //    Covs { cov_indexs, vals }
    //}

    /// if cov_name="" return Covs with 0 covariates
    //pub fn new(
    //    cov_buf: Option<&[u8]>,
    //    fin_genot: &GenotFile,
    //    cov_name: &str,
    //    //n_in: usize,
    //    sample_id_to_n: &HashMap<String, usize>,
    //    //use_samples: &[bool],
    //) -> Covs {
    //    Self::new_buf(cov_buf, fin_genot, cov_name, sample_id_to_n)
    //}

    /*
    // TODO: deprecate
    fn load_from_file(
        fin_var: Option<&Path>,
        cov_name: &str,
        sample_id_to_n: &HashMap<String, usize>,
    ) -> Covs {
        //let cols_name = textfile::load_table_header(fin_var.unwrap()).unwrap();
        //let cov_name_v = parse_cov_name_in_cols(cov_name, &cols_name);

        let buf = textfile::read_file_to_end(fin_var.unwrap(), None).unwrap();

        Self::load_from_file_buf(&buf[..], cov_name, sample_id_to_n)

        /*
        let valss = textfile::load_table(fin_var.unwrap(), true).unwrap();

        let n = sample_id_to_n.keys().len();

        let mut covs_indexs: Vec<CovId> = Vec::with_capacity(n);
        let mut covs_vals: Vec<Vec<f64>> = vec![];
        //let mut vars: Vec<Var> = Vec::with_capacity(n);

        //let col_fid = 0;
        // col=0 could be fid or iid but both are fine
        let col_iid = 0;

        for pi in 0..cols_name.len() {
            if pi == col_iid {
                continue;
            }

            let col_name = cols_name[pi].clone();
            if !cov_name_v.contains(&col_name) {
                continue;
            }

            let cov_index = CovId::new_cov(col_name);
            //let cov_index = CovId::new_cov(cols_name[pi].clone());

            //let mut var = Var::construct_var(kind, cols_name[pi].clone());

            let vals = cov::vals_align_id(
                &valss[pi],
                &valss[col_iid],
                //(&valss[col_fid], &valss[col_iid]),
                sample_id_to_n,
                n,
            );

            //var.set_vals(vals);
            //vars.push(var);

            covs_indexs.push(cov_index);
            covs_vals.push(vals);
        }
        let covs = Covs {
            cov_indexs: covs_indexs,
            vals: covs_vals,
        };
        //Some(covs)
        covs */
    } */

    /// if cov_name="" return Covs with 0 covariates
    pub fn new_buf(
        cov_buf: Option<&[u8]>,
        fin_genot: &GenotFile,
        cov_name: &str,
        sample_id_to_n: &HashMap<String, usize>,
    ) -> Self {
        //TODO: use csv

        let cov_buf_v;
        let cov_buf = match cov_buf {
            Some(x) => x,
            None => match fin_genot {
                GenotFile::Plink1(_) => {
                    panic!("Do not use --phe-name for plink1.");
                }
                GenotFile::Plink2(_) | GenotFile::Plink2Vzs(_) => {
                    let fin_fam = fin_genot.sample_file_exist_chrom().unwrap();
                    //let fin_fam = genot_io::fname_fam_exist_chrom(fin, gfmt).unwrap();
                    cov_buf_v = textfile::read_file_to_end(&fin_fam, None).unwrap();
                    &cov_buf_v[..]
                }
            },
        };

        let valss = textfile::load_table_buf(cov_buf, true);
        let cols_name = textfile::load_table_header_buf(cov_buf);
        let cov_name_v = parse_cov_name_in_cols(cov_name, &cols_name);

        let n = sample_id_to_n.keys().len();

        let mut covs_indexs: Vec<CovId> = Vec::with_capacity(n);
        let mut covs_vals: Vec<Vec<f64>> = vec![];
        //let mut vars: Vec<Var> = Vec::with_capacity(n);

        //let col_fid = 0;
        // col=0 could be fid or iid but both are fine
        let col_iid = 0;

        if !vec::has_unique_elements(&valss[col_iid]) {
            panic!("First column in --file-phe should be unique.")
        };

        for pi in 0..cols_name.len() {
            if pi == col_iid {
                continue;
            }

            let col_name = cols_name[pi].clone();
            if !cov_name_v.contains(&col_name) {
                continue;
            }
            log::debug!("pi, col_name: {}, {}", pi, col_name);

            let cov_index = CovId::new_cov(col_name);
            //let cov_index = CovId::new_cov(cols_name[pi].clone());

            let vals =
                sample::vals_align_id_type::<f64>(&valss[pi], &valss[col_iid], sample_id_to_n);

            //let vals = cov::vals_align_id_f64(
            //    &valss[pi],
            //    &valss[col_iid],
            //    //(&valss[col_fid], &valss[col_iid]),
            //    sample_id_to_n,
            //);

            //var.set_vals(vals);
            //vars.push(var);

            covs_indexs.push(cov_index);
            covs_vals.push(vals);
        }

        let covs = Self {
            cov_indexs: covs_indexs,
            vals: covs_vals,
        };
        //Some(covs)
        covs
    }

    // for test
    //pub fn new_vals(vals: Vec<Vec<f64>>) -> Self {
    //    Self {
    //        cov_indexs: vec![CovId::default(); vals.len()],
    //        vals,
    //    }
    //}

    /// convert form row major to column major and return
    /// samples x var
    /// mainly for smartcore::DenseMatrix
    pub fn vals_row_major(&self) -> Vec<Vec<f64>> {
        let vals = self.vals().unwrap();

        let vals_col = vec::convert_vec2d_to_row_major(vals);

        //let samples_n = vals[0].len();
        ////let covs_n = self.covs_n();
        //let mut vals_col = Vec::new();
        //for _ in 0..samples_n {
        //    vals_col.push(Vec::new());
        //}
        //for v in vals.iter() {
        //    for (x, v_row) in v.iter().zip(vals_col.iter_mut()) {
        //        v_row.push(*x);
        //    }
        //}

        log::debug!(
            "row_major: sample {} x var {}",
            vals_col.len(),
            vals_col[0].len()
        );
        vals_col
    }

    /// convert form col major to row major vec
    /// mainly for smartcore::DenseMatrix, linfa and Array2D
    //pub fn vals_row_major_vec(&self) -> (Vec<f64>, usize, usize) {
    //    let vals = self.vals().unwrap();
    //    let samples_n = vals[0].len();
    //    let covs_n = self.covs_n();

    //    let mut vals_col = vec![0.0f64; samples_n * covs_n];
    //    // col_v: vec of samples
    //    for (ci, col_v) in vals.iter().enumerate() {
    //        for (ri, x) in col_v.iter().enumerate() {
    //            vals_col[ci + ri * covs_n] = *x
    //        }
    //    }
    //    (vals_col, samples_n, covs_n)
    //}

    // uncomment when use
    /*     ///
    pub fn vals_column_major_vec(&self) -> (Vec<f64>, usize, usize) {
        let vals = self.vals().unwrap();
        let samples_n = vals[0].len();
        let covs_n = self.covs_n();

        let mut vals_col = vec![0.0f64; samples_n * covs_n];
        // col_v: vec of samples
        for (ci, col_v) in vals.iter().enumerate() {
            for (ri, x) in col_v.iter().enumerate() {
                vals_col[ri + ci * samples_n] = *x
            }
        }
        (vals_col, samples_n, covs_n)
    } */

    pub fn vals_row_major_norm(&self) -> Vec<Vec<f64>> {
        vec::convert_vec2d_to_row_major(&vec::norm_vec2d(self.vals().unwrap().clone()))

        //let mut vals = self.vals_row_major();
        //let means = self.covs().unwrap().means();
        //let stds = self.covs().unwrap().stds();
        //for row_v in vals.iter_mut() {
        //    for (ci, x) in row_v.iter_mut().enumerate() {
        //        let mean = means[ci];
        //        let var = stds[ci];
        //        *x = (*x - mean) / var;
        //    }
        //}
        //vals
    }

    // for linfa
    // (vec1d, (row, col))
    // row=sample_n
    // col=cov_n
    pub fn vals_norm_vec1d_row_major(&self) -> (Vec<f64>, (usize, usize)) {
        let v_norm = self.vals_row_major_norm();
        vec::convert_vec2d_to_vec1(v_norm)
    }

    pub fn vals_id(&self, name: &str) -> &[f64] {
        /*
        log::debug!("name {}", name);
        for x in self.cov_indexs().unwrap().iter() {
            log::debug!("x {}", x.name())
        }
         */

        //log::debug!("name {:?}", name);
        //log::debug!("cov_indes {:?}", self.cov_indexs());

        let index = self
            .cov_indexs()
            .unwrap()
            .iter()
            .position(|x| x.name() == name)
            .unwrap_or_else(|| panic!("name is not in Covs: {}", name));
        //.unwrap();
        &self.vals().unwrap()[index]
    }

    pub fn means(&self) -> Vec<f64> {
        vec::mean_v2(self.vals().unwrap())

        //let vals = self.vals().unwrap();
        //let mut means = Vec::new();
        //for v in vals.iter() {
        //    let mean = v.iter().sum::<f64>() / (v.len() as f64);
        //    means.push(mean);
        //}
        //means
    }

    pub fn stds(&self) -> Vec<f64> {
        vec::std_v2(self.vals().unwrap())
        //let vals = self.vals().unwrap();
        //let means = self.means();
        //let mut stds = Vec::new();
        //for (v, mean) in vals.iter().zip(means.iter()) {
        //    let sq = v.iter().map(|x| x * x).sum::<f64>() / (v.len() as f64);
        //    let var = sq - mean * mean;
        //    let std = var.sqrt();
        //    stds.push(std);
        //}
        //stds
    }
}

// mainly for Wgt
#[derive(Debug, Clone, Default)]
pub struct CovId {
    // const's name must be "const"
    name: String,
    kind: CovKind,
}

impl CovId {
    pub fn new_const() -> CovId {
        CovId {
            name: "const".to_owned(),
            kind: CovKind::Const,
        }
    }
    pub fn new_cov(name: String) -> CovId {
        CovId {
            name,
            kind: CovKind::Cov,
        }
    }
    // should not use since accessing CovKind is unnecessary
    pub fn new(name: String, kind: CovKind) -> CovId {
        CovId { name, kind }
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn kind(&self) -> CovKind {
        self.kind
    }
}

fn parse_cov_name_in_cols(cov_name: &str, cols: &[String]) -> Vec<String> {
    let cov_name_v: Vec<&str> = cov_name.split(',').collect();

    fn index_in_col(x: &str, cols: &[String]) -> usize {
        cols.iter()
            .position(|y| y == x)
            .unwrap_or_else(|| panic!("Cov_name {} is not in the column of fin_phe", x))
    }

    let mut cov_name_in_col: Vec<String> = vec![];
    for x in cov_name_v.iter() {
        if x.contains('-') {
            // add all the cols in the range
            let x_split: Vec<&str> = x.split('-').collect();
            if x_split.len() != 2 {
                panic!("Format of cov_name is wrong on {}", x);
            }
            let x1 = x_split[0];
            let x2 = x_split[1];
            let x1_idx = index_in_col(x1, cols);
            let x2_idx = index_in_col(x2, cols);
            for idx in x1_idx..=x2_idx {
                cov_name_in_col.push(cols[idx].clone());
            }
        } else {
            if cols.contains(&x.to_string()) {
                cov_name_in_col.push(x.to_string());
            } else {
                panic!("Cov_name {} is not in the column of fin_phe", x)
            }
        }
    }
    cov_name_in_col
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cov_name_in_cols() {
        let cols = vec!["IID", "sex", "cov", "PC1", "PC2", "PC3", "t2d"];
        let cols: Vec<String> = cols.iter().map(|x| x.to_string()).collect();
        let cov_name = "cov,sex,PC1-PC3";

        let cov_use = parse_cov_name_in_cols(cov_name, &cols);

        // orders follow cov_name
        assert_eq!(cov_use, vec!["cov", "sex", "PC1", "PC2", "PC3"]);
    }
}
