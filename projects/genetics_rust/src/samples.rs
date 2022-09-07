//! Dataset
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait

pub mod base_phe;
pub mod iterator;
pub mod phe;
pub mod prelude;

use std::collections::HashMap;
use std::path::Path;

use crate::cov;
use crate::{plink, text};
use crate::{CovKind, Var, B8};
use cmatrix::BaseCMatrix;
use prelude::*;

// maybe change to (String, String)
type Name = (String, String);

// SampleTrait: CovsTrait; using subtrait is troublesome
// indicate both

#[derive(Clone)]
pub struct Samples {
    // change to Labels for several phenotypes
    phe: Phe,
    names: Option<Vec<Name>>,
    // move to DatasetBiNew? -> put here for a while because covs also needs names
    // also, Samples only can be input of fn
    covs: Option<Covs>,
    // fields: Option<SamplesFields> // input
    //qc: Option<SamplesQc>,
    samples_n: usize,
}

impl Samples {
    pub fn new(
        phe: Vec<bool>,
        names: Option<Vec<Name>>,
        covs: Option<Covs>,
        samples_n: usize,
    ) -> Samples {
        let phe_struct = Phe::new(phe);
        Samples {
            phe: phe_struct,
            names,
            covs,
            samples_n,
        }
    }
    pub fn new_data(phe: Vec<bool>, names: Option<Vec<Name>>, samples_n: usize) -> Samples {
        let phe_struct = Phe::new(phe);
        Samples {
            phe: phe_struct,
            names,
            covs: None,
            samples_n,
        }
    }
    pub fn new_empty() -> Samples {
        Samples {
            phe: Phe::new_empty(0),
            names: None,
            covs: None,
            samples_n: 0,
        }
    }

    pub fn add_covs_with_vars(self, vars: Vec<Var>) -> Samples {
        let covs = Covs::new_data_from_vars(vars);
        Samples {
            phe: self.phe,
            names: self.names,
            covs: Some(covs),
            samples_n: self.samples_n,
        }
    }

    pub fn samples(&self) -> &Samples {
        self
    }

    pub fn samples_n(&self) -> usize {
        self.samples_n
    }
    pub fn phe(&self) -> &Phe {
        &self.phe
    }
    pub fn ys(&self) -> &[B8] {
        self.phe.phe_inner().inner()
        //&self.phe
    }
    /*     pub fn ys_f64(&self) -> Vec<f64> {
           let ys = self.ys();
           (0..self.samples_n())
               .map(|ni| bit::bget(ys, ni) as i32 as f64)
               .collect()
       }
    */
    pub fn names(&self) -> &[Name] {
        self.names.as_ref().unwrap()
    }
}

impl CovsTrait for Samples {
    fn covs(&self) -> Option<&Covs> {
        self.covs.as_ref()
    }
    //fn vals(&self) -> Option<&Vec<Vec<f64>>> {
    //    self.covs().unwrap().vals()
    //}
    //fn cov_indexs(&self) -> Option<&Vec<CovIndex>> {
    //    self.covs().unwrap().cov_indexs()
    //}
    //fn covs_n(&self)->usize {
    //    self.cov_indexs().unwrap().len()
    //}
}

/* /// SamplesSingle
/// one value (label) for onebool
pub struct SamplesSingle {
    // change to Labels for several phenotypes
    //ys: Y,
    ys: Vec<bool>,
    //ys: Vec<u8>,
    names: Option<Vec<Name>>,
    //names: Option<Vec<String>>,
    // move to DatasetBiNew? -> put here for a while because covs also needs names
    // also, Samples only can be input of fn
    covs: Option<Covs>,
    // fields: Option<SamplesFields> // input
    //qc: Option<SamplesQc>,
    samples_n: usize,
}

impl CovsTrait for SamplesSingle {
    fn covs(&self) -> Option<&Covs> {
        self.covs.as_ref()
    }
    //fn vals(&self) -> Option<&Vec<Vec<f64>>> {
    //    self.covs().unwrap().vals()
    //}
    //fn cov_indexs(&self) -> Option<&Vec<CovIndex>> {
    //    self.covs().unwrap().cov_indexs()
    //}
    //fn covs_n(&self)->usize {
    //    self.cov_indexs().unwrap().len()
    //}
}

impl SamplesSingle {
    pub fn new_empty() -> SamplesSingle {
        SamplesSingle {
            ys: vec![],
            names: None,
            covs: None,
            samples_n: 0,
        }
    }

    pub fn samples(&self) -> &SamplesSingle {
        self
    }
    pub fn samples_n(&self) -> usize {
        self.samples_n
    }
    //fn ys(&self) -> &[u8] {
    pub fn ys(&self) -> &[bool] {
        &self.ys
    }
} */

/*
// calc by program
pub struct SamplesQc {
    call_rate: Vec<f64>,
    // ...
}
*/

// TODO: unnecessary? since Covs only should not be used;
// usually accompanied with Samples
// should implement these in Samples?
pub trait CovsTrait {
    fn covs(&self) -> Option<&Covs>;
    //fn vals(&self) -> Option<&Vec<Vec<f64>>>;
    //fn cov_indexs(&self) -> Option<&Vec<CovIndex>>;
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

// not Vec<Cov> since one field in Covs cannot be stored
// TODO: Covs must not be None since "const" should be in
// -> or "const " should not be here?? since "const" is not covariates
#[derive(Clone)]
pub struct Covs {
    cov_indexs: Vec<CovId>,
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
    //pub fn new() -> Covs {
    //    Covs {
    //        cov_indexs: vec![],
    //        vals: vec![vec![]],
    //    }
    //}
    //pub fn new_with_names(fin_covs: &str, names: &[Name]) -> Covs {
    //    // in order of names
    //    Covs {
    //        cov_indexs: vec![],
    //        vals: vec![vec![]],
    //    }
    //}
    pub fn new_data_from_vars(vars: Vec<Var>) -> Covs {
        let cov_indexs = (&vars)
            .iter()
            .map(|x| CovId::new(x.name().to_owned(), x.kind()))
            .collect();
        let vals = vars.into_iter().map(|x| x.vals_consume()).collect();
        Covs { cov_indexs, vals }
    }

    pub fn new(
        fin_var: Option<&Path>,
        n_in: usize,
        sample_id_to_n: &HashMap<String, usize>,
        //use_samples: &[bool],
    ) -> Option<Covs> {
        if fin_var.is_none() {
            // or return const-only Covs
            return None;
        }

        let cols_name = text::load_table_header(fin_var.unwrap()).unwrap();

        let valss = text::load_table(fin_var.unwrap(), true).unwrap();

        let n = sample_id_to_n.keys().len();

        // TODO: cov could be extracted or non-extracted
        // -> fix for whole code
        // -> usually non-extracted?
        assert!(
            (valss[0].len() == n_in) | (valss[0].len() == n),
            "vals.len(), n_in, n: {}, {}, {}",
            valss[0].len(),
            n_in,
            n
        );

        let mut covs_indexs: Vec<CovId> = Vec::with_capacity(n);
        let mut covs_vals: Vec<Vec<f64>> = vec![];
        //let mut vars: Vec<Var> = Vec::with_capacity(n);

        //let p = cols_name.len() - 2;

        let col_fid = 0;
        let col_iid = 1;

        for pi in 0..cols_name.len() {
            if (pi == col_fid) || (pi == col_iid) {
                continue;
            }

            let cov_index = CovId::new_cov(cols_name[pi].clone());

            //let mut var = Var::construct_var(kind, cols_name[pi].clone());

            let vals = cov::vals_align_id(
                &valss[pi],
                (&valss[col_fid], &valss[col_iid]),
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
        Some(covs)
    }

    /// convert form row major to column major and return
    /// mainly for smartcore::DenseMatrix
    pub fn vals_column_major(&self) -> Vec<Vec<f64>> {
        let vals = self.vals().unwrap();
        let samples_n = vals[0].len();
        //let covs_n = self.covs_n();
        let mut vals_col = Vec::new();
        for _ in 0..samples_n {
            vals_col.push(Vec::new());
        }
        for v in vals.iter() {
            for (x, v_col) in v.iter().zip(vals_col.iter_mut()) {
                v_col.push(*x)
            }
        }
        println!(
            "col_major: sample {} x var {}",
            vals_col.len(),
            vals_col[0].len()
        );
        vals_col
    }

    /// convert form row major to col major vec
    /// mainly for smartcore::DenseMatrix
    pub fn vals_column_major_vec(&self) -> (Vec<f64>, usize, usize) {
        let vals = self.vals().unwrap();
        let samples_n = vals[0].len();
        let covs_n = self.covs_n();

        let mut vals_col = vec![0.0f64; samples_n * covs_n];
        for (ci, col_v) in vals.iter().enumerate() {
            for (ri, x) in col_v.iter().enumerate() {
                vals_col[ci + ri * covs_n] = *x
            }
        }
        (vals_col, samples_n, covs_n)
    }

    pub fn vals_id(&self, name: &str) -> &[f64] {
        /*
        println!("name {}", name);
        for x in self.cov_indexs().unwrap().iter() {
            println!("x {}", x.name())
        }
         */

        //println!("name {:?}", name);
        //println!("cov_indes {:?}", self.cov_indexs());

        let index = self
            .cov_indexs()
            .unwrap()
            .iter()
            .position(|x| x.name() == name)
            .unwrap();
        &self.vals().unwrap()[index]
    }
}

// mainly for Wgt
#[derive(Debug, Clone)]
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

#[inline]
//pub fn sample_id(fid: String, iid: &str) -> String {
pub fn sample_id(fid: impl Into<String>, iid: &str) -> String {
    fid.into() + ":" + iid
}

/// sample_id -> n
pub fn create_sample_id_to_n(fin: &Path, use_samples: &[bool]) -> HashMap<String, usize> {
    let fin_fam = plink::fname_fam_exist_chrom(fin).unwrap();
    //let fin_fam = plink::get_fname_fam(fin, Some(1));
    // fid, iid
    let cols = [0usize, 1];
    let samples_in: Vec<Vec<String>> = text::load_table_cols(&fin_fam, &cols, false).unwrap();

    assert_eq!(samples_in[0].len(), use_samples.len());

    /*
    let mut samples: Vec<(String, String)> = Vec::with_capacity(n_use);
    // load fid and iid
    for vi in 0..samples_in[0].len() {
        samples.push((samples_in[0][vi].clone(), samples_in[1][vi].clone()));
    }
    */

    let n = use_samples.iter().filter(|&v| *v).count();

    let mut sample_id_to_n: HashMap<String, usize> = HashMap::with_capacity(n);

    let mut ni = 0;
    for (n_in_i, v) in use_samples.iter().enumerate() {
        if *v {
            sample_id_to_n.insert(
                sample_id(samples_in[0][n_in_i].clone(), &samples_in[1][n_in_i]),
                ni,
            );
            ni += 1;
        }
    }

    assert_eq!(ni, n);

    sample_id_to_n

    /*
    for n_in in 0..n {
        sample_id_to_n.insert(sample_id(samples_in[0][ni].clone(), &samples_in[1][ni]), ni);
        //sample_in_to_index.insert(samples_in[0][si].clone() + ":" + &samples_in[1][si], si);
    }
    */
}

#[cfg(test)]
mod tests {
    //use super::*;
}
