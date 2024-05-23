//!
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait
//!
//! TODO: assume FID==IID for plink1 and plink2
//!

pub mod base_phe;
pub mod iterator;
pub mod phe;
pub mod prelude;

use std::collections::HashMap;

use crate::{Covs, CovsTrait};
//use crate::{genot_io, vec, GenotFormat, Var, B8};
use crate::{genot_io, vec, GenotFile, B8};
//use crate::{genot_io, vec, GenotFile, Cov, B8};
use cmatrix::prelude::*;
//use cmatrix::BaseCMatrix;
pub use prelude::*;

// maybe change to (String, String)
//type Name = (String, String);
type Name = String;

// SampleTrait: CovsTrait; using subtrait is troublesome
// indicate both

#[derive(Clone)]
pub struct Samples {
    // change to Labels for several phenotypes
    phe: Option<Phe>,
    names: Option<Vec<Name>>,
    // move to Dataset? -> put here for a while because covs also needs names
    // also, Samples only can be input of fn
    // also, you don't need to modify Dataset when adding another field
    covs: Option<Covs>,
    // fields: Option<SamplesFields> // input
    //qc: Option<SamplesQc>,
    samples_n: usize,
}

impl Samples {
    pub fn new_plink(
        fin_genot: &GenotFile,
        phe_buf: Option<&[u8]>,
        phe_name: Option<&str>,
        cov_name: Option<&str>,
        use_samples: Option<&[bool]>,
        n: Option<usize>,
    ) -> Self {
        let n = match n {
            Some(x) => x,
            None => match use_samples {
                Some(x) => vec::count_true(x),
                None => genot_io::compute_num_sample(fin_genot).unwrap(),
            },
        };

        let samples_id = genot_io::load_samples_id(fin_genot, use_samples);

        let sample_id_to_n = create_sample_id_to_n_samples(&samples_id);
        //let sample_id_to_n = create_sample_id_to_n(fin_genot, use_samples);

        let ys: Option<Vec<bool>> =
            genot_io::load_ys_buf(fin_genot, phe_buf, phe_name, &sample_id_to_n);

        let covs: Option<Covs> =
            cov_name.map(|x| Covs::new_buf(phe_buf, fin_genot, x, &sample_id_to_n));

        //let samples_id = genot_io::load_samples_id(fin_genot, use_samples);

        let phe = ys.map(|x| Phe::new(&x));
        Self::new(phe, Some(samples_id), covs, n)

        //Samples::new_bool(ys.as_deref(), Some(samples_id), covs, n)
    }

    //pub fn new_bool(
    //    ys: Option<&[bool]>,
    //    names: Option<Vec<Name>>,
    //    covs: Option<Covs>,
    //    samples_n: usize,
    //) -> Self {
    //    let phe = ys.map(|x| Phe::new(x));
    //    Self::new(phe, names, covs, samples_n)
    //    //Self {
    //    //    phe,
    //    //    names,
    //    //    covs,
    //    //    samples_n,
    //    //}
    //}

    pub fn new(
        phe: Option<Phe>,
        names: Option<Vec<Name>>,
        covs: Option<Covs>,
        samples_n: usize,
    ) -> Self {
        Self {
            phe,
            names,
            covs,
            samples_n,
        }
    }

    //pub fn new_data(ys: &[bool], names: Option<Vec<Name>>, samples_n: usize) -> Self {
    //    let phe = Phe::new(ys);
    //    Self {
    //        phe: Some(phe),
    //        names,
    //        covs: None,
    //        samples_n,
    //    }
    //}

    /// for score
    //pub fn new_nophe(names: Option<Vec<Name>>, covs: Option<Covs>, samples_n: usize) -> Self {
    //    Self::new(None, names, covs, samples_n)
    //    //let phe = Phe::new_empty(samples_n);
    //    //Self {
    //    //    //phe,
    //    //    phe: None,
    //    //    names,
    //    //    covs,
    //    //    samples_n,
    //    //}
    //}

    //pub fn new_empty() -> Self {
    //    Self::new(None, None, None, 0)
    //}

    //pub fn add_covs_with_vars(self, vars: Vec<Cov>) -> Self {
    //    let covs = Covs::new_data_from_vars(vars);
    //    Self {
    //        phe: self.phe,
    //        names: self.names,
    //        covs: Some(covs),
    //        samples_n: self.samples_n,
    //    }
    //}

    pub fn samples(&self) -> &Self {
        self
    }

    pub fn samples_n(&self) -> usize {
        self.samples_n
    }

    pub fn phe(&self) -> Option<&Phe> {
        self.phe.as_ref()
        //self.phe.as_ref().unwrap()
        //&self.phe
    }

    pub fn phe_unwrap(&self) -> &Phe {
        self.phe().unwrap()
        //self.phe.as_ref().unwrap()
        //&self.phe
    }

    pub fn ys_unwrap(&self) -> &[B8] {
        self.phe_unwrap().phe_inner().inner()
        //self.phe_unwrap().phe_inner().inner()
        //self.phe.phe_inner().inner()
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

/*
// calc by program
pub struct SamplesQc {
    call_rate: Vec<f64>,
    // ...
}
*/

/// sample_id -> n_i (not n_in_i)
pub fn create_sample_id_to_n(
    fin_genot: &GenotFile,
    use_samples: Option<&[bool]>,
) -> HashMap<String, usize> {
    let samples_id: Vec<String> = genot_io::load_samples_id(fin_genot, use_samples);

    create_sample_id_to_n_samples(&samples_id)

    //let n = samples_in.len();
    //let mut sample_id_to_n: HashMap<String, usize> = HashMap::with_capacity(n);
    ////let mut sample_id_to_n: HashMap<String, usize> = HashMap::new();

    //for n_i in 0..samples_in.len() {
    //    if sample_id_to_n.contains_key(&samples_in[n_i]) {
    //        panic!("Sample id is duplicated: {}", samples_in[n_i])
    //    }
    //    sample_id_to_n.insert(samples_in[n_i].clone(), n_i);
    //}

    //sample_id_to_n
}

pub fn create_sample_id_to_n_samples(samples_id: &[String]) -> HashMap<String, usize> {
    let n = samples_id.len();
    let mut sample_id_to_n: HashMap<String, usize> = HashMap::with_capacity(n);
    //let mut sample_id_to_n: HashMap<String, usize> = HashMap::new();

    for n_i in 0..samples_id.len() {
        if sample_id_to_n.contains_key(&samples_id[n_i]) {
            panic!("Sample id is duplicated: {}", samples_id[n_i])
        }
        sample_id_to_n.insert(samples_id[n_i].clone(), n_i);
    }

    // this is for n_in_i
    /*     if use_samples.is_some() {
        let mut ni = 0;
        for (n_in_i, v) in use_samples.unwrap().iter().enumerate() {
            if *v {
                sample_id_to_n.insert(
                    samples_in[n_in_i].clone(),
                    //sample_id(samples_in[0][n_in_i].clone(), &samples_in[1][n_in_i]),
                    ni,
                );
                ni += 1;
            }
        }
        assert_eq!(ni, use_samples.unwrap().iter().filter(|&v| *v).count());
    } else {
        for n_in_i in 0..samples_in.len() {
            sample_id_to_n.insert(samples_in[n_in_i].clone(), n_in_i);
        }
    }; */

    sample_id_to_n
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_sample_id_to_n_samples() {
        let samples_in = vec!["aaa".to_string(), "bbb".to_string()];

        let m = create_sample_id_to_n_samples(&samples_in);

        let m_ans = HashMap::from([("aaa".to_string(), 0), ("bbb".to_string(), 1)]);
        //let mut m_ans = HashMap::new();
        //m_ans.insert("aaa".to_string(), 0);
        //m_ans.insert("bbb".to_string(), 1);

        assert_eq!(m, m_ans);
    }
}
