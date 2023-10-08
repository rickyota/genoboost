//!
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait
//!
//! TODO: assume FID==IID for plink1 or plink2
//!

pub mod base_phe;
pub mod covs;
pub mod iterator;
pub mod phe;
pub mod prelude;

use std::collections::HashMap;
use std::path::Path;

use crate::{io_genot, GenotFormat};
use crate::{Var, B8};
use cmatrix::BaseCMatrix;
pub use covs::{CovId, Covs, CovsTrait};
use prelude::*;

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
    pub fn new(
        ys: Option<&[bool]>,
        names: Option<Vec<Name>>,
        covs: Option<Covs>,
        samples_n: usize,
    ) -> Samples {
        let phe = ys.map(|x| Phe::new(x));
        Samples {
            phe,
            names,
            covs,
            samples_n,
        }
    }
    pub fn new_data(ys: &[bool], names: Option<Vec<Name>>, samples_n: usize) -> Samples {
        let phe = Phe::new(ys);
        Samples {
            phe: Some(phe),
            names,
            covs: None,
            samples_n,
        }
    }

    /// for score
    pub fn new_nophe(names: Option<Vec<Name>>, covs: Option<Covs>, samples_n: usize) -> Samples {
        //let phe = Phe::new_empty(samples_n);
        Samples {
            //phe,
            phe: None,
            names,
            covs,
            samples_n,
        }
    }

    pub fn new_empty() -> Samples {
        Samples {
            phe: None,
            //phe: Phe::new_empty(0),
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

    // TMP: merge to phe_unwrap();
    //pub fn phe(&self) -> &Phe {
    //    self.phe_unwrap()
    //    //self.phe.as_ref().unwrap()
    //    //&self.phe
    //}

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

/* #[inline]
//pub fn sample_id(fid: String, iid: &str) -> String {
pub fn sample_id(fid: impl Into<String>, iid: &str) -> String {
    fid.into() + ":" + iid
} */

/// sample_id -> n
pub fn create_sample_id_to_n(
    fin: &Path,
    gfmt: GenotFormat,
    use_samples: Option<&[bool]>,
) -> HashMap<String, usize> {
    //let fin_fam = io_genot::fname_fam_exist_chrom(fin, gfmt).unwrap();
    //let fin_fam = plink::get_fname_fam(fin, Some(1));
    // fid, iid
    //let cols = [0usize, 1];
    //let samples_in: Vec<Vec<String>> = textfile::load_table_cols(&fin_fam, &cols, false).unwrap();

    /*     let cols = [0usize];
    let samples_in: Vec<Vec<String>> = textfile::load_table_cols(&fin_fam, &cols, false).unwrap(); */

    let samples_in: Vec<String> = io_genot::load_samples_id(fin, gfmt, use_samples);

    /*
    let mut samples: Vec<(String, String)> = Vec::with_capacity(n_use);
    // load fid and iid
    for vi in 0..samples_in[0].len() {
        samples.push((samples_in[0][vi].clone(), samples_in[1][vi].clone()));
    }
    */

    //let n = use_samples.iter().filter(|&v| *v).count();
    //let mut sample_id_to_n: HashMap<String, usize> = HashMap::with_capacity(n);
    let mut sample_id_to_n: HashMap<String, usize> = HashMap::new();

    for n_in_i in 0..samples_in.len() {
        if sample_id_to_n.contains_key(&samples_in[n_in_i]) {
            panic!("Sample id is duplicated: {}", samples_in[n_in_i])
        }
        sample_id_to_n.insert(samples_in[n_in_i].clone(), n_in_i);
    }

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
