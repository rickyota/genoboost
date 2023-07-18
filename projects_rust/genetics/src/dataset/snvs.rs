//! Dataset
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait

use crate::{SnvId, GenotFormat};
use crate::{io_genot, snv};
use std::path::Path;

#[derive(Clone)]
pub struct Snvs {
    snv_indexs: Vec<SnvId>,
    mafs: Option<Vec<f64>>, // for ss: indicated by user or for mean value of genotype
    snvs_n: usize,
    //qc: Option<SnvsQc>, //later
}

impl Snvs {
    //pub fn new(fin: &str) -> Snvs {
    //    Snvs {
    //        snv_indexs: vec![],
    //        mafs: None,
    //        snvs_n: 1,
    //    }
    //}
    /*
    pub fn new_data_from_snv(snvs: Vec<Snv>) -> Snvs {
        let snvs_n = snvs.len();
        let snv_indexs = snvs.into_iter().map(|x| x.snv_index_consume()).collect();
        Snvs {
            snv_indexs,
            mafs: None,
            snvs_n,
        }
    }
     */
    pub fn new_data_from_snv_index(snvs: Vec<SnvId>) -> Snvs {
        let snvs_n = snvs.len();
        let snv_indexs = snvs.into_iter().collect();
        Snvs {
            snv_indexs,
            mafs: None,
            snvs_n,
        }
    }

    pub fn new_plink(fin: &Path, gfmt: GenotFormat) -> Snvs {
        let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        let snvs_in = io_genot::load_snvs(fin, gfmt);
        let use_snvs = vec![true; m_in];
        let snv_indexs = snv::extract_snvs_consume(snvs_in, &use_snvs, m_in);
        let snvs = Snvs::new_data_from_snv_index(snv_indexs);
        snvs
    }

    // tmp
    pub fn new_empty() -> Snvs {
        Snvs {
            snv_indexs: vec![],
            mafs: None,
            snvs_n: 0,
        }
    }

    // or set_maf(&self, mafs)->Self{}
    // -> this way makes difficult to set_maf if in Dataset{Snvs}
    pub fn set_maf(&mut self, mafs: Vec<f64>) {
        self.mafs = Some(mafs);
    }

    pub fn extract_snvs(self, use_snvs: &[bool]) -> Snvs {
        let Snvs {
            snv_indexs,
            mafs,
            snvs_n,
        } = self;

        assert_eq!(snvs_n, use_snvs.len());

        let snv_indexs_use: Vec<SnvId> = snv_indexs
            .iter()
            .zip(use_snvs.iter())
            .filter(|(_, b)| **b)
            .map(|(x, _)| x.clone())
            .collect();

        let mafs_use: Option<Vec<f64>> = if let Some(x) = mafs {
            Some(
                x.iter()
                    .zip(use_snvs.iter())
                    .filter(|(_, b)| **b)
                    .map(|(x, _)| x.clone())
                    .collect(),
            )
        } else {
            None
        };

        let n = snv_indexs_use.len();

        Snvs {
            snv_indexs: snv_indexs_use,
            mafs: mafs_use,
            snvs_n: n,
        }
    }

    pub fn snvs_n(&self) -> usize {
        self.snvs_n
    }

    pub fn snvs(&self) -> &Snvs {
        &self
    }
    pub fn mafs(&self) -> Option<&Vec<f64>> {
        self.mafs.as_ref()
    }
    pub fn snv_indexs(&self) -> &[SnvId] {
        &self.snv_indexs
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
    //use crate::Snv;
    //use genetics_rust::Snv;
}
