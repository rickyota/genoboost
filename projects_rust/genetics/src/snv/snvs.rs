//! Snvs
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait

//use crate::{genot_io, snv, vec};
//use crate::{GenotFormat, SnvId};
//use std::path::Path;
use crate::{genot_io, snv, vec, Chrom, GenotFile, SnvId};

#[derive(Clone, Default)]
pub struct Snvs {
    snv_ids: Vec<SnvId>,
    mafs: Option<Vec<f64>>,
    // mafs_outside: Option<> // for missing value and pass maf file as argument
    snvs_n: usize, // not always len(snv_ids) since sometimes snv_ids are empty
                   //qc: Option<SnvsQc>, //later
}

impl Snvs {
    pub fn snvs_n(&self) -> usize {
        self.snvs_n
    }

    pub fn snvs(&self) -> &Self {
        &self
    }
    pub fn mafs(&self) -> Option<&Vec<f64>> {
        self.mafs.as_ref()
    }
    pub fn snv_ids(&self) -> &[SnvId] {
        &self.snv_ids
    }

    pub fn snv_ids_mut(&mut self) -> &mut [SnvId] {
        &mut self.snv_ids
    }

    pub fn new_plink_use_snvs(
        fin_genot: &GenotFile,
        use_snvs: Option<&[bool]>,
        m: Option<usize>, // avoid loading again
    ) -> Self {
        let m = match m {
            Some(x) => x,
            None => match use_snvs {
                Some(x) => vec::count_true(x),
                None => genot_io::compute_num_snv(&fin_genot).unwrap(),
            },
        };

        //let m_in: usize = genot_io::compute_num_snv(fin, gfmt).unwrap();

        let snvs_in = genot_io::load_snvs(fin_genot);
        let snv_ids = match use_snvs {
            Some(x) => snv::extract_snvs_consume(snvs_in, x, m),
            None => snvs_in,
        };
        //let use_snvs = vec![true; m_in];
        //let snv_indexs = snv::extract_snvs_consume(snvs_in, &use_snvs, m_in);
        let snvs = Self::new_from_snv_ids(snv_ids);
        snvs
    }

    pub fn new_from_snv_ids(snv_ids: Vec<SnvId>) -> Self {
        Self::new(snv_ids, None)
        //let snvs_n = snvs.len();
        //// what for?
        //let snv_indexs = snvs.into_iter().collect();
        //Self {
        //    snv_indexs,
        //    mafs: None,
        //    snvs_n,
        //}
    }

    // tmp
    pub fn new_empty() -> Self {
        Self::new(vec![], None)
        //Self {
        //    snv_indexs: vec![],
        //    mafs: None,
        //    snvs_n: 0,
        //}
    }

    pub fn new(snv_ids: Vec<SnvId>, mafs: Option<Vec<f64>>) -> Self {
        let snvs_n = snv_ids.len();
        Self {
            snv_ids,
            mafs,
            snvs_n,
        }
    }

    // or set_maf(&self, mafs)->Self{}
    // -> this way makes difficult to set_maf if in Dataset{Snvs}
    pub fn set_maf(&mut self, mafs: Vec<f64>) {
        self.mafs = Some(mafs);
    }

    pub fn extract_snvs(self, use_snvs: &[bool]) -> Self {
        let Self {
            snv_ids,
            mafs,
            snvs_n,
        } = self;

        assert_eq!(snvs_n, use_snvs.len());

        let snv_ids_use: Vec<SnvId> = snv_ids
            .iter()
            .zip(use_snvs.iter())
            .filter(|(_, b)| **b)
            .map(|(x, _)| x.clone())
            .collect();

        let mafs_use: Option<Vec<f64>> = mafs.map(|x| {
            x.iter()
                .zip(use_snvs.iter())
                .filter(|(_, b)| **b)
                .map(|(x, _)| x.clone())
                .collect()
        });

        let n = snv_ids_use.len();

        Self {
            snv_ids: snv_ids_use,
            mafs: mafs_use,
            snvs_n: n,
        }
    }

    pub fn use_snvs_chrom(&self, chrom: &Chrom) -> Vec<bool> {
        self.snv_ids().iter().map(|x| x.chrom() == chrom).collect()
    }

    //pub fn extract_chrom_indexs(&self, chrom: &Chrom) -> Vec<usize> {
    //    self.snv_ids()
    //        .iter()
    //        .enumerate()
    //        .filter(|(_, x)| x.chrom() == chrom)
    //        .map(|(i, _)| i)
    //        .collect()
    //}

    pub fn positions(&self) -> Vec<usize> {
        self.snv_ids().iter().map(|x| x.pos()).collect()
    }

    // use extrac_chrom_indexs() and positions()
    //pub fn positions(&self, chrom: &Chrom) -> Vec<usize> {
    //    self.snv_ids()
    //        .iter()
    //        .filter(|x| x.chrom() == chrom)
    //        .map(|x| x.pos())
    //        .collect()
    //}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_snvs() {
        let snv_ids = vec![SnvId::default(); 3];
        let snvs = Snvs::new_from_snv_ids(snv_ids);
        let use_snvs = [true, false, true];
        let snvs_use = snvs.extract_snvs(&use_snvs);

        //let snv_use_exp = Snvs::new_from_snv_index(vec![SnvId::default(); 2]);

        assert_eq!(snvs_use.snv_ids().len(), 2);
    }
}
