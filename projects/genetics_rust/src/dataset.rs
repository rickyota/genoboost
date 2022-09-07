//! Dataset
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait
//! For now, do not use generics, use in dataset_future.rs

use crate::{plink, samples, snv, wgt::WgtTrait, Covs, Genot, Samples, Snvs};
use std::{path::Path, time::Instant};

// 8 x 1 bit
//type B8 = u8;

#[derive(Clone)]
pub struct Dataset {
    genot: Genot,
    snvs: Snvs,
    samples: Samples,
}

impl Dataset {
    pub fn new(
        fin: &Path,
        fin_snv: Option<&Path>,
        fin_sample: Option<&Path>,
        fin_cov: Option<&Path>,
        use_missing: bool,
    ) -> Dataset {
        let start = Instant::now();

        plink::check_valid_fin(fin);
        let m_in: usize = plink::compute_num_snv(fin).unwrap();
        println!("m_in: {}", m_in);
        let n_in: usize = plink::compute_num_sample(fin).unwrap();
        println!("n_in: {}", n_in);
        // load snvs
        let snvs_in = plink::load_snvs(fin, m_in);
        let (m, use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);
        if m == 0 {
            panic!("Using SNVs are zero. Please check fin_snv.")
        }
        println!("m: {}", m);
        let snv_indexs = snv::extract_snvs_consume(snvs_in, &use_snvs, m);
        let (n, use_samples) = plink::make_use_samples(fin_sample, fin, n_in);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }
        println!("n: {}", n);

        let genot =
            plink::load::generate_genot(fin, m, n, &use_snvs, Some(&use_samples), use_missing);

        let ys: Vec<bool> = plink::load_ys(fin, n, &use_samples);
        println!("ys {}", ys[0]);
        let sample_id_to_n = samples::create_sample_id_to_n(fin, &use_samples);
        let covs: Option<Covs> = Covs::new(fin_cov, n_in, &sample_id_to_n);

        let snvs = Snvs::new_data_from_snv_index(snv_indexs);
        let samples = Samples::new(ys, None, covs, n);

        println!(
            "It took {} seconds to create Dataset.",
            start.elapsed().as_secs()
        );

        let dataset = Dataset {
            genot,
            snvs,
            samples,
        };
        dataset.check();
        dataset
    }

    /// for prune SNVs by loss
    pub fn extract_snvs(self, use_snvs: &[bool]) -> Dataset {
        let Dataset {
            genot,
            snvs,
            samples,
        } = self;

        let genot_use = genot.extract_snvs(use_snvs);
        let snvs_use = snvs.extract_snvs(use_snvs);

        let dataset = Dataset {
            genot: genot_use,
            snvs: snvs_use,
            samples,
        };
        dataset.check();
        dataset
    }

    pub fn genot(&self) -> &Genot {
        &self.genot
    }
    pub fn snvs(&self) -> &Snvs {
        &self.snvs
    }
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    // in SnvTrait
    //fn snvs(&self) -> &Snvs {
    //    &self.snvs
    //}
    //fn samples(&self) -> &Samples {
    //    &self.samples
    //}
    // move to impl DatasetBiNew{}
    fn check(&self) {
        // check if samples_n, snvs_n are the same
    }

    pub fn new_score<W: WgtTrait>(
        fin: &Path,
        fin_sample: Option<&Path>,
        fin_cov: Option<&Path>,
        wgts: &mut [W],
        //wgts: &[WgtBoost],
    ) -> Dataset {
        let m_in: usize = plink::compute_num_snv(fin).unwrap();
        println!("{}", m_in);
        let n_in: usize = plink::compute_num_sample(fin).unwrap();
        println!("{}", n_in);

        let (n, use_samples) = plink::make_use_samples(fin_sample, fin, n_in);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        // TODO: make (string, string) as key of hashmap
        //let samples_id = plink::load_samples_id(fin, &use_samples);

        let sample_id_to_n = samples::create_sample_id_to_n(fin, &use_samples);
        let ys: Vec<bool> = plink::load_ys(fin, n, &use_samples);
        let covs: Option<Covs> = Covs::new(fin_cov, n_in, &sample_id_to_n);
        //// covs could not exist
        //let covs: Option<Vec<Var>> =
        //    cov::load_vars(fin_cov, n_in, &sample_id_to_n, cov::CovKind::Cov);

        let samples = Samples::new(ys, None, covs, n);

        //let mut wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(fin_wgt);
        //wgt_boost::io::set_covs(&mut wgts, covs.as_deref(), n);

        // set genotype index in wgt
        let use_missing = true;
        let genot = plink::load_score::load_genotypes_for_score(
            fin,
            wgts,
            n,
            Some(&use_samples),
            use_missing,
        );

        Dataset {
            genot,
            // unncessesary?
            snvs: Snvs::new_empty(),
            samples,
        }
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
