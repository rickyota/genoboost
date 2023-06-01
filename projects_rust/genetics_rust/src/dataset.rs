//! Dataset
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait
//! For now, do not use generics, use in dataset_future.rs

use rayon::iter::{ParallelBridge, ParallelIterator};

use crate::{plink, samples, snv, vec, wgt::WgtTrait, Covs, Genot, Samples, Snvs, Wgts, genot::{BaseGenot, BaseGenotSnv, BaseGenotMut}};
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
    /// filt_snv: for loading part of SNVs mainly for pruning snv
    pub fn new(
        fin: &Path,
        fin_phe: Option<&Path>,
        phe_name: Option<&str>,
        fin_snv: Option<&Path>,
        fin_sample: Option<&Path>,
        fin_cov: Option<&Path>,
        use_missing: bool,
        // TODO: merge filt_snv and fin_snv into use_snvs,
        filt_snv: Option<&[bool]>,
        snvs_train: Option<&Snvs> //for use_missing of vali
    ) -> Self {
        let start = Instant::now();

        plink::check_valid_fin(fin);
        let m_in: usize = plink::compute_num_snv(fin).unwrap();
        log::debug!("m_in: {}", m_in);
        let n_in: usize = plink::compute_num_sample(fin).unwrap();
        log::debug!("n_in: {}", n_in);
        // load snvs
        let snvs_in = plink::load_snvs(fin, m_in);
        let (mut m, mut use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);
        //let (m, use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);

        if let Some(filt_snv_) = filt_snv {
            log::debug!("filt_snv before m: {}", m);
            assert_eq!(filt_snv_.len(), m);
            use_snvs = vec::and_in_bool_vec(&use_snvs, filt_snv_);
            m = vec::count_true(&use_snvs);
            log::debug!("filt_snv m: {}", m);
        }

        if m == 0 {
            panic!("Using SNVs are zero. Please check fin_snv.")
        }
        log::debug!("m: {}", m);
        let snv_indexs = snv::extract_snvs_consume(snvs_in, &use_snvs, m);
        let (n, use_samples) = plink::make_use_samples(fin_sample, fin, n_in);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }
        log::debug!("n: {}", n);

        // do not fill missing here; below
        let genot =
            plink::load::generate_genot(fin, m, n, &use_snvs, Some(&use_samples), true);
            //plink::load::generate_genot(fin, m, n, &use_snvs, Some(&use_samples), use_missing);

        let ys: Vec<bool> = plink::load_ys(fin, fin_phe, phe_name, n, &use_samples);
        log::debug!("ys {}", ys[0]);
        let sample_id_to_n = samples::create_sample_id_to_n(fin, &use_samples);
        let covs: Option<Covs> = Covs::new(fin_cov, &sample_id_to_n);

        let snvs = Snvs::new_data_from_snv_index(snv_indexs);
        let samples = Samples::new(&ys, None, covs, n);

        log::debug!(
            "It took {} seconds to create Dataset.",
            start.elapsed().as_secs()
        );

        let mut dataset = Dataset {
            genot,
            snvs,
            samples,
        };
        dataset.check();


        dataset.compute_maf();


        // FIXME: flip A1, A2 for lims2
        // TOOD: when s2 is not minor homo but major homo
        // if force_a1_minor{}


        // fill missing here
        // for training; use maf assuming hwe
        // for validation; use the training maf
        if !use_missing{
            // for training
            if let Some(snvs_train)=snvs_train{
                // for val
                // unwrap to raise error when None
                let mafs=snvs_train.mafs().unwrap();
                dataset.fill_missing_maf(Some(mafs));
            } else{
                // for training 
                dataset.fill_missing_maf(None);
            }
        }

        dataset
    }

    /// for prune SNVs by loss
    pub fn extract_snvs(self, use_snvs: &[bool]) -> Self {
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
    pub fn genot_mut(&mut self) -> &mut Genot {
        &mut self.genot
    }
    pub fn snvs(&self) -> &Snvs {
        &self.snvs
    }
    pub fn snvs_mut(&mut self) -> &mut Snvs {
        &mut self.snvs
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

    pub fn fill_missing_maf(&mut self, mafs: Option<&Vec<f64>>){
        // TODO: how to avoid clone?
        let mafs=if let Some(mafs_)=mafs{
            // for validation
            mafs_.clone()
        }else{
            // for training
            self.snvs().mafs().unwrap().clone()
        };
        let genot =self.genot_mut();
        genot.iter_snv_mut()
            .zip(mafs.iter())
            .par_bridge()
            .for_each(|(mut g_snv, maf)| plink::load::fill_missing_maf(&mut g_snv, *maf));
        //genot.iter_snv_mut()
        //    .par_bridge()
        //    .for_each(|mut g_snv| plink::load::fill_missing_maf(&mut g_snv));
    }

    pub fn compute_maf(&mut self){
        let genot=self.genot();
        let m=genot.m();
        let mut mafs = vec![f64::NAN; m];
        mafs.iter_mut().zip(
        genot.iter_snv())
            .par_bridge()
            .for_each(|(maf,g_snv)| *maf=g_snv.maf());

        self.snvs_mut().set_maf(mafs)

    }

    pub fn new_score<W: WgtTrait>(
        fin: &Path,
        fin_phe: Option<&Path>,
        phe_name: Option<&str>,
        fin_sample: Option<&Path>,
        fin_cov: Option<&Path>,
        wgts: &mut [W],
        //wgts: &[WgtBoost],
        use_missing: bool, // use WgtBoosts and use wgts.use_missing()
    ) -> Self {
        let m_in: usize = plink::compute_num_snv(fin).unwrap();
        log::debug!("{}", m_in);
        let n_in: usize = plink::compute_num_sample(fin).unwrap();
        log::debug!("{}", n_in);

        let (n, use_samples) = plink::make_use_samples(fin_sample, fin, n_in);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        // TODO: make (string, string) as key of hashmap
        //let samples_id = plink::load_samples_id(fin, &use_samples);

        let sample_id_to_n = samples::create_sample_id_to_n(fin, &use_samples);
        let ys: Vec<bool> = plink::load_ys(fin, fin_phe, phe_name, n, &use_samples);
        let covs: Option<Covs> = Covs::new(fin_cov, &sample_id_to_n);
        //// covs could not exist
        //let covs: Option<Vec<Var>> =
        //    cov::load_vars(fin_cov, n_in, &sample_id_to_n, cov::CovKind::Cov);

        let samples = Samples::new(&ys, None, covs, n);

        //let mut wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(fin_wgt);
        //wgt_boost::io::set_covs(&mut wgts, covs.as_deref(), n);

        // set genotype index in wgt
        //let use_missing = true;
        let genot = plink::load_score::load_genotypes_for_score(
            fin,
            wgts,
            n,
            Some(&use_samples),
            use_missing,
        );

        Dataset {
            genot,
            // unnecessary since index is in WgtKInd
            snvs: Snvs::new_empty(),
            samples,
        }
    }

    // for genetics::score()
    pub fn new_score_genetics(
        fin: &Path,
        fin_phe: Option<&Path>,
        phe_name: Option<&str>,
        fin_sample: Option<&Path>,
        fin_cov: Option<&Path>,
        // TODO: &mut[&mut[Wgt]],
        wgts: &mut [Wgts],
        //wgts: &[WgtBoost],
    ) -> Dataset {
        let m_in: usize = plink::compute_num_snv(fin).unwrap();
        log::debug!("{}", m_in);
        let n_in: usize = plink::compute_num_sample(fin).unwrap();
        log::debug!("{}", n_in);

        let (n, use_samples) = plink::make_use_samples(fin_sample, fin, n_in);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        // TODO: make (string, string) as key of hashmap
        //let samples_id = plink::load_samples_id(fin, &use_samples);

        let sample_id_to_n = samples::create_sample_id_to_n(fin, &use_samples);

        let ys: Vec<bool> = plink::load_ys(fin, fin_phe, phe_name, n, &use_samples);
        let covs: Option<Covs> = Covs::new(fin_cov, &sample_id_to_n);
        //// covs could not exist
        //let covs: Option<Vec<Var>> =
        //    cov::load_vars(fin_cov, n_in, &sample_id_to_n, cov::CovKind::Cov);

        let samples = Samples::new(&ys, None, covs, n);

        //let mut wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(fin_wgt);
        //wgt_boost::io::set_covs(&mut wgts, covs.as_deref(), n);

        // set genotype index in wgt
        // TODO: argparse
        let use_missing = false;
        let genot = plink::load_score::load_genotypes_for_score_multiwgts(
            fin,
            wgts,
            //&mut wgts,
            n,
            Some(&use_samples),
            use_missing,
        );

        Dataset {
            genot,
            // unnecessary since index is in WgtKInd
            snvs: Snvs::new_empty(),
            samples,
        }
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
