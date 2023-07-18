//! Dataset
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait
//! For now, do not use generics, use in dataset_future.rs
//!

pub mod io_genot;
pub mod samples;
pub mod snvs;

use crate::{
    genot::{BaseGenot, BaseGenotMut, BaseGenotSnv},
    sample,
    snv, vec,
    wgt::WgtTrait,
    Covs, Genot, GenotFormat, Samples, Snvs, Wgts,
};
use rayon::iter::{ParallelBridge, ParallelIterator};
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
    // Do not use R since buf should be used twice and File::open() cannot be copied...
    //pub fn new<R: std::io::Read + Copy>(
    ///
    /// Input is one of following
    /// 1. plink2 + fin_phe
    /// 2. plink2 (cov and phe in .psam) : phe_buf is option
    /// 3. plink1 (phe) + fin_phe (cov) : phe_name is option
    /// 4. plink1 + fin_phe (cov and phe)
    pub fn new(
        fin: &Path,
        gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        //phe_buf: &[u8],
        //fin_phe: Option<&Path>,
        phe_name: Option<&str>,
        cov_name: &str,
        extract_snv_buf: Option<&[u8]>,
        //fin_snv: Option<&Path>,
        extract_sample_buf: Option<&[u8]>,
        //fin_sample: Option<&Path>,
        //fin_cov: Option<&Path>, //to deprecate
        use_missing: bool,
        // TODO: merge filt_snv and fin_snv into use_snvs,
        filt_snv: Option<&[bool]>,
        snvs_train: Option<&Snvs>, //for use_missing of vali
    ) -> Self {
        let start = Instant::now();

        io_genot::check_valid_fin(fin, gfmt);
        let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        log::debug!("m_in: {}", m_in);
        let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        log::debug!("n_in: {}", n_in);
        // load snvs
        let snvs_in = io_genot::load_snvs(fin, gfmt);
        //let snvs_in = plink::load_snvs(fin, m_in);
        let (mut m, mut use_snvs) = snv::make_use_snvs_buf(extract_snv_buf, &snvs_in);
        //let (mut m, mut use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
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
        //let snv_indexs = snv::extract_snvs_consume(snvs_in, &use_snvs, m);

        let (n, use_samples) = sample::make_use_samples_buf(extract_sample_buf, fin, gfmt);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }
        log::debug!("n: {}", n);

        // do not fill missing here; below
        let genot =
            io_genot::load::generate_genot(fin, gfmt, m, n, &use_snvs, Some(&use_samples), true);
        //plink::load::generate_genot(fin, m, n, &use_snvs, Some(&use_samples), use_missing);

        let sample_id_to_n = samples::create_sample_id_to_n(fin, gfmt, Some(&use_samples));

        let ys: Vec<bool> = io_genot::load_ys_buf(fin, gfmt, phe_buf, phe_name, &sample_id_to_n);
        //let ys: Vec<bool> = io_genot::load_ys_buf(fin, gfmt, phe_buf, phe_name, &use_samples);
        log::debug!("ys {}", ys[0]);
        let covs: Covs = Covs::new(phe_buf, fin, gfmt, cov_name, &sample_id_to_n);
        //println!("covs: {:?}", covs.clone().unwrap().cov_indexs());

        let snv_indexs = snv::extract_snvs_consume(snvs_in, &use_snvs, m);
        let snvs = Snvs::new_data_from_snv_index(snv_indexs);
        // TODO: add names here, instead of using iid in score
        let samples = Samples::new(&ys, None, Some(covs), n);

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
        if !use_missing {
            // for training
            if let Some(snvs_train) = snvs_train {
                // for val
                // unwrap to raise error when None
                let mafs = snvs_train.mafs().unwrap();
                dataset.fill_missing_maf(Some(mafs));
            } else {
                // for training
                dataset.fill_missing_maf(None);
            }
        }

        dataset
    }

    /*     /// AVOID using this func since fin_snv might be called several time, and this fials when input from /std/in
    /// use new()
    ///
    /// filt_snv: for loading part of SNVs mainly for pruning snv
    pub fn new_old(
        fin: &Path,
        gfmt: GenotFormat,
        fin_phe: Option<&Path>,
        phe_name: Option<&str>,
        cov_name: Option<&str>,
        fin_snv: Option<&Path>,
        fin_sample: Option<&Path>,
        fin_cov: Option<&Path>,
        use_missing: bool,
        // TODO: merge filt_snv and fin_snv into use_snvs,
        filt_snv: Option<&[bool]>,
        snvs_train: Option<&Snvs>, //for use_missing of vali
    ) -> Self {
        let start = Instant::now();

        io_genot::check_valid_fin(fin, gfmt);
        let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        log::debug!("m_in: {}", m_in);
        let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        log::debug!("n_in: {}", n_in);
        // load snvs
        let snvs_in = io_genot::load_snvs(fin, gfmt);
        //let snvs_in = plink::load_snvs(fin, m_in);
        let (mut m, mut use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
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
        let (n, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }
        log::debug!("n: {}", n);

        // do not fill missing here; below
        let genot =
            io_genot::load::generate_genot(fin, gfmt, m, n, &use_snvs, Some(&use_samples), true);
        //plink::load::generate_genot(fin, m, n, &use_snvs, Some(&use_samples), use_missing);

        let ys: Vec<bool> = io_genot::load_ys(fin, gfmt, fin_phe, phe_name, &use_samples);
        log::debug!("ys {}", ys[0]);
        let sample_id_to_n = samplesm::create_sample_id_to_n(fin, gfmt, Some(&use_samples));
        let covs: Option<Covs> = if fin_cov.is_none() {
            Covs::new_old(fin_cov, cov_name, &sample_id_to_n)
        } else if fin_phe.is_none() {
            Covs::new_old(fin_phe, cov_name, &sample_id_to_n)
        } else {
            None
        };

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
        if !use_missing {
            // for training
            if let Some(snvs_train) = snvs_train {
                // for val
                // unwrap to raise error when None
                let mafs = snvs_train.mafs().unwrap();
                dataset.fill_missing_maf(Some(mafs));
            } else {
                // for training
                dataset.fill_missing_maf(None);
            }
        }

        dataset
    } */

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

    pub fn fill_missing_maf(&mut self, mafs: Option<&Vec<f64>>) {
        // TODO: how to avoid clone?
        // create let mafs:Vec<f64>;
        // mafs_ref:&[f64]=if ..
        let mafs = if let Some(mafs_) = mafs {
            // for validation
            mafs_.clone()
        } else {
            // for training
            self.snvs().mafs().unwrap().clone()
        };
        let genot = self.genot_mut();
        genot
            .iter_snv_mut()
            .zip(mafs.iter())
            .par_bridge()
            .for_each(|(mut g_snv, maf)| io_genot::load::fill_missing_maf(&mut g_snv, *maf));
        //genot.iter_snv_mut()
        //    .par_bridge()
        //    .for_each(|mut g_snv| plink::load::fill_missing_maf(&mut g_snv));
    }

    pub fn compute_maf(&mut self) {
        let genot = self.genot();
        let m = genot.m();
        let mut mafs = vec![f64::NAN; m];
        mafs.iter_mut()
            .zip(genot.iter_snv())
            .par_bridge()
            .for_each(|(maf, g_snv)| *maf = g_snv.maf());

        self.snvs_mut().set_maf(mafs)
    }

    pub fn new_score<W: WgtTrait>(
        fin: &Path,
        gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        //fin_phe: Option<&Path>,
        phe_name: Option<&str>,
        cov_name: Option<&str>,
        extract_sample_buf: Option<&[u8]>,
        //fin_sample: Option<&Path>,
        //fin_cov: Option<&Path>,
        wgts: &mut [W],
        //wgts: &[WgtBoost],
        use_missing: bool, // use WgtBoosts and use wgts.use_missing()
    ) -> Self {
        let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        log::debug!("m_in {}", m_in);
        let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        log::debug!("n_in {}", n_in);

        let (n, use_samples) = sample::make_use_samples_buf(extract_sample_buf, fin, gfmt);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        // TODO: make (string, string) as key of hashmap
        //let samples_id = plink::load_samples_id(fin, &use_samples);

        let sample_id_to_n = samples::create_sample_id_to_n(fin, gfmt, Some(&use_samples));
        let ys: Vec<bool> = io_genot::load_ys_buf(fin, gfmt, phe_buf, phe_name, &sample_id_to_n);
        //let ys: Vec<bool> = io_genot::load_ys(fin, gfmt, fin_phe, phe_name, &use_samples);
        let covs: Option<Covs> = match cov_name {
            Some(cov_name) => Some(Covs::new(phe_buf, fin, gfmt, cov_name, &sample_id_to_n)),
            None => None,
        };
        //let covs: Option<Covs> = Covs::new_old(fin_cov, cov_name, &sample_id_to_n);
        //// covs could not exist
        //let covs: Option<Vec<Var>> =
        //    cov::load_vars(fin_cov, n_in, &sample_id_to_n, cov::CovKind::Cov);

        let samples = Samples::new(&ys, None, covs, n);

        //let mut wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(fin_wgt);
        //wgt_boost::io::set_covs(&mut wgts, covs.as_deref(), n);

        // set genotype index in wgt
        //let use_missing = true;
        let genot = io_genot::load_score::load_genotypes_for_score(
            fin,
            gfmt,
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
        gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        //fin_phe: Option<&Path>,
        phe_name: Option<&str>,
        cov_name: Option<&str>,
        extract_sample_buf: Option<&[u8]>,
        //fin_sample: Option<&Path>,
        //fin_cov: Option<&Path>,
        // TODO: &mut[&mut[Wgt]],
        wgts: &mut [Wgts],
        //wgts: &[WgtBoost],
    ) -> Dataset {
        let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        log::debug!("m_in {}", m_in);
        let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        log::debug!("n_in {}", n_in);

        let (n, use_samples) = sample::make_use_samples_buf(extract_sample_buf, fin, gfmt);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        // TODO: make (string, string) as key of hashmap
        //let samples_id = plink::load_samples_id(fin, &use_samples);

        let sample_id_to_n = samples::create_sample_id_to_n(fin, gfmt, Some(&use_samples));

        let ys: Vec<bool> = io_genot::load_ys_buf(fin, gfmt, phe_buf, phe_name, &sample_id_to_n);
        //let ys: Vec<bool> = io_genot::load_ys(fin, gfmt, fin_phe, phe_name, &use_samples);
        let covs: Option<Covs> = match cov_name {
            Some(cov_name) => Some(Covs::new(phe_buf, fin, gfmt, cov_name, &sample_id_to_n)),
            None => None,
        };
        //let ys: Vec<bool> = io_genot::load_ys(fin, gfmt, fin_phe, phe_name, &use_samples);
        //let covs: Option<Covs> = Covs::new_old(fin_cov, cov_name, &sample_id_to_n);
        //// covs could not exist
        //let covs: Option<Vec<Var>> =
        //    cov::load_vars(fin_cov, n_in, &sample_id_to_n, cov::CovKind::Cov);

        let samples = Samples::new(&ys, None, covs, n);

        //let mut wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(fin_wgt);
        //wgt_boost::io::set_covs(&mut wgts, covs.as_deref(), n);

        // set genotype index in wgt
        // TODO: argparse
        let use_missing = false;
        let genot = io_genot::load_score::load_genotypes_for_score_multiwgts(
            fin,
            gfmt,
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
