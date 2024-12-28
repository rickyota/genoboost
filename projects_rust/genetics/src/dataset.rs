//! Dataset
//! It is important that all vector values can be extract without constructing. ex. mafs: Vec<f64> not Vec<Info>
//! Must use fn to access data so that it is easy to use Trait
//! For now, do not use generics, use in dataset_future.rs
//!
//!

//pub mod genot_io;
//pub mod samples;
//pub mod snvs;

use rayon::iter::{ParallelBridge, ParallelIterator};
use std::collections::HashMap;
use std::time::Instant;

use crate::dataset_file::DatasetFile;
use crate::genot::prelude::*;
use crate::genot_io;
use crate::{sample, snv, vec, wgt::WgtTrait, GenotFile, Samples, Snvs, Wgts};

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

    // add assertion for training
    pub fn new_datasetfile_training(
        dfile: &DatasetFile, // should be immutable for simplicity
        use_sample_val: bool,
        filt_snv: Option<&[bool]>,
        fill_missing: bool,
        make_major_a2_train: bool,
        snvs_train: Option<&Snvs>, //for use_missing of vali
        //fill_missing: bool,
        //filt_snv: Option<&[bool]>,
        //snvs_train: Option<&Snvs>, //for use_missing of vali
        //make_major_a2_train: bool,
        mem: Option<usize>,
    ) -> Self {
        if fill_missing && use_sample_val && snvs_train.is_none() {
            panic!("fill_missing and use_sample_val require snvs_train for maf");
        }

        let dataset = Self::new_datasetfile(
            dfile,
            use_sample_val,
            filt_snv,
            fill_missing,
            make_major_a2_train,
            snvs_train,
            mem,
        );

        // check if ys exist
        if dataset.samples().phe().is_none() {
            panic!("Could not load sample phenotype.");
        }

        dataset
    }

    fn new_datasetfile(
        dfile: &DatasetFile, // should be immutable for simplicity
        use_sample_val: bool,
        filt_snv: Option<&[bool]>,
        fill_missing: bool,
        make_major_a2_train: bool,
        snvs_train: Option<&Snvs>,
        mem: Option<usize>,
    ) -> Self {
        let sample_buf = if use_sample_val {
            dfile.sample_val_buf()
        } else {
            dfile.sample_buf()
        };
        Self::new(
            dfile.fin_genot(),
            //dfile.fin(),
            //dfile.gfmt(),
            dfile.phe_buf(),
            dfile.phe_name(),
            dfile.cov_name(),
            dfile.snv_buf(),
            dfile.snv_set_buf(),
            sample_buf,
            //dfile.sample_buf(),
            filt_snv,
            fill_missing,
            make_major_a2_train,
            snvs_train,
            mem,
        )
    }

    /// Input is one of the following
    /// 1. plink2 + fin_phe
    /// 2. plink2 (phe) + fin_phe (cov) : phe_name is in psam
    /// 3. plink2 (cov and phe in .psam) : phe_buf is option
    /// 4. plink1 + fin_phe (cov and phe)
    /// 5. plink1 (phe) + fin_phe (cov) : phe_name is option
    fn new(
        fin_genot: &GenotFile,
        //fin: &Path,
        //gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        phe_name: Option<&str>,
        cov_name: Option<&str>,
        snv_buf: Option<&[u8]>,
        group_snv_buf: Option<&[u8]>,
        sample_buf: Option<&[u8]>,
        // TODO: merge filt_snv and fin_snv into use_snvs,
        filt_snv: Option<&[bool]>,
        fill_missing: bool,
        make_major_a2_train: bool,
        snvs_train: Option<&Snvs>, //for use_missing of vali
        mem: Option<usize>,
    ) -> Self {
        fin_genot.check_valid_open();
        //genot_io::check_valid_fin(fin, gfmt);

        //let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        //log::debug!("m_in: {}", m_in);
        //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        //log::debug!("n_in: {}", n_in);
        // load snvs
        let snvs_in = genot_io::load_snvs(fin_genot);
        //let snvs_in = plink::load_snvs(fin, m_in);
        let (mut use_snvs, mut m_snv) = snv::make_use_snvs_buf(snv_buf, &snvs_in);
        log::debug!("m_snv: {}", m_snv);

        // filter snv only for use_snvs not for set_snvs
        if let Some(filt_snv_) = filt_snv {
            log::debug!("filt_snv before m: {}", m_snv);
            assert_eq!(filt_snv_.len(), m_snv);
            use_snvs = vec::and_in_bool_vec(&use_snvs, filt_snv_);
            m_snv = vec::count_true(&use_snvs);
            log::debug!("filt_snv m_snv: {}", m_snv);
        }

        // load snv set
        //let (use_set_snvs, m_set) = snv::make_set_snvs_buf(set_snv_buf, &snvs_in);
        //assert_eq!(use_set_snvs.len(), use_snvs.len());
        //log::debug!("m_set: {}", m_set);

        let (group_to_m_in, m_set) = snv::make_group_to_m_in_buf(group_snv_buf, &snvs_in);
        //let set_snvs_use = snv::load_set_snvs_buf(set_snv_buf.unwrap());
        //let set_to_m_in = snv::load_set_to_m_in(set_snvs_use, &snvs_in);

        // later
        //if m_set != 0 {
        //    //  set snvs exist
        //    use_snvs = vec::or_bool_vec(&use_snvs, &use_set_snvs);
        //    m_snv = vec::count_true(&use_snvs);
        //    log::debug!("m_snv after merge with snv set: {}", m_snv);
        //}

        //if m_snv == 0 {
        if (m_snv == 0) && (m_set == 0) {
            panic!("Using SNVs are zero. Please check fin_snv.")
        }
        log::debug!("m_snv: {}", m_snv);

        let (use_samples, n) = sample::make_use_samples_buf(sample_buf, fin_genot);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }
        log::debug!("n: {}", n);

        Self::new_use_vec(
            fin_genot,
            phe_buf,
            phe_name,
            cov_name,
            group_snv_buf,
            Some(&use_snvs),
            //Some(&use_set_snvs),
            Some(&use_samples),
            fill_missing,
            make_major_a2_train,
            snvs_train,
            mem,
            group_to_m_in,
        )
    }

    pub fn new_use_vec(
        fin_genot: &GenotFile,
        phe_buf: Option<&[u8]>,
        phe_name: Option<&str>,
        cov_name: Option<&str>,
        set_snv_buf: Option<&[u8]>,
        use_snvs: Option<&[bool]>,
        //use_set_snvs: Option<&[bool]>,
        use_samples: Option<&[bool]>,
        fill_missing: bool,
        make_major_a2_train: bool,
        snvs_train: Option<&Snvs>, //for fill_missing and make_major_a2_train of vali
        mem: Option<usize>,
        group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    ) -> Self {
        let start = Instant::now();

        //let use_snvs_or_set = match use_snvs {
        //    Some(use_snvs) => match use_set_snvs {
        //        // use both
        //        Some(use_set_snvs) => Some(vec::or_bool_vec(use_snvs, use_set_snvs)),
        //        // use_snvs only
        //        None => Some(use_snvs.to_vec()),
        //    },
        //    // use all snvs
        //    None => None,
        //};

        let m = match use_snvs {
            Some(x) => vec::count_true(x),
            None => genot_io::compute_num_snv(fin_genot).unwrap(),
        };

        let n = match use_samples {
            Some(x) => vec::count_true(x),
            None => genot_io::compute_num_sample(fin_genot).unwrap(),
        };

        // do not fill missing here since validation dataset requires maf to fill missing; below
        let genot = genot_io::load::generate_genot(
            fin_genot,
            m,
            n,
            use_snvs,
            //use_snvs_or_set.as_deref(),
            use_samples,
            false,
            mem,
            group_to_m_in,
        );
        //genot_io::load::generate_genot(fin_genot, m, n, use_snvs, use_samples, false, mem);

        let samples = Samples::new_plink(
            fin_genot,
            phe_buf,
            phe_name,
            cov_name,
            use_samples,
            Some(n),
            true,
        );

        let snvs = Snvs::new_plink_use_snvs_and_set(fin_genot, use_snvs, None, set_snv_buf);
        // Snvs::new_plink_use_snvs_and_set(fin_genot, use_snvs, None, set_snv_buf.unwrap());
        //let snvs = Snvs::new_plink_use_snvs(fin_genot, use_snvs, Some(m_snv));

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

        // fill missing to mode here
        // for training; use maf assuming hwe
        // for validation; use the training maf
        if fill_missing {
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

        if make_major_a2_train {
            if let Some(snvs_train) = snvs_train {
                // for val
                dataset.set_major_a2(Some(snvs_train));
            } else {
                // for training
                dataset.set_major_a2(None);
            }
        }

        dataset
    }

    // for test
    pub fn new_field_phe(genot: Genot, samples: Samples, snvs: Snvs) -> Self {
        Self {
            genot,
            snvs,
            samples,
        }
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

    pub fn fill_missing_maf(&mut self, mafs: Option<&Vec<f64>>) {
        let mafs_v: Vec<f64>;
        let mafs = if let Some(mafs) = mafs {
            // for validation
            mafs
        } else {
            // for training
            // if not using clone(), genot_mut() is error
            mafs_v = self.snvs().mafs().unwrap().clone();
            &mafs_v
            // self.snvs().mafs().unwrap()
        };
        //let genot = self.genot_mut();
        self.genot_mut()
            .iter_snv_mut()
            .zip(mafs.iter())
            .par_bridge()
            .for_each(|(mut g_snv, maf)| g_snv.fill_missing_mode_maf(*maf));
        //.for_each(|(mut g_snv, maf)| io_genot::load::fill_missing_maf(&mut g_snv, *maf));
    }

    pub fn set_major_a2(&mut self, snvs: Option<&Snvs>) {
        // 1. set major as a2 in snvs
        // 2. same in genot
        // 3. update maf

        let mafs_v: Vec<f64>;
        let mafs = if let Some(snvs) = snvs {
            // for validation
            snvs.mafs().unwrap()
        } else {
            // for training
            mafs_v = self.snvs().mafs().unwrap().clone();
            &mafs_v
        };

        // reverse snvs
        self.snvs_mut()
            .snv_ids_mut()
            .iter_mut()
            .zip(mafs.iter())
            .filter(|(_, maf)| **maf > 0.5)
            .par_bridge()
            .for_each(|(snv_id, _)| snv_id.reverse_alleles());

        // reverse genot
        self.genot_mut()
            .iter_snv_mut()
            .zip(mafs.iter())
            .filter(|(_, maf)| **maf > 0.5)
            .par_bridge()
            .for_each(|(mut g_snv, _)| g_snv.reverse_allele());

        // update maf
        self.compute_maf();

        // check all maf<0.5 for training?
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

    pub fn new_datasetfile_score<W: WgtTrait>(
        dfile: &DatasetFile,
        wgts: &mut [W],
        //wgts: &[WgtBoost],
        fill_missing_in_dataset: bool,
        //fill_missing: bool,
        allow_nonexist_snv: bool,
        use_snv_pos: bool,
        mem: Option<usize>,
    ) -> Self {
        Self::new_score(
            dfile.fin_genot(),
            //dfile.fin(),
            //dfile.gfmt(),
            dfile.phe_buf(),
            dfile.cov_name(),
            dfile.sample_buf(),
            wgts,
            fill_missing_in_dataset,
            //fill_missing,
            allow_nonexist_snv,
            use_snv_pos,
            mem,
        )
    }

    // for boosting
    // merge partly to new()
    pub fn new_score<W: WgtTrait>(
        fin_genot: &GenotFile,
        //fin: &Path,
        //gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        //fin_phe: Option<&Path>,
        //phe_name: Option<&str>,
        cov_name: Option<&str>,
        sample_buf: Option<&[u8]>,
        //extract_sample_buf: Option<&[u8]>,
        //fin_sample: Option<&Path>,
        //fin_cov: Option<&Path>,
        wgts: &mut [W],
        //wgts: &[WgtBoost],
        fill_missing_in_dataset: bool,
        //use_missing: bool, // use WgtBoosts and use wgts.use_missing()
        // fill_missing: bool,
        allow_nonexist_snv: bool,
        use_snv_pos: bool,
        mem: Option<usize>,
    ) -> Self {
        //let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        //log::debug!("m_in {}", m_in);
        //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        //log::debug!("n_in {}", n_in);

        let (use_samples, n) = sample::make_use_samples_buf(sample_buf, fin_genot);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        let samples = Samples::new_plink(
            fin_genot,
            phe_buf,
            None,
            cov_name,
            Some(&use_samples),
            Some(n),
            false,
        );

        // TODO: make (string, string) as key of hashmap
        //let samples_id = plink::load_samples_id(fin, &use_samples);

        //let sample_id_to_n = samples::create_sample_id_to_n(fin, gfmt, Some(&use_samples));
        ////let ys: Vec<bool> = io_genot::load_ys_buf(fin, gfmt, phe_buf, phe_name, &sample_id_to_n);
        ////let ys: Vec<bool> = io_genot::load_ys(fin, gfmt, fin_phe, phe_name, &use_samples);

        //let covs: Option<Covs> =
        //    cov_name.map(|x| Covs::new(phe_buf, fin, gfmt, x, &sample_id_to_n));
        ////let covs: Option<Covs> = match cov_name {
        ////    Some(cov_name) => Some(Covs::new(phe_buf, fin, gfmt, cov_name, &sample_id_to_n)),
        ////    None => None,
        ////};
        ////let covs: Option<Covs> = Covs::new_old(fin_cov, cov_name, &sample_id_to_n);
        ////// covs could not exist
        ////let covs: Option<Vec<Var>> =
        ////    cov::load_vars(fin_cov, n_in, &sample_id_to_n, cov::CovKind::Cov);

        //let samples_id = genot_io::load_samples_id(fin, gfmt, Some(&use_samples));

        //let sample_id_to_n = samples::create_sample_id_to_n(fin_genot, Some(&use_samples));
        ////let ys: Vec<bool> = io_genot::load_ys_buf(fin, gfmt, phe_buf, phe_name, &sample_id_to_n);
        ////let ys: Vec<bool> = io_genot::load_ys(fin, gfmt, fin_phe, phe_name, &use_samples);

        //let covs: Option<Covs> =
        //    cov_name.map(|x| Covs::new(phe_buf, fin_genot, x, &sample_id_to_n));
        ////let covs: Option<Covs> = match cov_name {
        ////    Some(cov_name) => Some(Covs::new(phe_buf, fin, gfmt, cov_name, &sample_id_to_n)),
        ////    None => None,
        ////};
        ////let covs: Option<Covs> = Covs::new_old(fin_cov, cov_name, &sample_id_to_n);
        ////// covs could not exist
        ////let covs: Option<Vec<Var>> =
        ////    cov::load_vars(fin_cov, n_in, &sample_id_to_n, cov::CovKind::Cov);

        //let samples_id = genot_io::load_samples_id(fin_genot, Some(&use_samples));

        //let samples = Samples::new_nophe(Some(samples_id), covs, n);
        //let samples = Samples::new(&ys, None, covs, n);

        //let mut wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(fin_wgt);
        //wgt_boost::io::set_covs(&mut wgts, covs.as_deref(), n);

        // set genotype index in wgt
        //let use_missing = true;
        let genot = genot_io::load_score::load_genotypes_for_score(
            fin_genot,
            //gfmt,
            wgts,
            n,
            Some(&use_samples),
            // fill_missing,
            fill_missing_in_dataset,
            allow_nonexist_snv,
            use_snv_pos,
            mem,
        );

        Dataset {
            genot,
            // unnecessary since index is in WgtKInd
            snvs: Snvs::new_empty(),
            samples,
        }
    }

    /// For multiple wgts
    pub fn new_datasetfile_score_genetics(
        dfile: &DatasetFile,
        wgts: &mut [Wgts],
        fill_missing_in_dataset: bool,
        allow_nonexist_snv: bool,
        use_snv_pos: bool,
        mem: Option<usize>,
    ) -> Self {
        Self::new_score_genetics(
            dfile.fin_genot(),
            //fin_genot,
            dfile.phe_buf(),
            //phe_buf,
            dfile.cov_name(),
            //cov_name,
            dfile.sample_buf(),
            //sample_buf,
            dfile.freq_buf(),
            wgts,
            fill_missing_in_dataset,
            allow_nonexist_snv,
            use_snv_pos,
            mem,
        )
    }

    // Usually not fill missing here. fill_missing_in_dataset is for backward compatibility.
    // for genetics::score()
    pub fn new_score_genetics(
        fin_genot: &GenotFile,
        //fin: &Path,
        //gfmt: GenotFormat,
        phe_buf: Option<&[u8]>,
        cov_name: Option<&str>,
        sample_buf: Option<&[u8]>,
        freq_buf: Option<&[u8]>,
        wgts: &mut [Wgts],
        //wgts: &[WgtBoost],
        fill_missing_in_dataset: bool,
        allow_nonexist_snv: bool,
        use_snv_pos: bool,
        mem: Option<usize>,
    ) -> Self {
        //let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
        //log::debug!("m_in {}", m_in);
        //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
        //log::debug!("n_in {}", n_in);

        //let (n, use_samples) = sample::make_use_samples_buf(sample_buf, fin, gfmt);
        let (use_samples, n) = sample::make_use_samples_buf(sample_buf, fin_genot);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
        if n == 0 {
            panic!("Using samples are zero. Please check fin_sample.")
        }

        // TODO: make (string, string) as key of hashmap
        //let samples_id = plink::load_samples_id(fin, &use_samples);

        let samples = Samples::new_plink(
            fin_genot,
            phe_buf,
            None,
            cov_name,
            Some(&use_samples),
            Some(n),
            false,
        );

        // set genotype index in wgt
        // TODO: argparse
        //let use_missing = false;
        // TOFIX: fill_missing using MAF (in wgt or) provided.
        //let fill_missing_in_test = true;
        let genot = genot_io::load_score::load_genotypes_for_score_multiwgts(
            fin_genot,
            wgts,
            //&mut wgts,
            freq_buf,
            n,
            Some(&use_samples),
            fill_missing_in_dataset,
            allow_nonexist_snv,
            use_snv_pos,
            mem,
        );

        // add freq to wgts
        //let freqs_in = genot_io::load_freq(freq_buf);
        // index
        // is_reversed

        // TMP
        //for snv in genot.iter_snv() {
        //    println!("snv {:?}", &snv.vals()[..10]);
        //}

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
