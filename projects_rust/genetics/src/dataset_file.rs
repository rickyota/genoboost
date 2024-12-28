//! Dataset File Struct

//use crate::genot_io;
use crate::genot_io::GenotFile;
use crate::textfile;
//use crate::GenotFormat;
use std::path::{Path, PathBuf};

// FIXME: not ideal for cross_validations since sample_buf is changed every time
// how about implement update_sample_buf()? -> requires copy of toher buf, but work fine.
#[derive(Clone)]
pub struct DatasetFile {
    fin_genot: GenotFile,
    fin_phe: Option<PathBuf>,
    phe_name: Option<String>,
    cov_name: Option<String>,
    fin_snv: Option<PathBuf>,
    fin_sample: Option<PathBuf>,
    //fin_cov: Option<&Path>,
    fin_sample_val: Option<PathBuf>,
    fin_score_start: Option<PathBuf>,
    fin_initial_snv: Option<PathBuf>,
    fin_set_snv: Option<PathBuf>,
    fin_freq: Option<PathBuf>,
    //option
    // could make mistakes
    //make_major_a2_train: bool,
    // buffer
    phe_buf: Option<Vec<u8>>,
    snv_buf: Option<Vec<u8>>,
    set_snv_buf: Option<Vec<u8>>,
    // sample_buf might not correspond to fin_sample for cross validation
    sample_buf: Option<Vec<u8>>,
    sample_val_buf: Option<Vec<u8>>,
    // unnecessary -> necessary for pipe input
    score_start_buf: Option<Vec<u8>>,
    freq_buf: Option<Vec<u8>>,
}

impl DatasetFile {
    pub fn check_valid_fin(&self) {
        //genot_io::check_valid_fin(self.fin(), self.gfmt());
        self.fin_genot().check_valid_open();
    }

    // FIXME: check if in .fam or .phe, col=0 is not duplicated
    //pub fn check_sample(&self){ }

    // use genot_io::compute_num_sample() instead
    //pub fn compute_num_sample(&self) -> Option<usize> {
    //    genot_io::compute_num_sample(self.fin(), self.gfmt())
    //}

    // use genot_io::load_samples_id() instead
    //pub fn load_samples_id(&self) -> Vec<String> {
    //    genot_io::load_samples_id(self.fin(), self.gfmt(), None)
    //}

    pub fn new(
        fin_genot: GenotFile,
        fin_phe: Option<PathBuf>,
        phe_name: Option<String>,
        cov_name: Option<String>,
        fin_snv: Option<PathBuf>,
        fin_sample: Option<PathBuf>,
        fin_sample_val: Option<PathBuf>,
        //make_major_a2_train: bool,
    ) -> Self {
        Self {
            fin_genot,
            fin_phe,
            phe_name,
            cov_name,
            fin_snv,
            fin_sample,
            fin_sample_val,
            fin_score_start: None,
            fin_initial_snv: None,
            fin_set_snv: None,
            fin_freq: None,
            //make_major_a2_train,
            phe_buf: None,
            snv_buf: None,
            set_snv_buf: None,
            sample_buf: None,
            sample_val_buf: None,
            score_start_buf: None,
            freq_buf: None,
        }
    }

    // avoid clone; use update_file_initial_snv()
    //pub fn set_file_initial_snv(self, fin_initial_snv: Option<PathBuf>) -> Self {
    //    let mut dfile = self.clone();
    //    dfile.update_file_initial_snv(fin_initial_snv);
    //    dfile
    //}

    pub fn update_file_freq(&mut self, fin_freq: Option<PathBuf>) {
        self.fin_freq = fin_freq;
    }

    pub fn update_file_group_snv(&mut self, fin_set_snv: Option<PathBuf>) {
        self.fin_set_snv = fin_set_snv;
    }

    pub fn update_file_score_start(&mut self, fin_score_start: Option<PathBuf>) {
        self.fin_score_start = fin_score_start;
    }

    pub fn update_file_initial_snv(&mut self, fin_initial_snv: Option<PathBuf>) {
        self.fin_initial_snv = fin_initial_snv;
    }

    pub fn fin_genot(&self) -> &GenotFile {
        &self.fin_genot
    }
    pub fn fin_phe(&self) -> Option<&Path> {
        self.fin_phe.as_ref().map(|x| x.as_path())
    }
    pub fn phe_name(&self) -> Option<&str> {
        self.phe_name.as_ref().map(|x| x.as_str())
    }
    pub fn cov_name(&self) -> Option<&str> {
        self.cov_name.as_ref().map(|x| x.as_str())
    }
    pub fn fin_snv(&self) -> Option<&Path> {
        self.fin_snv.as_ref().map(|x| x.as_path())
    }
    pub fn fin_snv_set(&self) -> Option<&Path> {
        self.fin_set_snv.as_ref().map(|x| x.as_path())
    }
    pub fn fin_sample(&self) -> Option<&Path> {
        self.fin_sample.as_ref().map(|x| x.as_path())
    }
    // DO NOT use fin_sample_val() to check if is_monitor is true or false;
    // use sample_val_buf() instead
    // since sample_val_buf could be not None even when fin_sample_val is None
    pub fn fin_sample_val(&self) -> Option<&Path> {
        self.fin_sample_val.as_ref().map(|x| x.as_path())
    }
    pub fn fin_score_start(&self) -> Option<&Path> {
        self.fin_score_start.as_ref().map(|x| x.as_path())
    }

    pub fn fin_initial_snv(&self) -> Option<&Path> {
        self.fin_initial_snv.as_ref().map(|x| x.as_path())
    }

    pub fn fin_freq(&self) -> Option<&Path> {
        self.fin_freq.as_ref().map(|x| x.as_path())
    }

    //pub fn make_major_a2_train(&self) -> bool {
    //    self.make_major_a2_train
    //}
    pub fn phe_buf(&self) -> Option<&[u8]> {
        self.phe_buf.as_deref()
    }
    pub fn snv_buf(&self) -> Option<&[u8]> {
        self.snv_buf.as_deref()
    }
    pub fn snv_set_buf(&self) -> Option<&[u8]> {
        self.set_snv_buf.as_deref()
    }
    pub fn sample_buf(&self) -> Option<&[u8]> {
        self.sample_buf.as_deref()
    }
    pub fn sample_val_buf(&self) -> Option<&[u8]> {
        self.sample_val_buf.as_deref()
    }
    pub fn score_start_buf(&self) -> Option<&[u8]> {
        self.score_start_buf.as_deref()
    }
    pub fn freq_buf(&self) -> Option<&[u8]> {
        self.freq_buf.as_deref()
    }

    // Assume to call read_...() only once at first, and make unmutable.
    pub fn reads(&mut self) {
        self.read_phe();
        self.read_snv();
        self.read_snv_set();
        self.read_sample();
        self.read_sample_val();
        self.read_score_start();
        self.read_freq();
    }

    // do not call directly
    fn read_phe(&mut self) {
        if self.phe_buf.is_none() {
            if let Some(fin_phe) = self.fin_phe() {
                let phe_buf = textfile::read_file_to_end(fin_phe, None).unwrap_or_else(|_| {
                    panic!("Cannot read file: {:?}", fin_phe);
                });
                self.phe_buf = Some(phe_buf);
            }
        }
        // otherwise, do nothing
    }

    fn read_snv(&mut self) {
        if self.snv_buf.is_none() {
            if let Some(fin_snv) = self.fin_snv() {
                let snv_buf = textfile::read_file_to_end(fin_snv, None).unwrap_or_else(|_| {
                    panic!("Cannot read file: {:?}", fin_snv);
                });
                self.snv_buf = Some(snv_buf);
            }
        }
        // otherwise, do nothing
    }

    fn read_snv_set(&mut self) {
        if self.set_snv_buf.is_none() {
            if let Some(fin_snv_set) = self.fin_snv_set() {
                let snv_set_buf =
                    textfile::read_file_to_end(fin_snv_set, None).unwrap_or_else(|_| {
                        panic!("Cannot read file: {:?}", fin_snv_set);
                    });
                self.set_snv_buf = Some(snv_set_buf);
            }
        }
        // otherwise, do nothing
    }

    fn read_sample(&mut self) {
        if self.sample_buf.is_none() {
            if let Some(fin_sample) = self.fin_sample() {
                let sample_buf =
                    textfile::read_file_to_end(fin_sample, None).unwrap_or_else(|_| {
                        panic!("Cannot read file: {:?}", fin_sample);
                    });
                self.sample_buf = Some(sample_buf);
            }
        }
        // otherwise, do nothing
    }

    fn read_sample_val(&mut self) {
        if self.sample_val_buf.is_none() {
            if let Some(fin_sample_val) = self.fin_sample_val() {
                let sample_val_buf = textfile::read_file_to_end(fin_sample_val, None)
                    .unwrap_or_else(|_| {
                        panic!("Cannot read file: {:?}", fin_sample_val);
                    });
                self.sample_val_buf = Some(sample_val_buf);
            }
        }
    }

    fn read_freq(&mut self) {
        if self.freq_buf.is_none() {
            if let Some(fin_freq) = self.fin_freq() {
                let freq_buf = textfile::read_file_to_end(fin_freq, None).unwrap_or_else(|_| {
                    panic!("Cannot read file: {:?}", fin_freq);
                });
                self.freq_buf = Some(freq_buf);
            }
        }
        // otherwise, do nothing
    }

    pub fn update_sample_buf(&mut self, sample_buf: Vec<u8>, sample_val_buf: Vec<u8>) {
        self.sample_buf = Some(sample_buf);
        self.sample_val_buf = Some(sample_val_buf);
    }

    //pub fn update_sample_buf(&self, sample_buf: Vec<u8>, sample_val_buf: Vec<u8>) -> Self {
    //    let mut dfile = self.clone();
    //    dfile.sample_buf = Some(sample_buf);
    //    dfile.sample_val_buf = Some(sample_val_buf);
    //    dfile
    //}

    // TODO: merge with score:score_struct:new_from()
    // or create load_score()
    // 
    // directly run in boosting_batch()
    fn read_score_start(&mut self) {
        if self.score_start_buf.is_none() {
            if let Some(fin_score_start) = self.fin_score_start() {
                let score_start_buf = if fin_score_start.extension().unwrap() == "zst" {
                    textfile::read_file_to_end(fin_score_start, Some("zst")).unwrap_or_else(|_| {
                        panic!("Cannot read file: {:?}", fin_score_start);
                    })
                } else {
                    textfile::read_file_to_end(fin_score_start, None).unwrap_or_else(|_| {
                        panic!("Cannot read file: {:?}", fin_score_start);
                    })
                };
                self.score_start_buf = Some(score_start_buf);
            }
        }
    }
}
