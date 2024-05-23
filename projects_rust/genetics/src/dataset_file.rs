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
    //option
    // could create mistakes
    //make_major_a2_train: bool,
    // buffer
    phe_buf: Option<Vec<u8>>,
    snv_buf: Option<Vec<u8>>,
    // sample_buf might not correspond to fin_sample for cross validation
    sample_buf: Option<Vec<u8>>,
    sample_val_buf: Option<Vec<u8>>,
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
        //args.file_score_start,
        fin_score_start: Option<PathBuf>,
    ) -> Self {
        Self {
            fin_genot,
            fin_phe,
            phe_name,
            cov_name,
            fin_snv,
            fin_sample,
            fin_sample_val,
            fin_score_start,
            fin_initial_snv: None,
            //make_major_a2_train,
            phe_buf: None,
            snv_buf: None,
            sample_buf: None,
            sample_val_buf: None,
        }
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
    pub fn fin_sample(&self) -> Option<&Path> {
        self.fin_sample.as_ref().map(|x| x.as_path())
    }
    pub fn fin_sample_val(&self) -> Option<&Path> {
        self.fin_sample_val.as_ref().map(|x| x.as_path())
    }
    pub fn fin_score_start(&self) -> Option<&Path> {
        self.fin_score_start.as_ref().map(|x| x.as_path())
    }

    pub fn fin_initial_snv(&self) -> Option<&Path> {
        self.fin_initial_snv.as_ref().map(|x| x.as_path())
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
    pub fn sample_buf(&self) -> Option<&[u8]> {
        self.sample_buf.as_deref()
    }
    pub fn sample_val_buf(&self) -> Option<&[u8]> {
        self.sample_val_buf.as_deref()
    }

    // Assume to call read_...() only once at first, and make unmutable.
    pub fn reads(&mut self) {
        self.read_phe();
        self.read_snv();
        self.read_sample();
        self.read_sample_val();
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
}
