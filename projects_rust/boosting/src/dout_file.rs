//! Dataset File Struct

use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use crate::BoostType;
use genetics::textfile;

/// For wgt output dir
/// For dout_cv: use fn dout_cv() -> DoutFile
#[derive(Clone, Debug)]
pub struct DoutFile {
    dout: PathBuf,
}

impl DoutFile {
    pub fn new(dout: PathBuf) -> Self {
        Self { dout }
    }

    pub fn dout(&self) -> &Path {
        &self.dout
    }

    pub fn dout_cv(&self, cvi: usize) -> Self {
        let mut d = self.dout().to_owned();
        let dpara = String::from("cv-") + &cvi.to_string();
        d.push(dpara);

        Self::new(d)
    }

    pub fn dout_para_integrate(&self, learning_rate: f64, boost_type: BoostType) -> DoutParaFile {
        // TODO: use .join
        let mut d = self.dout().to_owned();
        let mut dpara = String::from("para");
        dpara += &(String::from("_type-") + &boost_type.genetic_model());
        //if let Some(lr) = learning_rate {
        dpara += &(String::from("_lr-") + &learning_rate.to_string());
        //dpara += &(String::from(".lr") + &learning_rate.to_string());
        //}
        d.push(dpara);

        DoutParaFile::new(d)
    }

    // when increase para, add them here as argument
    pub fn dout_para_lr(&self, learning_rate: f64) -> DoutParaFile {
        // TODO: use .join
        let mut d = self.dout.to_owned();
        //let mut d = dout.to_owned();
        let mut dpara = String::from("para");
        dpara += &(String::from("_lr-") + &learning_rate.to_string());
        d.push(dpara);

        DoutParaFile::new(d)
        //DoutParaFile::new(self.get_dname_para(learning_rate))
    }

    // when increase para, add them here as argument
    //fn get_dname_para(&self, learning_rate: f64) -> PathBuf {
    //    // TODO: use .join
    //    let mut d = self.dout.to_owned();
    //    //let mut d = dout.to_owned();
    //    let mut dpara = String::from("para");
    //    dpara += &(String::from("_lr-") + &learning_rate.to_string());
    //    d.push(dpara);
    //    d
    //}

    // use this in DoutParaFile
    //pub fn get_file_wgt(&self, learning_rate: f64) -> PathBuf {
    //    let mut p = self.get_dname_para(learning_rate);
    //    p.push("boosting.wgt");
    //    p

    //    //get_fname_wgt(&p)
    //}

    pub fn get_fname_wgt_best_para(&self) -> PathBuf {
        let mut p = self.dout().to_owned();
        p.push("boosting.wgt.para");
        p
    }

    pub fn get_fname_wgt_best(&self) -> PathBuf {
        let mut p = self.dout().to_owned();
        p.push("boosting.wgt");
        p
    }

    //pub fn check_file_wgt_not_exist(&self, learning_rate: f64) {
    //    let fwgt = self.get_file_wgt(learning_rate);
    //    //let fwgt = get_file_wgt(dout, learning_rate);
    //    let exist_fwgt = textfile::exist_file(&fwgt);
    //    if exist_fwgt {
    //        panic!(
    //            "Weight file already exists: {:?}. Delete it or use --resume option.",
    //            &fwgt
    //        );
    //    }
    //}

    pub fn create_dout(&self) {
        create_dir(self.dout());
    }
}

// too complicated
//#[derive(Clone)]
//pub struct DoutCvFile {
//    dout_cv: PathBuf,
//}
//
//pub struct DoutIntegrateFile {
//    dout_integrate: PathBuf,
//}
//
//trait DoutPara {}

// pass this to boosting_train::boosting()
#[derive(Clone)]
pub struct DoutParaFile {
    dout_para: PathBuf,
}

impl DoutParaFile {
    pub fn new(dout_para: PathBuf) -> Self {
        Self { dout_para }
    }

    pub fn dout_para(&self) -> &Path {
        &self.dout_para
    }

    pub fn get_file_wgt(&self) -> PathBuf {
        let mut p = self.dout_para().to_owned();
        p.push("boosting.wgt");
        p

        //get_fname_wgt(&p)
    }

    pub fn check_file_wgt_not_exist(&self) {
        let fwgt = self.get_file_wgt();

        let exist_fwgt = textfile::exist_file(&fwgt);
        if exist_fwgt {
            panic!(
                "Weight file already exists: {:?}. Delete it or use --resume option.",
                &fwgt
            );
        }
    }

    pub fn create_dout_para(&self) {
        create_dir(self.dout_para());
    }

    pub fn bufwriter_fwgt_append(&self) -> BufWriter<File> {
        let fwgt = self.get_file_wgt();
        create_dir(&fwgt.parent().unwrap());
        textfile::bufwriter_append(&fwgt)

        //let f = match OpenOptions::new()
        //    .write(true)
        //    .append(true)
        //    .create(true)
        //    .open(&fwgt)
        //{
        //    Ok(f) => f,
        //    Err(e) => panic!("Error opening file {:?}: {}", &fwgt, e),
        //};
        //BufWriter::new(f)
    }

    pub fn floss(&self, ti: usize) -> PathBuf {
        let mut f = self.dout_para().to_owned();
        f.push("loss/");
        let fname = "iter-".to_string() + &ti.to_string() + ".loss";
        f.push(fname);
        f
    }

    pub fn floss_adjmax(&self, ti: usize) -> PathBuf {
        let mut f = self.dout_para().to_owned();
        f.push("loss/");
        let fname = "iter-".to_string() + &ti.to_string() + ".adjmax.loss";
        f.push(fname);
        f
    }

    pub fn floss_initial(&self, file_i: usize, finitial: Option<&Path>) -> PathBuf {
        let mut f = self.dout_para().to_owned();
        f.push("loss/");
        if let Some(finitial) = finitial {
            let initial_for_dir = finitial.file_stem().unwrap().to_str().unwrap();
            f.push(format!("{}/", initial_for_dir));
        }
        let fname = "initial-".to_string() + &file_i.to_string() + ".loss";
        //let fname = "initial.loss";
        f.push(fname);
        f
    }

    pub fn floss_initial_single(&self) -> PathBuf {
        let mut f = self.dout_para().to_owned();
        f.push("loss/");
        let fname = "initial".to_string() + ".single.loss";
        //let fname = "initial.loss";
        f.push(fname);
        f
    }

    pub fn bufwriter_floss(&self, ti: usize) -> BufWriter<File> {
        let floss = self.floss(ti);
        create_dir(&floss.parent().unwrap());
        textfile::bufwriter(&floss)
    }

    pub fn bufwriter_floss_adjmax(&self, ti: usize) -> BufWriter<File> {
        let floss = self.floss_adjmax(ti);
        create_dir(&floss.parent().unwrap());
        textfile::bufwriter(&floss)
    }

    pub fn bufwriter_floss_initial(
        &self,
        file_i: usize,
        finitial: Option<&Path>,
    ) -> BufWriter<File> {
        let floss = self.floss_initial(file_i, finitial);
        create_dir(&floss.parent().unwrap());
        textfile::bufwriter(&floss)
    }

    //pub fn bufwriter_floss_initial(&self, file_i: usize) -> BufWriter<File> {
    //    let floss = self.floss_initial(file_i);
    //    create_dir(&floss.parent().unwrap());
    //    textfile::bufwriter(&floss)
    //}

    pub fn bufwriter_floss_initial_single(&self) -> BufWriter<File> {
        let floss = self.floss_initial_single();
        create_dir(&floss.parent().unwrap());
        textfile::bufwriter(&floss)
    }

    pub fn facc(&self) -> PathBuf {
        let mut f = self.dout_para().to_owned();
        let fname = "monitor.acc".to_string();
        f.push(fname);
        f
    }

    pub fn bufwriter_facc(&self) -> BufWriter<File> {
        let facc = self.facc();
        create_dir(&facc.parent().unwrap());
        textfile::bufwriter(&facc)
    }

    pub fn is_exist_wgt_nonzero(&self) -> bool {
        let fwgt = self.get_file_wgt();
        let is_exist_wgt = textfile::exist_file(&fwgt);
        let is_non_zero = textfile::is_nonzero(&fwgt);

        is_exist_wgt && is_non_zero
    }
}

#[derive(Clone, Debug)]
pub struct DoutScoreFile {
    dout_score: PathBuf,
}

impl DoutScoreFile {
    pub fn new(dout_score: PathBuf) -> Self {
        Self { dout_score }
    }

    pub fn dout_score(&self) -> &Path {
        &self.dout_score
    }

    pub fn dir_score_cv(&self, cvi: usize) -> Self {
        let mut d = self.dout_score().to_owned();
        let dpara = String::from("cv-") + &cvi.to_string();
        d.push(dpara);

        Self::new(d)
    }

    pub fn dout_score_para(&self) -> DoutScoreParaFile {
        DoutScoreParaFile::new(self.dout_score().to_owned())
    }

    pub fn dout_score_para_lr(&self, learning_rate: f64) -> DoutScoreParaFile {
        // merge with DoutFile::dout_para()
        let mut d = self.dout_score().to_owned();
        let mut dpara = String::from("para");
        //if let Some(lr) = learning_rate {
        dpara += &(String::from("_lr-") + &learning_rate.to_string());
        //dpara += &(String::from(".lr") + &learning_rate.to_string());
        //}
        d.push(dpara);

        DoutScoreParaFile::new(d)
    }
}

#[derive(Clone, Debug)]
pub struct DoutScoreParaFile {
    dout_score_para: PathBuf,
}

impl DoutScoreParaFile {
    pub fn new(dout_score_para: PathBuf) -> Self {
        Self { dout_score_para }
    }

    pub fn dout_score_para(&self) -> &Path {
        &self.dout_score_para
    }

    pub fn fname_score_concat_createdir(
        &self,
        nocov: bool,
        is_nsnv: bool,
        integrate: bool,
    ) -> PathBuf {
        let f = if integrate {
            self.fname_score_integrate(nocov)
        } else {
            self.fname_score_concat(nocov, is_nsnv)
        };
        create_dir(&f.parent().unwrap());
        f
    }

    pub fn fname_score_integrate(&self, nocov: bool) -> PathBuf {
        let mut d = self.dout_score_para().to_owned();
        let fname = if nocov {
            "boosting.score"
            //"score.tsv"
        } else {
            "boosting.scorecov"
            //"score.withcov.tsv"
        };
        d.push(fname);
        d
    }

    pub fn fname_score_concat(&self, nocov: bool, is_nsnv: bool) -> PathBuf {
        let para = if is_nsnv {
            "n".to_string()
        } else {
            "iter".to_string()
        };

        let mut d = self.dout_score_para().to_owned();
        let fname = if nocov {
            "boosting_".to_string() + &para + ".score"
        } else {
            "boosting_".to_string() + &para + ".scorecov"
        };
        d.push(fname);
        d
    }
}

pub enum WgtDoutOrFile {
    Dout(DoutFile),
    File(PathBuf),
}

impl WgtDoutOrFile {
    pub fn new(dout: Option<DoutFile>, file: Option<PathBuf>) -> Self {
        if dout.is_some() == file.is_some() {
            panic!("Either dout or file should be Some");
        }

        if let Some(dout) = dout {
            WgtDoutOrFile::Dout(dout)
        } else if let Some(file) = file {
            WgtDoutOrFile::File(file)
        } else {
            panic!("Either dout or file should be Some");
        }
    }

    pub fn new_path(dout: Option<PathBuf>, file: Option<PathBuf>) -> Self {
        let dout = dout.map(|x| DoutFile::new(x));
        Self::new(dout, file)
    }

    pub fn dout_wgt(&self) -> &DoutFile {
        match self {
            WgtDoutOrFile::Dout(dout) => dout,
            WgtDoutOrFile::File(_) => panic!("WgtDoutOrFile::File(PathBuf) not implemented"),
        }
    }

    pub fn is_file(&self) -> bool {
        match self {
            WgtDoutOrFile::Dout(_) => false,
            WgtDoutOrFile::File(_) => true,
        }
    }

    //pub fn get_file_wgt(&self) -> PathBuf {
    //	match self{
    //		WgtDoutOrFile::Dout(dout) => dout.get_file_wgt(),
    //		WgtDoutOrFile::File(f) => f.to_owned(),
    //	}
    //}
}

pub fn create_dir(dout: &Path) {
    fs::create_dir_all(&dout).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_dname_para() {
        let dout = DoutFile::new(PathBuf::from("./abc"));
        let lr = 0.1;
        let dout_para = dout.dout_para_lr(lr);
        assert_eq!(dout_para.dout_para(), &PathBuf::from("./abc/para_lr-0.1/"));
    }

    #[test]
    fn test_get_dname_para_lrnone() {
        let dout = DoutFile::new(PathBuf::from("./abc"));
        let lr = 1.0;
        let dout_para = dout.dout_para_lr(lr);
        assert_eq!(dout_para.dout_para(), &PathBuf::from("./abc/para_lr-1/"));

        //let dout = PathBuf::from("./abc");
        //let lr = 1.0;
        //let dout_para = get_dname_para(&dout, lr);
        //assert_eq!(dout_para, PathBuf::from("./abc/para_lr-1/"));
        //assert_eq!(dout_para, PathBuf::from("./abc/para/"));
    }
}
