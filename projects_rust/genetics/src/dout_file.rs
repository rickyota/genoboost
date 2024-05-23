//use std::fs;
//use std::fs::File;
//use std::io::BufWriter;
use std::path::{Path, PathBuf};

//use crate::textfile;
use crate::wgt;

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

    // TODO: better
    pub fn fscore_from_fwgt(&self, fwgt: &Path) -> PathBuf {
        let fname_score = fwgt.file_stem().unwrap().to_str().unwrap().to_owned() + ".score";
        let fout_score = self.dout_score().join(fname_score);
        fout_score
    }

    pub fn fscore_concat(&self, para: &str, fwgt: &Path) -> PathBuf {
        // fwgt to get method name
        // ex. 'clump_p-0.1_n-100' -> 'clump_p-0.1_n.score'

        let method = fwgt
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .split(&("_".to_string() + para + "-"))
            .collect::<Vec<&str>>()[0]
            .to_string();

        // [method]_n.score
        //let method = fwgt
        //    .file_name()
        //    .unwrap()
        //    .to_str()
        //    .unwrap()
        //    .split("_")
        //    .collect::<Vec<&str>>()[0]
        //    .to_string();
        let fname_score = method + "_" + para + ".score";
        let fout_score = self.dout_score().join(fname_score);
        //let fout_score = dscore.join(fname_score);
        fout_score
    }

    pub fn dir_score_cv(&self, cvi: usize) -> Self {
        let mut d = self.dout_score().to_owned();
        let dpara = String::from("cv-") + &cvi.to_string();
        d.push(dpara);

        Self::new(d)
    }
}

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

    pub fn get_files_wgt(&self) -> Vec<PathBuf> {
        wgt::io::get_files_wgt(&self.dout())
        //let mut files = vec![];
        //for entry in fs::read_dir(self.dout()).unwrap() {
        //    let entry = entry.unwrap();
        //    let path = entry.path();
        //    if path.extension().unwrap() == "wgt" {
        //        files.push(path);
        //    }
        //}
        //files
    }

    //pub fn dout_cv(&self, cvi: usize) -> Self {
    //    let mut d = self.dout().to_owned();
    //    let dpara = String::from("cv-") + &cvi.to_string();
    //    d.push(dpara);
    //    Self::new(d)
    //}

    // when increase para, add them here as argument
    //fn get_dname_para(&self, learning_rate: f64) -> PathBuf {
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

    //pub fn get_fname_wgt_best_para(&self) -> PathBuf {
    //    let mut p = self.dout().to_owned();
    //    p.push("boosting.wgt.para");
    //    p
    //}

    //pub fn get_fname_wgt_best(&self) -> PathBuf {
    //    let mut p = self.dout().to_owned();
    //    p.push("boosting.wgt");
    //    p
    //}

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

    //pub fn create_dout(&self) {
    //    create_dir(self.dout());
    //}
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

pub fn para_from_fwgt(fwgt: &Path, para: &str) -> String {
    let para = fwgt
        .file_stem() // exclude .wgt here
        .unwrap()
        .to_str()
        .unwrap()
        .split(&("_".to_string() + para + "-"))
        .collect::<Vec<&str>>()[1]
        .to_string();
    para
}

pub fn is_fwgt_concat(fwgt: &Path, concat_para: &str) -> bool {
    let is_concat = fwgt
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .contains(&("_".to_string() + concat_para + "-"));

    if is_concat {
        // check if file stem name end with _(para)-*
        if fwgt
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .split(&("_".to_string() + concat_para + "-"))
            .collect::<Vec<&str>>()
            .last()
            .unwrap()
            .contains("_")
        {
            panic!("File stem name should end with _(para)-*. ");
        }
    };
    is_concat
}

//pub fn create_dir(dout: &Path) {
//    fs::create_dir_all(&dout).unwrap();
//}

//#[cfg(test)]
//mod tests {
//    use super::*;
//
//    #[test]
//    fn test_get_dname_para() {
//        let dout = DoutFile::new(PathBuf::from("./abc"));
//        let lr = 0.1;
//        let dout_para = dout.dout_para_lr(lr);
//        assert_eq!(dout_para.dout_para(), &PathBuf::from("./abc/para_lr-0.1/"));
//    }
//
//    #[test]
//    fn test_get_dname_para_lrnone() {
//        let dout = DoutFile::new(PathBuf::from("./abc"));
//        let lr = 1.0;
//        let dout_para = dout.dout_para_lr(lr);
//        assert_eq!(dout_para.dout_para(), &PathBuf::from("./abc/para_lr-1/"));
//
//        //let dout = PathBuf::from("./abc");
//        //let lr = 1.0;
//        //let dout_para = get_dname_para(&dout, lr);
//        //assert_eq!(dout_para, PathBuf::from("./abc/para_lr-1/"));
//        //assert_eq!(dout_para, PathBuf::from("./abc/para/"));
//    }
//}
