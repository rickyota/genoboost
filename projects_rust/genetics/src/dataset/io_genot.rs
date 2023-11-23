//! TODO: make fgenot class
//! since for plink2, columns are not fixed and troublesome to detect them every time.
//! when creating fgenot class, detect columns etc.

pub mod file_genot;
pub mod load;
pub mod load_score;

use crate::sample;
use crate::{textfile, vec, Chrom, SnvId};

pub use file_genot::FileGenot;
use serde::Deserialize;
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::ffi::OsStr;
use std::ffi::OsString;
use std::fs::File;
use std::io::prelude::*; // for seek
use std::io::Read;
use std::io::SeekFrom;
use std::path::{Path, PathBuf};

// 8 x 1 bit
//type B8 = u8;
// 4 x 2 bit
//type B8_2 = u8;

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub enum GenotFormat {
    Plink1,
    Plink2,
    Plink2Vzs,
}

impl GenotFormat {
    // use match; try not to use this
    //fn is_plink2(self) -> bool {
    //    match self {
    //        Self::Plink => false,
    //        Self::Plink2 | Self::Plink2Vzs => true,
    //    }
    //}
    //fn is_plink(self) -> bool {
    //    match self {
    //        Self::Plink => true,
    //        Self::Plink2 | Self::Plink2Vzs => false,
    //    }
    //}
}

// error:
//pub fn judge_split_chrom<P: AsRef<Path>>(fin: P) -> bool {
pub fn judge_split_chrom(fin: &Path) -> bool {
    //fin.contains("%")
    let fname = fin.file_name().unwrap().to_owned().into_string().unwrap();
    fname.contains("%")
}

//pub fn replace_fname(fin: &Path, chrom: Option<&Chrom>) -> String {
pub fn replace_fname(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
    match chrom {
        Some(chrom) => {
            let mut dname = fin.parent().unwrap().to_path_buf();
            let fname = fin.file_name().unwrap().to_owned().into_string().unwrap();
            let fname = fname.replace("%", &chrom.to_string());
            dname.push(fname);
            dname
        }
        None => fin.to_owned(),
    }
    //fin.replace("%", &chrom.to_string())
}

pub fn fname_chrom(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
    replace_fname(fin, chrom)
}

// pathbuf ver. of fout+"_abc"
pub fn add_name(f: &Path, ext: impl AsRef<OsStr>) -> PathBuf {
    let mut fnew: OsString = f.into();
    fnew.push(ext.as_ref());
    fnew.into()
}

// https://internals.rust-lang.org/t/pathbuf-has-set-extension-but-no-add-extension-cannot-cleanly-turn-tar-to-tar-gz/14187/11
pub fn add_ext(f: &Path, ext: impl AsRef<OsStr>) -> PathBuf {
    let mut fnew: OsString = f.into();
    fnew.push(".");
    fnew.push(ext.as_ref());
    fnew.into()
}

pub fn fname_plinks_sample(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> PathBuf {
    let f = replace_fname(fin, chrom);
    let ext = match gfmt {
        GenotFormat::Plink1 => "fam",
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => "psam",
    };
    add_ext(&f, ext)
}

//pub fn fname_sample(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
//    let f = replace_fname(fin, chrom);
//    add_ext(&f, "fam")
//    //replace_fname(fin, chrom) + ".fam"
//}

pub fn fname_plinks_snv(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> PathBuf {
    let f = replace_fname(fin, chrom);
    let ext = match gfmt {
        GenotFormat::Plink1 => "bim",
        GenotFormat::Plink2 => "pvar",
        GenotFormat::Plink2Vzs => "pvar.zst",
    };
    add_ext(&f, ext)
}

pub fn fname_plinks_genot(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> PathBuf {
    let f = replace_fname(fin, chrom);
    let ext = match gfmt {
        GenotFormat::Plink1 => "bed",
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => "pgen",
    };
    add_ext(&f, ext)
}

pub fn bed_per_snv_size(n: usize) -> usize {
    (n + 3) / 4
}

pub fn calculate_bed_size_genotype(m: usize, n: usize) -> usize {
    m * bed_per_snv_size(n)
    //m * ((n + 3) / 4)
}

pub fn calculate_bed_size(m: usize, n: usize) -> usize {
    3 + calculate_bed_size_genotype(m, n)
    //3 + m * ((n + 3) / 4)
}

/// get smallest number of chrom whose plink exists.
/// mainly used for get .fam
/// return if fin is not split chrom or all chrom does not exist
pub fn choose_a_chrom_exist(fin: &Path, gfmt: GenotFormat) -> Option<Chrom> {
    if !judge_split_chrom(fin) {
        None
    } else {
        for chrom_i in Chrom::variants().iter() {
            let fin_fam = fname_plinks_sample(fin, gfmt, Some(chrom_i));
            if textfile::able_open_file(&fin_fam) {
                return Some(chrom_i.clone());
            }
        }
        None
    }
}

pub fn fname_fam_exist_chrom(fin: &Path, gfmt: GenotFormat) -> Option<PathBuf> {
    if !judge_split_chrom(fin) {
        Some(fname_plinks_sample(fin, gfmt, None))
    } else {
        let chrom = choose_a_chrom_exist(fin, gfmt);
        chrom.map(|chrom_| fname_plinks_sample(fin, gfmt, Some(&chrom_)))
    }
}

// return None for missing chrom file
pub fn compute_num_sample_chrom(
    fin: &Path,
    gfmt: GenotFormat,
    chrom: Option<&Chrom>,
) -> Option<usize> {
    let fin_sample = fname_plinks_sample(&fin, gfmt, chrom);
    // replace chrom if any
    //let fin = replace_fname(fin, chrom);
    //let fin_sample = fname_plinks_sample(&fin, gfmt, None);
    let linen = match gfmt {
        GenotFormat::Plink1 => textfile::compute_num_line_text(&fin_sample, None),
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
            textfile::compute_num_line_text(&fin_sample, None)
        }
    };
    let has_header = match gfmt {
        GenotFormat::Plink1 => false,
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => has_bim_header_plink2(&fin_sample, None),
    };

    let linen = linen.ok();

    if has_header {
        linen.map(|x| x - 1)
    } else {
        linen
    }
}

/// assume only part of chrom files could exist
pub fn compute_num_sample(fin: &Path, gfmt: GenotFormat) -> Option<usize> {
    // use just one of valid chrom
    if !judge_split_chrom(fin) {
        compute_num_sample_chrom(fin, gfmt, None)
        //let fin_fam = fname_plinks_sample(fin, gfmt, None);
        //textfile::compute_num_line_text(&fin_fam, None)
    } else {
        let chrom = choose_a_chrom_exist(fin, gfmt);
        match chrom {
            Some(chrom_) => compute_num_sample_chrom(fin, gfmt, Some(&chrom_)),
            None => None,
        }
    }
}

pub fn compute_num_snv_file_chrom(
    fin: &Path,
    gfmt: GenotFormat,
    chrom: Option<&Chrom>,
) -> Option<usize> {
    let fin_snv = fname_plinks_snv(&fin, gfmt, chrom);
    // replace chrom if any
    //let fin = replace_fname(fin, chrom);
    //let fin_snv = fname_plinks_snv(&fin, gfmt, None);
    let linen = match gfmt {
        GenotFormat::Plink1 => textfile::compute_num_line_text(&fin_snv, None),
        GenotFormat::Plink2 => textfile::compute_num_line_text(&fin_snv, None),
        GenotFormat::Plink2Vzs => textfile::compute_num_line_text(&fin_snv, Some("zst")),
    };
    let has_header = match gfmt {
        GenotFormat::Plink1 => false,
        GenotFormat::Plink2 => has_bim_header_plink2(&fin_snv, None),
        GenotFormat::Plink2Vzs => has_bim_header_plink2(&fin_snv, Some("zst")),
    };

    let linen = linen.ok();

    if has_header {
        linen.map(|x| x - 1)
        //Some(linen - 1)
    } else {
        linen
        //Some(linen)
    }
}

/*
// merge above
pub fn compute_num_snv_file_chrom(
    fin: &Path,
    gfmt: GenotFormat,
    chrom: Option<&Chrom>,
) -> Option<usize> {
    // replace chrom if any
    let fin = replace_fname(fin, chrom);
    compute_num_snv_file(&fin, gfmt)
} */

/// assume only part of chrom files could exist
pub fn compute_num_snv(fin: &Path, gfmt: GenotFormat) -> Option<usize> {
    if !judge_split_chrom(fin) {
        compute_num_snv_file_chrom(&fin, gfmt, None)
        //compute_num_snv_file(&fin, gfmt)
    } else {
        let mut num_snv = 0;
        //for chrom_i in 1..=22 {
        for chrom_i in Chrom::variants().iter() {
            let num_snv_chrom = compute_num_snv_file_chrom(fin, gfmt, Some(chrom_i));
            //let fin_bim = get_fname_bim(fin, Some(chrom_i));
            //let num_snv_chrom = match text::compute_num_line(&fin_bim) {
            let num_snv_chrom = match num_snv_chrom {
                Some(v) => v,
                None => 0,
            };
            num_snv += num_snv_chrom;
        }

        if num_snv == 0 {
            return None;
        } else {
            return Some(num_snv);
        }
    }
}

// TODO: shoudld use `check_open_file()`
// here, should be check_exist
/// Check ...
/// If chrom is split,
///     fin_chrom.bim, .fam, .bed exists
/// allow some of chorom not exist
pub fn check_valid_fin(fin: &Path, gfmt: GenotFormat) {
    if !judge_split_chrom(fin) {
        // check three files exist
        textfile::check_exist_file(&fname_plinks_sample(fin, gfmt, None));
        textfile::check_exist_file(&fname_plinks_snv(fin, gfmt, None));
        textfile::check_exist_file(&fname_plinks_genot(fin, gfmt, None));
    } else {
        for chrom_i in Chrom::variants().iter() {
            let mut is_exists = vec![false, false, false];
            is_exists[0] = textfile::exist_file(&fname_plinks_sample(fin, gfmt, Some(chrom_i)));
            is_exists[1] = textfile::exist_file(&fname_plinks_snv(fin, gfmt, Some(chrom_i)));
            is_exists[2] = textfile::exist_file(&fname_plinks_genot(fin, gfmt, Some(chrom_i)));

            // check if all elements are the same
            if !is_exists.iter().all(|v| *v == is_exists[0]) {
                if !is_exists[0] {
                    panic!(
                        "Fam file does not exist: {:?}",
                        fname_plinks_sample(fin, gfmt, Some(chrom_i))
                    );
                    //panic!("Fam file does not exist even though ~ file exists.");
                } else if !is_exists[1] {
                    panic!(
                        "Bim file does not exist: {:?}",
                        fname_plinks_snv(fin, gfmt, Some(chrom_i))
                    );
                } else if !is_exists[2] {
                    panic!(
                        "Bed file does not exist: {:?}",
                        fname_plinks_genot(fin, gfmt, Some(chrom_i))
                    );
                }
            }
        }
    }
}

//  TODO: should unwrap here?
/// return bed_size if valid, error otherwise
//pub fn check_valid_bed(fin: &str, n: usize, m: usize) -> Result<usize, String> {
pub fn check_valid_bed(
    fin: &Path,
    gfmt: GenotFormat,
    chrom: Option<&Chrom>,
    m: usize,
    n: usize,
) -> Result<usize, Box<dyn Error>> {
    let fin_bed = fname_plinks_genot(fin, gfmt, chrom);
    // check if open
    let mut reader = File::open(fin_bed)?;

    // check if size is correct
    let f_end: usize = reader.seek(SeekFrom::End(0)).unwrap() as usize;
    log::debug!("file end {}", f_end);
    let bed_size = calculate_bed_size(m, n);
    if f_end != bed_size {
        return Err(format!(
            "File size of .bed is wrong: {} vs correct {}.",
            f_end, bed_size
        )
        .into());
    }

    // check if the first 3 bytes are correct.
    reader.seek(SeekFrom::Start(0)).unwrap();
    let mut buf: Vec<u8> = Vec::with_capacity(n);
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    //log::debug!("{:?}", buf);
    if buf != vec![0x6cu8, 0x1b, 0x01] {
        return Err("Magic number of .bed file is wrong.".into());
    }
    Ok(bed_size)
}

fn load_fam(fin: &Path, gfmt: GenotFormat) -> Vec<u8> {
    match gfmt {
        GenotFormat::Plink1 => {
            panic!("Do not use --phe-name for plink1.");
        }
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
            let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
            textfile::read_file_to_end(&fin_fam, None).unwrap()
            //&phe_buf_v[..]
            //phe_buf_v
        }
    }
}

// If phe_buf=None, phe_name=None, gfmt!=plink1, then return None
pub fn load_ys_buf_option(
    fin: &Path,
    gfmt: GenotFormat,
    phe_buf: Option<&[u8]>,
    phe_name: Option<&str>,
    sample_id_to_n: &HashMap<String, usize>,
    //use_samples: &[bool],
) -> Option<Vec<bool>> {
    let phe_buf_v: Vec<u8>;
    let valss = if let Some(phe_name) = phe_name {
        let phe_buf = match phe_buf {
            None => {
                // phe is in .psam
                phe_buf_v = load_fam(fin, gfmt);
                &phe_buf_v[..]
                //match gfmt {
                //    GenotFormat::Plink1 => {
                //        panic!("Do not use --phe-name for plink1.");
                //    }
                //    GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
                //        let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
                //        phe_buf_v = textfile::read_file_to_end(&fin_fam, None).unwrap();
                //        &phe_buf_v[..]
                //    }
                //}
            }
            Some(x) => {
                // if phe_name is in phe_buf, use the col
                // elif phe in .psam, use the col
                if textfile::coli_of_header_buf(x, phe_name).is_some() {
                    x
                } else {
                    phe_buf_v = load_fam(fin, gfmt);
                    if textfile::coli_of_header_buf(&phe_buf_v[..], phe_name).is_some() {
                        &phe_buf_v[..]
                    } else {
                        panic!("phe_name is not in fin-phe or psam file.");
                    }
                }
            }
        };
        //TODO
        let col_iid = 0;
        let col_phe = textfile::coli_of_header_buf(phe_buf, phe_name)
            .unwrap_or_else(|| panic!("phe_name is not in fin-phe or psam."));
        textfile::load_table_cols_buf(&phe_buf[..], &[col_iid, col_phe], true)
    } else {
        match gfmt {
            GenotFormat::Plink1 => {
                //TODO
                let col_iid = 1;
                let col_y = 5;
                let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
                textfile::load_table_cols(&fin_fam, &[col_iid, col_y], false)
            }
            GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
                return None;
                //panic!("Use fin_phe for plink2.")
            }
        }
    };

    // TODO: use text::coli_of_header_buf()
    //let col_phe=textfile::coli_of_header_buf(buf, col)

    let vals = sample::vals_align_id(&valss[1], &valss[0], sample_id_to_n);

    log::debug!("vals[0]: {}", vals[0]);

    let uniq = vec::uniq_string(&vals);
    let code_type = if uniq == HashSet::from_iter([String::from("0"), String::from("1")]) {
        "01"
    } else if uniq == HashSet::from_iter([String::from("1"), String::from("2")]) {
        "12"
    } else {
        panic!("Unknown coding {:?}", uniq);
    };

    fn decode_phe(val: &str, code_type: &str) -> bool {
        let code: u8 = (*val).parse::<u8>().unwrap();
        if code_type == "01" {
            match code {
                0 | 1 => code != 0,
                z => panic!("Unknown phenotype included: {}.", z),
            }
        } else if code_type == "12" {
            match code {
                1 | 2 => (code - 1) != 0,
                z => panic!("Unknown phenotype included: {}.", z),
            }
        } else {
            // TODO: match
            panic!("Unknown code_type {:?}", code_type);
        }
    }

    let ys = vals.iter().map(|x| decode_phe(x, code_type)).collect();

    Some(ys)
}

// TODO: what if fin.fam includes phe=9?
// -> let user exclude samples with --fin-sample
pub fn load_ys_buf(
    fin: &Path,
    gfmt: GenotFormat,
    phe_buf: Option<&[u8]>,
    phe_name: Option<&str>,
    sample_id_to_n: &HashMap<String, usize>,
    //use_samples: &[bool],
) -> Vec<bool> {
    //let mut ys: Vec<bool> = Vec::with_capacity(len_n);

    load_ys_buf_option(fin, gfmt, phe_buf, phe_name, sample_id_to_n).unwrap()

    /*     let phe_buf_v: Vec<u8>;
    let valss = if let Some(phe_name) = phe_name {
        let phe_buf = match phe_buf {
            None => {
                // phe is in .psam
                match gfmt {
                    GenotFormat::Plink1 => {
                        panic!("Do not use --phe-name for plink1.");
                    }
                    GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
                        let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
                        phe_buf_v = textfile::read_file_to_end(&fin_fam, None).unwrap();
                        &phe_buf_v[..]
                    }
                }
            }
            Some(x) => x,
        };
        //TODO
        let col_iid = 0;
        let col_phe = textfile::coli_of_header_buf(phe_buf, phe_name)
            .unwrap_or_else(|| panic!("phe_name is not in fin-phe."));
        //let col = vec![phe_name];
        textfile::load_table_cols_buf(&phe_buf[..], &[col_iid, col_phe], true)
    } else {
        match gfmt {
            GenotFormat::Plink1 => {
                //TODO
                let col_iid = 1;
                let col_y = 5;
                let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
                textfile::load_table_cols(&fin_fam, &[col_iid, col_y], false)
            }
            GenotFormat::Plink2 | GenotFormat::Plink2Vzs => panic!("Use fin_phe for plink2."),
        }
    };

    // TODO: use text::coli_of_header_buf()
    //let col_phe=textfile::coli_of_header_buf(buf, col)

    let vals = sample::vals_align_id(&valss[1], &valss[0], sample_id_to_n);

    log::debug!("vals[0]: {}", vals[0]);

    let uniq = vec::uniq_string(&vals);
    let code_type = if uniq == HashSet::from_iter([String::from("0"), String::from("1")]) {
        "01"
    } else if uniq == HashSet::from_iter([String::from("1"), String::from("2")]) {
        "12"
    } else {
        panic!("Unknown coding {:?}", uniq);
    };

    fn decode_phe(val: &str, code_type: &str) -> bool {
        let code: u8 = (*val).parse::<u8>().unwrap();
        if code_type == "01" {
            match code {
                0 | 1 => code != 0,
                z => panic!("Unknown phenotype included: {}.", z),
            }
        } else if code_type == "12" {
            match code {
                1 | 2 => (code - 1) != 0,
                z => panic!("Unknown phenotype included: {}.", z),
            }
        } else {
            // TODO: match
            panic!("Unknown code_type {:?}", code_type);
        }
    }

    let ys = vals.iter().map(|x| decode_phe(x, code_type)).collect();

    ys */
}

/* // TODO: n should be Option or remove
// TODO: what if fin.fam includes phe=9?
pub fn load_ys(
    fin: &Path,
    gfmt: GenotFormat,
    fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    use_samples: &[bool],
) -> Vec<bool> {
    //let mut ys: Vec<bool> = Vec::with_capacity(len_n);

    let vals: Vec<String> = if let Some(fin_phe) = fin_phe {
        let col = vec![phe_name.unwrap()];
        let valss = textfile::load_table_cols_by_headers(fin_phe, &col);
        valss[0].clone()
        //(coli,fin_phe.to_owned())
    } else {
        let col_y = 5;
        match gfmt {
            GenotFormat::Plink => {
                let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
                textfile::load_table_col(&fin_fam, col_y, false)
            }
            GenotFormat::Plink2 | GenotFormat::Plink2Vzs => panic!("Use fin_phe for plink2."),
        }
    };

    log::debug!("vals[0]: {}", vals[0]);

    let uniq = vec::uniq_string(&vals);

    let code_type = if uniq == HashSet::from_iter([String::from("0"), String::from("1")]) {
        "01"
    } else if uniq == HashSet::from_iter([String::from("1"), String::from("2")]) {
        "12"
    } else {
        panic!("Unknown coding {:?}", uniq);
    };

    fn decode_phe(val: &str, code_type: &str) -> bool {
        let code: u8 = (*val).parse::<u8>().unwrap();
        if code_type == "01" {
            match code {
                0 | 1 => code != 0,
                z => panic!("Unknown phenotype included: {}.", z),
            }
        } else if code_type == "12" {
            match code {
                1 | 2 => (code - 1) != 0,
                z => panic!("Unknown phenotype included: {}.", z),
            }
        } else {
            // TODO: match
            panic!("Unknown code_type {:?}", code_type);
        }
    }

    let mut ys: Vec<bool> = Vec::new();
    for (val, use_i) in vals.iter().zip(use_samples) {
        //log::debug!("{}", val);
        if *use_i {
            let y = decode_phe(val, code_type);
            ys.push(y);

            //n_i += 1;
        }
    }

    //assert_eq!(ys.len(), n);

    ys
} */

// option for missing chrom
fn load_snvs_text_plink(
    fin: &Path,
    gfmt: GenotFormat,
    chrom: Option<&Chrom>,
) -> Option<Vec<SnvId>> {
    let mut snvs: Vec<SnvId> = vec![];
    //let mut snvs: Vec<SnvId> = Vec::with_capacity(m);

    // original order is chrom, rs, None, pos, A1, A2
    // rs, chrom, pos, A1, A2
    let cols = [1usize, 0, 3, 4, 5];

    let fin_bim = fname_plinks_snv(fin, gfmt, chrom);
    //let fin_bim = fname_plinks_snv(fin, gfmt, None);
    let buf = textfile::read_file_to_end(&fin_bim, None);
    if buf.is_err() {
        return None;
    }
    let buf = buf.unwrap();

    let vss = textfile::load_table_cols_buf(&buf[..], &cols, false);
    //let vss = textfile::load_table_cols(&fin_bim, &cols, false);
    for vi in 0..vss[0].len() {
        snvs.push(SnvId::construct_snv_index(
            vss[0][vi].clone(),
            &vss[1][vi],
            &vss[2][vi],
            vss[3][vi].clone(),
            vss[4][vi].clone(),
        ));
    }
    Some(snvs)
}

#[derive(Clone, Debug, Deserialize)]
struct SnvPlink2In {
    #[serde(alias = "ID", alias = "#ID")]
    id: String,
    #[serde(alias = "CHROM", alias = "#CHROM")]
    chrom: String,
    #[serde(alias = "POS", alias = "#POS")]
    pos: String,
    #[serde(alias = "REF")]
    refa: String,
    #[serde(alias = "ALT")]
    alta: String,
}

impl SnvPlink2In {
    pub fn new(id: String, chrom: String, pos: String, refa: String, alta: String) -> Self {
        if textfile::isin_string(",", &alta) {
            panic!("Cannot have multiple alt alleles.");
        }
        SnvPlink2In {
            id,
            chrom,
            pos,
            refa,
            alta,
        }
    }
}

// this will fail for a file without header but contain string "ID" as name, but should be fine
pub fn has_bim_header_plink2(fin_bim: &Path, compress: Option<&str>) -> bool {
    let buf = textfile::read_file_to_end(fin_bim, compress).unwrap();
    textfile::isin_header("ID", &buf[..])
}

fn load_snvs_tsv(
    fin: &Path,
    gfmt: GenotFormat,
    chrom: Option<&Chrom>,
    compress: Option<&str>,
) -> Option<Vec<SnvId>> {
    let fin_bim = fname_plinks_snv(fin, gfmt, chrom);
    //println!("fin_bim {:?}", fin_bim);

    let mut snvs_in: Vec<SnvPlink2In> = vec![];
    let buf = textfile::read_file_to_end(&fin_bim, compress);
    if buf.is_err() {
        return None;
    }
    let buf = buf.unwrap();

    if has_bim_header_plink2(&fin_bim, compress) {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(&buf[..]);

        for result in rdr.deserialize() {
            let record: SnvPlink2In =
                result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_bim));
            snvs_in.push(record)
        }
    } else {
        //.has_headers(false)
        // then should have exactly same cols as SnvIn, so avoid using Struct
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(&buf[..]);

        // original order is chrom, rs, None, pos, A1, A2
        // rs, chrom, pos, A1, A2
        let cols = [1usize, 0, 3, 4, 5];

        for result in rdr.records() {
            let record = result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_bim));
            println!("{:?}", record);
            snvs_in.push(SnvPlink2In::new(
                record[cols[0]].to_string(),
                record[cols[1]].to_string(),
                record[cols[2]].to_string(),
                record[cols[3]].to_string(),
                record[cols[4]].to_string(),
            ))
        }
    }

    let mut snvs: Vec<SnvId> = vec![];

    for record in snvs_in {
        let snv = SnvId::construct_snv_index(
            record.id.clone(),
            &record.chrom,
            &record.pos,
            record.alta.clone(),
            record.refa.clone(),
        );
        snvs.push(snv);
    }
    Some(snvs)
}

// option for missing chrom
fn load_snvs_chrom(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> Option<Vec<SnvId>> {
    match gfmt {
        GenotFormat::Plink1 => load_snvs_text_plink(fin, gfmt, chrom),
        GenotFormat::Plink2 => load_snvs_tsv(fin, gfmt, chrom, None),
        GenotFormat::Plink2Vzs => load_snvs_tsv(fin, gfmt, chrom, Some("zst")),
    }
}

// TODO: why not use_snvs here?
// since match snv_index is not simple?
pub fn load_snvs(fin: &Path, gfmt: GenotFormat) -> Vec<SnvId> {
    if !judge_split_chrom(fin) {
        //load_snvs_chrom(fin, gfmt, None).unwrap()

        load_snvs_chrom(fin, gfmt, None)
            .unwrap_or_else(|| panic!("Could not load snvs from fin: {:?}", fin))
    } else {
        let mut snvs: Vec<SnvId> = vec![];
        for chrom_i in Chrom::variants().iter() {
            let v = load_snvs_chrom(fin, gfmt, Some(chrom_i));
            if let Some(mut snvs_i) = v {
                snvs.append(&mut snvs_i);
            }
            // else continue
        }
        // all chrom returned None
        if snvs.len() == 0 {
            panic!("No snvs were loaded from {:?}", fin);
        }
        snvs
    }
}

#[derive(Clone, Debug, Deserialize)]
struct SamplePlink2In {
    // ng
    //#[serde(rename = ("FID","FID"))]
    //#[serde(alias = "FID", alias = "#FID")]
    //fid: Option<String>,
    #[serde(alias = "ID", alias = "IID", alias = "#ID", alias = "#IID")]
    iid: String,
}

impl SamplePlink2In {
    pub fn new(iid: String) -> Self {
        SamplePlink2In {
            //fid: Some(fid),
            iid,
        }
    }
}

// this will fail for a file without header but contain string "ID" as name, but should be fine
pub fn has_fam_header_plink2(fin_fam: &Path) -> bool {
    let buf = textfile::read_file_to_end(fin_fam, None).unwrap();
    textfile::isin_header("ID", &buf[..])
}

fn load_samples_id_tsv(fin: &Path, gfmt: GenotFormat, use_samples: Option<&[bool]>) -> Vec<String> {
    let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();
    let mut samples_in: Vec<SamplePlink2In> = vec![];
    if has_fam_header_plink2(&fin_fam) {
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&fin_fam)
            .unwrap();
        // cannot load properly
        //.quote(b'#')

        for result in rdr.deserialize() {
            let record: SamplePlink2In =
                result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_fam));
            samples_in.push(record)
        }
    } else {
        //.has_headers(false)
        // then should have exactly same cols as SnvIn, so avoid using Struct
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(&fin_fam)
            .unwrap();

        for result in rdr.records() {
            let record = result.unwrap_or_else(|_| panic!("Error while reading: {:?}", fin_fam));
            //println!("{:?}", record);
            // the first line sould be FID or IID and both are fine
            samples_in.push(SamplePlink2In::new(record[0].to_string()))
        }
    }

    //let n = use_samples.iter().filter(|&v| *v).count();
    //let mut samples_id: Vec<(String, String)> = Vec::with_capacity(n);

    let mut samples_id: Vec<String> = vec![];

    if use_samples.is_some() {
        for (record, b) in samples_in.into_iter().zip(use_samples.unwrap().iter()) {
            if *b {
                samples_id.push(record.iid);
            }
        }
    } else {
        for record in samples_in.into_iter() {
            samples_id.push(record.iid);
        }
    }

    samples_id
}

fn load_samples_id_text_plink(
    fin: &Path,
    gfmt: GenotFormat,
    use_samples: Option<&[bool]>,
) -> Vec<String> {
    let fin_fam = fname_fam_exist_chrom(fin, gfmt).unwrap();

    // fid, iid
    //let cols = [0usize, 1];
    // assume fid=iid
    let col_iid = 1;
    let cols = [col_iid];
    let samples_in: Vec<Vec<String>> = textfile::load_table_cols(&fin_fam, &cols, false);
    if use_samples.is_some() {
        assert_eq!(samples_in[0].len(), use_samples.unwrap().len());
    }

    //let n = use_samples.iter().filter(|&v| *v).count();
    //let mut sample_id_to_n: HashMap<String, usize> = HashMap::with_capacity(n);
    //let mut samples_id: Vec<(String, String)> = Vec::with_capacity(n);

    let mut samples_id: Vec<String> = vec![];

    if use_samples.is_some() {
        for (n_in_i, v) in use_samples.unwrap().iter().enumerate() {
            if *v {
                samples_id.push(samples_in[0][n_in_i].clone());
            }
        }
    } else {
        for n_in_i in 0..samples_in[0].len() {
            samples_id.push(samples_in[0][n_in_i].clone());
        }
    }

    samples_id
}

pub fn load_samples_id(fin: &Path, gfmt: GenotFormat, use_samples: Option<&[bool]>) -> Vec<String> {
    match gfmt {
        GenotFormat::Plink1 => load_samples_id_text_plink(fin, gfmt, use_samples),
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => load_samples_id_tsv(fin, gfmt, use_samples),
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_get_bed_size() {
        let m: usize = 3;
        let n: usize = 5;
        let bed_size: usize = calculate_bed_size(m, n);
        assert_eq!(bed_size, 9);

        let m: usize = 3;
        let n: usize = 4;
        let bed_size: usize = calculate_bed_size(m, n);
        assert_eq!(bed_size, 6);
    }

    #[test]
    fn test_check_valid_bed() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;

        let m: usize = compute_num_snv(&fin, gfmt).unwrap();
        let n: usize = compute_num_sample(&fin, gfmt).unwrap();
        log::debug!("m,n: {},{}", m, n);
        let bed_size = check_valid_bed(&fin, gfmt, None, m, n).unwrap();
        assert_eq!(bed_size, calculate_bed_size(m, n));
    }

    #[test]
    #[should_panic]
    fn test_check_valid_bed_panic_nofile() {
        // file does not exist
        let fin = PathBuf::from("./toy");
        let gfmt = GenotFormat::Plink1;

        let m: usize = compute_num_snv(&fin, gfmt).unwrap();
        let n: usize = compute_num_sample(&fin, gfmt).unwrap();
        check_valid_bed(&fin, gfmt, None, m, n).unwrap();

        // TODO; check violated first 3 bytes
    }

    #[test]
    #[should_panic]
    fn test_check_valid_bed_panic() {
        //pub fn check_valid_bed(fin_bed: &str, n: usize, m: usize) -> Result<usize, Box<dyn Error>> {

        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;

        let m: usize = compute_num_snv(&fin, gfmt).unwrap();
        let n: usize = compute_num_sample(&fin, gfmt).unwrap();
        // size is wrong
        check_valid_bed(&fin, gfmt, None, m, n - 1).unwrap();
    }

    #[test]
    fn test_load_samples_id_plink() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;

        let v = load_samples_id(&fin, gfmt, None);
        assert_eq!(v[0], "1");
    }

    #[test]
    fn test_load_samples_id_plink2() {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2Vzs;

        let v = load_samples_id(&fin, gfmt, None);
        assert_eq!(v[0], "1");
    }

    /*     #[test]
    fn test_load_ys() {
        let fin = String::from("../../test/data/toy1/genot");

        let n = compute_num_sample(&fin).unwrap();
        let (_, use_samples) = make_use_samples(None, &fin, n);

        let ys = load_ys_b8(&fin, n, &use_samples);
        let ys_expect = vec![true, true, false, false, false];
        //let ys_expect = vec![false, true, false, false, true]; -> this is sex

        for (i, exp) in ys_expect.iter().enumerate() {
            assert_eq!(bit::bget(&ys, i), *exp);
        }
    } */

    #[test]
    fn test_compute_num_sample() {
        let fin = PathBuf::from("../../test/data/toy3/genot");

        for gfmt in [
            GenotFormat::Plink1,
            GenotFormat::Plink2,
            GenotFormat::Plink2Vzs,
        ] {
            println!("gfmt {:?}", gfmt);
            let m = compute_num_sample(&fin, gfmt).unwrap();
            assert_eq!(m, 10)
        }
    }

    #[test]
    fn test_compute_num_snv() {
        let fin = PathBuf::from("../../test/data/toy3/genot");

        for gfmt in [
            GenotFormat::Plink1,
            GenotFormat::Plink2,
            GenotFormat::Plink2Vzs,
        ] {
            println!("gfmt {:?}", gfmt);
            let m = compute_num_snv(&fin, gfmt).unwrap();
            assert_eq!(m, 3)
        }
    }
}
