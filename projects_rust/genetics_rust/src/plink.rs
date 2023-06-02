pub mod load;
pub mod load_score;

use super::{samples, text, Chrom, SnvId};
use crate::B8_2;
use std::collections::HashMap;
use std::error::Error;
use std::ffi::OsStr;
use std::ffi::OsString;
use std::fs::File;
use std::io::prelude::*; // for seek
use std::io::Read;
use std::io::SeekFrom;
use std::path::{Path, PathBuf};

//use crate::alloc; //same
//mod alloc; // error
//use super::snv::Snv; // unnecessary

// 8 x 1 bit
//type B8 = u8;
// 4 x 2 bit
//type B8_2 = u8;

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

pub fn fname_bim(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
    let f = replace_fname(fin, chrom);
    add_ext(&f, "bim")
    //f.set_extension("bim");
    //f
    //replace_fname(fin, chrom) + ".bim"
}

pub fn fname_fam(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
    let f = replace_fname(fin, chrom);
    add_ext(&f, "fam")
    //replace_fname(fin, chrom) + ".fam"
}

pub fn fname_bed(fin: &Path, chrom: Option<&Chrom>) -> PathBuf {
    let f = replace_fname(fin, chrom);
    add_ext(&f, "bed")
    //replace_fname(fin, chrom) + ".bed"
}
/*
// TODO: change chrom into Option<Chrom>,
pub fn replace_fname(fin: &str, chrom: Option<usize>) -> String {
    match chrom {
        Some(chrom) => fin.replace("%", &chrom.to_string()),
        None => fin.to_owned(),
    }
    //fin.replace("%", &chrom.to_string())
}

pub fn get_fname_chrom(fin: &str, chrom: Option<usize>) -> String {
    replace_fname(fin, chrom)
}

pub fn get_fname_bim(fin: &str, chrom: Option<usize>) -> String {
    replace_fname(fin, chrom) + ".bim"
}

pub fn get_fname_fam(fin: &str, chrom: Option<usize>) -> String {
    replace_fname(fin, chrom) + ".fam"
}

pub fn get_fname_bed(fin: &str, chrom: Option<usize>) -> String {
    replace_fname(fin, chrom) + ".bed"
}
*/

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
pub fn choose_a_chrom_exist(fin: &Path) -> Option<Chrom> {
    if !judge_split_chrom(fin) {
        None
    } else {
        for chrom_i in Chrom::variants().iter() {
            let fin_fam = fname_fam(fin, Some(chrom_i));
            if text::able_open_file(&fin_fam) {
                return Some(chrom_i.clone());
            }
        }
        None
    }
}
/*
pub fn get_a_chrom_exist(fin: &str) -> Option<usize> {
    if !judge_split_chrom(fin) {
        None
    } else {
        for chrom_i in 1..=22 {
            let fin_fam = get_fname_fam(fin, Some(chrom_i));
            if text::able_open_file(&fin_fam) {
                return Some(chrom_i);
            }
        }
        None
    }
}
*/

pub fn fname_fam_exist_chrom(fin: &Path) -> Option<PathBuf> {
    if !judge_split_chrom(fin) {
        Some(fname_fam(fin, None))
        //Some(replace_fname(fin, None) + ".fam")
    } else {
        let chrom = choose_a_chrom_exist(fin);
        chrom.map(|chrom_| fname_fam(fin, Some(&chrom_)))
        //chrom.map(|chrom_| replace_fname(fin, Some(&chrom_)) + ".fam")
        /*
        match chrom {
            Some(chrom_) => Some(replace_fname(fin, Some(&chrom_)) + ".fam"),
            None => None,
        }
        */
    }
}

/// assume only part of chrom files could exist
pub fn compute_num_sample(fin: &Path) -> Option<usize> {
    // use just one of valid chrom
    if !judge_split_chrom(fin) {
        let fin_fam = fname_fam(fin, None);
        text::compute_num_line(&fin_fam)
    } else {
        let chrom = choose_a_chrom_exist(fin);
        match chrom {
            Some(chrom_) => {
                let fin_fam = fname_fam(fin, Some(&chrom_));
                // This should not be None
                return text::compute_num_line(&fin_fam);
            }
            None => {
                return None;
            }
        }
        /*
        for chrom_i in 1..=22 {
            let fin_fam = get_fname_fam(fin, Some(chrom_i));
            if let Some(v) = text::compute_num_line(&fin_fam) {
                return Some(v);
            }
            //     return text::compute_num_line(&fin_fam).unwrap();
        }
        */
        //return None;
        //panic!("No fam file of valid chrom exists.");
    }
    /*
    let fin_fam = get_fname_fam(fin, Some(1));
    text::compute_num_line(&fin_fam)
    */
}

pub fn compute_num_snv_chrom(fin: &Path, chrom: Option<&Chrom>) -> Option<usize> {
    let fin_bim = fname_bim(fin, chrom);

    text::compute_num_line(&fin_bim)
}

/// assume only part of chrom files could exist
pub fn compute_num_snv(fin: &Path) -> Option<usize> {
    if !judge_split_chrom(fin) {
        let fin_bim = fname_bim(fin, None);
        text::compute_num_line(&fin_bim)
    } else {
        let mut num_snv = 0;
        //for chrom_i in 1..=22 {
        for chrom_i in Chrom::variants().iter() {
            let num_snv_chrom = compute_num_snv_chrom(fin, Some(chrom_i));
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

/*
pub fn compute_num_bim(fin_bim: &str) -> usize {
    text::compute_num_line(&fin_bim)
}

pub fn compute_num_fam(fin_fam: &str) -> usize {
    text::compute_num_line(&fin_fam)
}
*/

// TODO: shoudld use `check_open_file()`
// here, should be check_exist
/// Check ...
/// If chrom is split,
///     fin_chrom.bim, .fam, .bed exists
/// allow some of chorom not exist
pub fn check_valid_fin(fin: &Path) {
    if !judge_split_chrom(fin) {
        // check three files exist
        text::check_exist_file(&fname_fam(fin, None));
        text::check_exist_file(&fname_bim(fin, None));
        text::check_exist_file(&fname_bed(fin, None));
    } else {
        for chrom_i in Chrom::variants().iter() {
            let mut is_exists = vec![false, false, false];
            is_exists[0] = text::exist_file(&fname_fam(fin, Some(chrom_i)));
            is_exists[1] = text::exist_file(&fname_bim(fin, Some(chrom_i)));
            is_exists[2] = text::exist_file(&fname_bed(fin, Some(chrom_i)));

            // check if all elements are the same
            if !is_exists.iter().all(|v| *v == is_exists[0]) {
                if !is_exists[0] {
                    panic!(
                        "Fam file does not exist: {:?}",
                        fname_fam(fin, Some(chrom_i))
                    );
                    //panic!("Fam file does not exist even though ~ file exists.");
                } else if !is_exists[1] {
                    panic!(
                        "Bim file does not exist: {:?}",
                        fname_bim(fin, Some(chrom_i))
                    );
                } else if !is_exists[2] {
                    panic!(
                        "Bed file does not exist: {:?}",
                        fname_bed(fin, Some(chrom_i))
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
    chrom: Option<&Chrom>,
    m: usize,
    n: usize,
) -> Result<usize, Box<dyn Error>> {
    let fin_bed = fname_bed(fin, chrom);
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

// {00: 2, 01: 3, 10: 1, 11: 0}
const CODE_TO_COUNT_AR: [u8; 4] = [2, 3, 1, 0];

// plink BED code to minor allele counts
// {00: 2, 01: 3, 10: 1, 11: 0}
#[inline]
fn code_to_count(v: u8) -> u8 {
    CODE_TO_COUNT_AR[v as usize]
}

#[inline]
fn byte_to_ped_code(v: B8_2, i: usize) -> u8 {
    (v & (0x03 << (i << 1))) >> (i << 1)
}

// extract ith code (= 2i~2i+1 bits from lower)
// ex. byte_to_count(0b11011000, 2) = code_to_count(0b10) = 1
#[inline]
fn byte_to_count(v: B8_2, i: usize) -> u8 {
    //#define count_pl(c, i) (code2count((((c) & (0x03 << (i << 1))) >> (i << 1))))
    code_to_count(byte_to_ped_code(v, i))
    //code_to_count((v & (0x03 << (i << 1))) >> (i << 1))
}

#[inline]
pub fn buf_to_ped_code(buf: &[B8_2], ni: usize) -> u8 {
    //buf[ni // 4], ni % 3
    byte_to_ped_code(buf[ni >> 2], ni & 3)
}

#[inline]
pub fn buf_to_count(buf: &[B8_2], ni: usize) -> u8 {
    //buf[ni // 4], ni % 3
    byte_to_count(buf[ni >> 2], ni & 3)
}

/* // TODO: what if fin.fam includes phe=9?
pub fn load_ys_b8(fin: &str, n: usize, use_samples: &[bool]) -> Vec<B8> {
    let len_n = n / 8 + 5;

    let mut ys: Vec<B8> = alloc::with_capacity_align_u8(len_n);
    for _ in 0..len_n {
        ys.push(0x00);
    }

    let col_y = 5;
    let fin_fam = fname_fam_exist_chrom(fin).unwrap();
    //let fin_fam = get_fname_fam(fin, Some(1));
    let vals: Vec<String> = text::load_table_col(&fin_fam, col_y, false).unwrap();
    log::debug!("vals: {}", vals[0]);

    let mut n_i = 0;
    for (val, use_i) in vals.iter().zip(use_samples) {
        //log::debug!("{}", val);
        if *use_i {
            let y: i8 = (*val).parse::<i8>().unwrap();
            //let phe: u8 = (*val).parse::<u8>().unwrap() - 1;
            // this seems to take mach time?
            match y {
                1 | 2 => bit::bset(&mut ys, n_i, (y - 1) != 0),
                z => panic!("Unknown phenotype included: {}.", z),
            }

            n_i += 1;
        }
    }

    assert_eq!(n_i, n);

    ys
}
 */

// TODO: n should be Option or remove 
// TODO: what if fin.fam includes phe=9?
pub fn load_ys(fin: &Path, 
        fin_phe: Option<&Path>,
        phe_name: Option<&str>,
    n: usize, use_samples: &[bool]) -> Vec<bool> {
    let len_n = n;

    let mut ys: Vec<bool> = Vec::with_capacity(len_n);


    //let vals: Vec<String> = text::load_table_col(&fin_fam, col_y, false).unwrap();

    let vals: Vec<String>=if let Some(fin_phe)=fin_phe{
        let col=vec![phe_name.unwrap()];
        let valss=text::load_table_cols_by_header(fin_phe,&col).unwrap();
        valss[0].clone()
        //(coli,fin_phe.to_owned())
    }else{
        let col_y = 5;
        let fin_fam = fname_fam_exist_chrom(fin).unwrap();
        //(5,fname_fam_exist_chrom(fin).unwrap())
        text::load_table_col(&fin_fam, col_y, false).unwrap()

    };



    log::debug!("vals[0]: {}", vals[0]);

    //let mut n_i = 0;
    for (val, use_i) in vals.iter().zip(use_samples) {
        //log::debug!("{}", val);
        if *use_i {
            let phe: u8 = (*val).parse::<u8>().unwrap();
            //let phe: u8 = (*val).parse::<u8>().unwrap() - 1;
            // this seems to take mach time?
            match phe {
                1 | 2 => ys.push((phe - 1) != 0),
                //1 | 2 => operate::bset(&mut ys, n_i, (phe - 1) != 0),
                z => panic!("Unknown phenotype included: {}.", z),
            }

            //n_i += 1;
        }
    }

    assert_eq!(ys.len(), n);

    ys
}

// TODO: get m inside
pub fn load_snvs(fin: &Path, m: usize) -> Vec<SnvId> {
    let mut snvs: Vec<SnvId> = Vec::with_capacity(m);

    // original order is chrom, rs, None, pos, A1, A2
    // rs, chrom, pos, A1, A2
    let cols = [1usize, 0, 3, 4, 5];

    if !judge_split_chrom(fin) {
        let fin_bim = fname_bim(fin, None);
        let vss: Vec<Vec<String>> = text::load_table_cols(&fin_bim, &cols, false).unwrap();

        for vi in 0..vss[0].len() {
            snvs.push(SnvId::construct_snv_index(
                vss[0][vi].clone(),
                &vss[1][vi],
                &vss[2][vi],
                vss[3][vi].clone(),
                vss[4][vi].clone(),
            ));
        }
        snvs
    } else {
        // iter through Chrom
        //for chrom_i in 1..=22 {
        for chrom_i in Chrom::variants().iter() {
            let fin_bim = fname_bim(fin, Some(chrom_i));
            //let vss: Vec<Vec<String>> = text::load_table_cols(&fin_bim, &cols, false).unwrap();
            if let Some(vss) = text::load_table_cols(&fin_bim, &cols, false) {
                for vi in 0..vss[0].len() {
                    snvs.push(SnvId::construct_snv_index(
                        vss[0][vi].clone(),
                        &vss[1][vi],
                        &vss[2][vi],
                        vss[3][vi].clone(),
                        vss[4][vi].clone(),
                    ));
                }
            }
            // else continue
        }
        snvs
    }
}

pub fn load_snvs_use(fin_snv: &Path) -> Vec<SnvId> {
    text::check_open_file(fin_snv);

    let m_use = text::compute_num_line(fin_snv).unwrap();

    let mut snvs: Vec<SnvId> = Vec::with_capacity(m_use);

    // TODO: for .snvs with header
    // original order is chrom, rs, None, pos, A1, A2
    // rs, chrom, pos, A1, A2
    //let cols = [1usize, 0, 3, 4, 5];
    let cols = [0usize];
    let vss: Vec<Vec<String>> = text::load_table_cols(&fin_snv, &cols, false).unwrap();

    for vi in 0..vss[0].len() {
        snvs.push(SnvId::construct_snv_index_rs(vss[0][vi].clone()));
    }
    snvs
}

/* pub fn load_snvs_use(fin_snv: &Path) -> Vec<SnvId> { text::check_open_file(fin_snv);

    let m_use = text::compute_num_line(fin_snv).unwrap();

    let mut snvs: Vec<SnvId> = Vec::with_capacity(m_use);

    // original order is chrom, rs, None, pos, A1, A2
    // rs, chrom, pos, A1, A2
    let cols = [1usize, 0, 3, 4, 5];
    let vss: Vec<Vec<String>> = text::load_table_cols(&fin_snv, &cols, false).unwrap();

    for vi in 0..vss[0].len() {
        snvs.push(SnvId::construct_snv_index_string(
            vss[0][vi].clone(),
            &vss[1][vi],
            &vss[2][vi],
            vss[3][vi].clone(),
            vss[4][vi].clone(),
        ));
    }
    snvs
} */

/// Use rs to extract snvs.
/// Format of fin_snv is the same as plink `--extract`.
/// FIXME: raise error? or do sth if SNVs are duplicated in fin_snv
pub fn make_use_snvs(fin_snv: Option<&Path>, snvs_in: &Vec<SnvId>) -> (usize, Vec<bool>) {
    let m_in: usize = snvs_in.len();

    if fin_snv.is_none() {
        let use_snvs = vec![true; m_in];
        return (m_in, use_snvs);
    }

    text::check_open_file(fin_snv.unwrap());

    let mut use_snvs = vec![false; m_in];

    // map: sid_in -> plink index
    let mut rs_in_to_index = HashMap::with_capacity(m_in);

    for (si, s) in snvs_in.iter().enumerate() {
        rs_in_to_index.insert(s.rs(), si);
        //rs_in_to_index.insert(s.sida(), si);
    }

    let snvs_use = load_snvs_use(fin_snv.unwrap());
    //let m_use = snvs_use.len();

    let mut m: usize = 0;
    for snv in snvs_use.iter() {
        //let sida: String = snvs_use[mui].get_sida();
        //let sida: String = snv.get_sida().to_owned();
        //let sida = snv.sida();
        let rs = snv.rs();

        match rs_in_to_index.get(rs) {
            Some(v) => {
                use_snvs[*v] = true;
                m += 1;
            }
            // ignore unfound SNVs
            None => {
                log::info!("SNV in fin_snv was not found in plink: {}.", rs);
            }
        }

        //if let Some(v) = rs_in_to_index.get(rs) {
        //    use_snvs[*v] = true;
        //    m += 1;
        //}
        // if you want to panic
        //use_snvs[sida_to_index[&sida]] = true;
    }
    (m, use_snvs)
}

/// return vec of tuple of (fid, iid)
/// TODO: or just return Vec<Vec<String>>?
pub fn load_samples_use(fin_sample: &Path) -> Vec<(String, String)> {
    text::check_open_file(fin_sample);

    let n_use = text::compute_num_line(fin_sample).unwrap();

    let mut samples: Vec<(String, String)> = Vec::with_capacity(n_use);
    //let mut snvs: Vec<Snv> = Vec::with_capacity(m_use);

    //let fin_bim = get_fname_bim(fin_snv, Some(1));

    // load fid and iid
    let cols = [0usize, 1];
    let vss: Vec<Vec<String>> = text::load_table_cols(&fin_sample, &cols, false).unwrap();

    for vi in 0..vss[0].len() {
        samples.push((vss[0][vi].clone(), vss[1][vi].clone()));
    }
    samples
}

pub fn load_samples_id(fin: &Path, use_samples: &[bool]) -> Vec<(String, String)> {
    let fin_fam = fname_fam_exist_chrom(fin).unwrap();
    //let fin_fam = plink::get_fname_fam(fin, Some(1));
    // fid, iid
    let cols = [0usize, 1];
    let samples_in: Vec<Vec<String>> = text::load_table_cols(&fin_fam, &cols, false).unwrap();

    assert_eq!(samples_in[0].len(), use_samples.len());

    /*
    let mut samples: Vec<(String, String)> = Vec::with_capacity(n_use);
    // load fid and iid
    for vi in 0..samples_in[0].len() {
        samples.push((samples_in[0][vi].clone(), samples_in[1][vi].clone()));
    }
    */

    let n = use_samples.iter().filter(|&v| *v).count();

    //let mut sample_id_to_n: HashMap<String, usize> = HashMap::with_capacity(n);
    let mut samples_id: Vec<(String, String)> = Vec::with_capacity(n);

    for (n_in_i, v) in use_samples.iter().enumerate() {
        if *v {
            samples_id.push((samples_in[0][n_in_i].clone(), samples_in[1][n_in_i].clone()));
        }
    }

    samples_id
}

pub fn make_use_samples(fin_sample: Option<&Path>, fin: &Path, n_in: usize) -> (usize, Vec<bool>) {
    if fin_sample.is_none() {
        //for _ in 0..n_in {
        //    use_samples.push(true);
        //}
        let use_samples = vec![true; n_in];
        return (n_in, use_samples);
    }

    text::check_open_file(fin_sample.unwrap());

    let mut use_samples = vec![false; n_in];
    // initialize with false
    //for _ in 0..n_in {
    //    use_samples.push(false);
    //}

    let fin_fam = fname_fam_exist_chrom(fin).unwrap();
    //let fin_fam = get_fname_fam(fin, Some(1));
    // fid, iid
    let cols = [0usize, 1];
    let samples_in: Vec<Vec<String>> = text::load_table_cols(&fin_fam, &cols, false).unwrap();
    //log::debug!("vals: {}", samples_in[0][0]);

    /*
    #[inline]
    fn sample_unique(fid: String, iid: &str) -> String {
        fid + ":" + iid
        //fid.clone() + ":" + iid
    }
    */

    // map: sample_in -> index
    let mut sample_in_to_index: HashMap<String, usize> = HashMap::with_capacity(n_in);
    for si in 0..samples_in[0].len() {
        sample_in_to_index.insert(
            samples::sample_id(samples_in[0][si].clone(), &samples_in[1][si]),
            si,
        );
        //sample_in_to_index.insert(samples_in[0][si].clone() + ":" + &samples_in[1][si], si);
    }

    let samples_use: Vec<(String, String)> = load_samples_use(fin_sample.unwrap());
    //let n_use = samples_use.len();

    let mut n: usize = 0;
    //for nui in 0..n_use {
    for sample in samples_use.into_iter() {
        let sample_id = samples::sample_id(sample.0, &sample.1);
        //let sample_id: String = sample.0.clone() + ":" + &sample.1;
        if let Some(v) = sample_in_to_index.get(&sample_id) {
            use_samples[*v] = true;
            n += 1;
        }
    }

    return (n, use_samples);
}

/// x: vector of minor allele (0,1,2,3)
/// convert 3 (missing) to mode
pub fn missing_to_mode(x: &mut [u8]) {
    // count 0,1,2,3
    let mut counts_allele = vec![0usize; 4];
    //let mut counts_allele = Vec::with_capacity(4);
    //for _ in 0..=3 {
    //    counts_allele.push(0);
    //}

    for x_v in x.iter() {
        counts_allele[*x_v as usize] += 1;
    }

    let mut mode: u8 = 4;
    let mut mode_counts = 0;
    for i in 0..=2 {
        if counts_allele[i] > mode_counts {
            mode_counts = counts_allele[i];
            mode = i as u8;
        }
    }

    assert_ne!(mode, 4);

    for x_v in x.iter_mut() {
        if *x_v == 3 {
            *x_v = mode;
        }
    }
}

// TODO: this might be clear if use .flat_map()
// no rayon here
pub fn load_x(x: &mut [u8], buf_mi: &[B8_2], use_samples: &[bool]) {
    let mut ni = 0;
    for (n_in_i, v) in use_samples.iter().enumerate() {
        if *v {
            x[ni] = buf_to_count(buf_mi, n_in_i);
            ni += 1;
        }
    }
    // ng: x could be larger thant n
    //assert_eq!(ni, x.len());

    missing_to_mode(x);
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

        let m: usize = compute_num_snv(&fin).unwrap();
        let n: usize = compute_num_sample(&fin).unwrap();
        log::debug!("m,n: {},{}", m, n);
        let bed_size = check_valid_bed(&fin, None, m, n).unwrap();
        assert_eq!(bed_size, calculate_bed_size(m, n));
    }

    #[test]
    #[should_panic]
    fn test_check_valid_bed_panic_nofile() {
        // file does not exist
        let fin = PathBuf::from("./toy");

        let m: usize = compute_num_snv(&fin).unwrap();
        let n: usize = compute_num_sample(&fin).unwrap();
        check_valid_bed(&fin, None, m, n).unwrap();

        // TODO; check violated first 3 bytes
    }

    #[test]
    #[should_panic]
    fn test_check_valid_bed_panic() {
        //pub fn check_valid_bed(fin_bed: &str, n: usize, m: usize) -> Result<usize, Box<dyn Error>> {

        let fin = PathBuf::from("../../test/data/toy1/genot");

        let m: usize = compute_num_snv(&fin).unwrap();
        let n: usize = compute_num_sample(&fin).unwrap();
        // size is wrong
        check_valid_bed(&fin, None, m, n - 1).unwrap();
    }

    #[test]
    fn test_code_to_count() {
        let vs_input = vec![0b00, 0b01, 0b10, 0b11];
        let expects = vec![2, 3, 1, 0];
        for (v, exp) in vs_input.iter().zip(expects.iter()) {
            assert_eq!(code_to_count(*v), *exp);
        }
    }

    #[test]
    #[should_panic]
    fn test_code_to_count_panic() {
        let v = 4;
        code_to_count(v);
    }

    #[test]
    fn test_byte_to_count() {
        let v = 0b00_01_10_11;
        let expects = vec![0, 1, 3, 2];

        for (i, exp) in expects.iter().enumerate() {
            assert_eq!(byte_to_count(v, i), *exp);
        }
    }

    #[test]
    fn test_buf_to_count() {
        let buf = vec![0b00_01_10_11, 0b00_00_11_01];
        // read the lowest of the first byte
        assert_eq!(buf_to_count(&buf, 0), 0);
        // read the higest of the first byte
        assert_eq!(buf_to_count(&buf, 3), 2);
        // read the lowest of the second byte
        assert_eq!(buf_to_count(&buf, 4), 3);
        // read the highest of the second byte
        assert_eq!(buf_to_count(&buf, 7), 2);
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
    fn test_missing_to_mode() {
        let mut x = vec![0, 1, 1, 2, 1, 3];
        let x_exp = vec![0, 1, 1, 2, 1, 1];

        missing_to_mode(&mut x);

        for (x_v, x_exp_v) in x.iter().zip(x_exp.iter()) {
            assert_eq!(*x_v, *x_exp_v);
        }

        // when missing is the mode
        let mut x = vec![0, 1, 1, 3, 3, 3];
        let x_exp = vec![0, 1, 1, 1, 1, 1];

        missing_to_mode(&mut x);

        for (x_v, x_exp_v) in x.iter().zip(x_exp.iter()) {
            assert_eq!(*x_v, *x_exp_v);
        }
    }
}
