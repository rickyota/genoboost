use super::{Coef, Model, Wgt, WgtKind};
use crate::wgt::SnvWgt;
use crate::{textfile, SnvId};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

//pub fn read_dir<P: AsRef<Path>>(path: P) -> Vec<String> {
pub fn read_dir(path: &Path) -> Vec<String> {
    //log::info!("read_dir {:?}", fs::read_dir(&path).unwrap().filter_map(|entry| Some(entry.file_name().to_string_lossy().into_owned())));

    //fs::read_dir(path)
    //.unwrap()
    fs::read_dir(path)
        .unwrap_or_else(|_| panic!("Directory does not exist: {}.", path.to_str().unwrap()))
        .filter_map(|entry| {
            let entry = entry.unwrap();

            // FIXME: is_file() should allow unbroken symlink...
            // This .is_symlink() allow broken symlink, so should be avoided
            // [ref](https://doc.rust-lang.org/std/path/struct.Path.html#method.is_file)
            //if entry.file_type().unwrap().is_file() {
            if entry.file_type().unwrap().is_file() || entry.file_type().unwrap().is_symlink() {
                Some(entry.file_name().to_string_lossy().into_owned())
            } else {
                None
            }
        })
        .collect()
}

pub fn get_files_wgt(dout: &Path) -> Vec<PathBuf> {
    let files_string = read_dir(dout);

    let mut files = Vec::new();
    for f in files_string.iter() {
        let mut d = dout.to_owned();
        d.push(f);
        if d.extension().unwrap() == "wgt" {
            files.push(d);
        }
    }
    files
}

pub fn check_file_wgt_exist(fwgt: &Path) {
    let exist_fwgt = textfile::exist_file(&fwgt);
    if !exist_fwgt {
        panic!("fwgt does not exist: {:?}.", &fwgt);
    }
}

//pub fn load_wgts_file(fwgt: &Path, use_snv_pos: bool, is_nonadd: bool) -> Vec<Wgt> {
pub fn load_wgts_file(fwgt: &Path, is_nonadd: bool) -> Vec<Wgt> {
    let mut wgts: Vec<Wgt> = Vec::new();

    // also load header
    let vss: Vec<Vec<String>> = textfile::load_table(fwgt, false);

    let col_n = vss.len();
    let mut columns: Vec<String> = Vec::new();
    for col_i in 0..col_n {
        columns.push(vss[col_i][0].clone());
    }

    let mut col_to_i: HashMap<String, usize> = HashMap::new();
    for col_i in 0..col_n {
        let col = columns[col_i].clone();

        // TODO: use hashmap
        let col_name: Option<String> = match col.as_str() {
            "vid" | "sid" | "rs" | "sida" | "id" => Some("var".to_string()),
            "chrom" | "chr" => Some("chrom".to_string()),
            "pos" | "position" => Some("pos".to_string()),
            "a1" | "A1" => Some("a1".to_string()),
            "a2" | "A2" => Some("a2".to_string()),
            "score0" => Some("score0".to_string()),
            "score1" => Some("score1".to_string()),
            "score2" => Some("score2".to_string()),
            "alpha" | "wgt" => Some("alpha".to_string()),
            "a1_frq" => Some("a1_frq".to_string()),
            "a2_frq" => Some("a2_frq".to_string()),
            _ => None,
            //z => z, // unuse col
            //_ => panic!("unknown column.")
        };

        if let Some(col_name) = col_name {
            if col_to_i.contains_key(&col_name) {
                panic!("Duplicated columns in wgt file.")
            } else {
                col_to_i.insert(col_name, col_i);
            }
        }
    }
    log::debug!("col_to_i {:?}", col_to_i);

    // check all variant cols are in wgt
    let cols_wgt_required = if is_nonadd {
        vec![
            "var", "chrom", "pos", "a1", "a2", "score0", "score1", "score2",
        ]
    } else {
        vec!["var", "chrom", "pos", "a1", "a2", "alpha"]
    };
    for col_wgt in cols_wgt_required.iter() {
        if !col_to_i.contains_key(*col_wgt) {
            panic!("Required column not in wgt: {}", col_wgt);
        }
    }

    for wgt_i in 1..vss[0].len() {
        // for snv
        let snv = SnvId::new(
            vss[col_to_i["var"]][wgt_i].clone(),
            &vss[col_to_i["chrom"]][wgt_i],
            &vss[col_to_i["pos"]][wgt_i],
            vss[col_to_i["a1"]][wgt_i].clone(),
            vss[col_to_i["a2"]][wgt_i].clone(),
        );

        let coef = if is_nonadd {
            let scores = (
                vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
                vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
                vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
            );
            Coef::Score3(scores)
        } else {
            Coef::Linear(vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap())
        };
        let model = Model::new_coef(coef);

        // TODO: cleaner
        // freq might not exist
        let maf = if !col_to_i.contains_key("a1_frq") {
            None
        } else if let Some(a1_frq) = vss[col_to_i["a1_frq"]][wgt_i].parse::<f64>().ok() {
            Some(a1_frq)
        } else if !col_to_i.contains_key("a2_frq") {
            None
        } else if let Some(a2_frq) = vss[col_to_i["a2_frq"]][wgt_i].parse::<f64>().ok() {
            Some(1.0 - a2_frq)
        } else {
            None
        };
        //let maf = if let Some(a1_frq) = vss[col_to_i["a1_frq"]][wgt_i].parse::<f64>().ok() {
        //    Some(a1_frq)
        //} else if let Some(a2_frq) = vss[col_to_i["a2_frq"]][wgt_i].parse::<f64>().ok() {
        //    Some(1.0 - a2_frq)
        //} else {
        //    None
        //};
        //let a2_frq = vss[col_to_i["a2_frq"]][wgt_i].parse::<f64>().ok();
        //let maf = a2_frq.map(|x| 1.0 - x);

        let snv_wgt = SnvWgt::new_score(snv, maf);
        let wgt = Wgt::new(WgtKind::new_snv(snv_wgt), model);
        //let wgt = Wgt::new(WgtKind::Snv(snv, None, None), model);

        wgts.push(wgt);
    }

    wgts
}

/*
#[cfg(test)]
mod tests {
    use super::*;

}
*/
