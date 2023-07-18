use super::{Coef, Model, Wgt, WgtKind};
use crate::{textfile, SnvId};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

pub fn read_dir<P: AsRef<Path>>(path: P) -> Vec<String> {
    //log::info!("read_dir {:?}", fs::read_dir(&path).unwrap().filter_map(|entry| Some(entry.file_name().to_string_lossy().into_owned())));

    fs::read_dir(path)
        .unwrap()
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

pub fn load_wgts_file(fwgt: &Path) -> Vec<Wgt> {
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

        let col_name = match col.as_str() {
            "vid" | "sid" | "rs" | "sida" => "var",
            "chrom" | "chr" => "chrom",
            "pos" | "position" => "pos",
            "A1" => "a1",
            "A2" => "a2",
            "wgt" => "alpha",
            z => z, // unuse col
                    //_ => panic!("unknown column.")
        }
        .to_string();

        col_to_i.insert(col_name, col_i);
    }
    log::debug!("col_to_i {:?}", col_to_i);

    for wgt_i in 1..vss[0].len() {
        // for snv
        let snv = SnvId::construct_snv_index(
            vss[col_to_i["var"]][wgt_i].clone(),
            &vss[col_to_i["chrom"]][wgt_i],
            &vss[col_to_i["pos"]][wgt_i],
            vss[col_to_i["a1"]][wgt_i].clone(),
            vss[col_to_i["a2"]][wgt_i].clone(),
        );

        let coef = Coef::Linear(vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap());
        let model = Model::new_coef(coef);

        let wgt = Wgt::construct_wgt(WgtKind::Snv(snv, None, None), model);

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
