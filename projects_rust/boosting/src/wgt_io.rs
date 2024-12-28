use genetics::WgtKind;

use crate::textfile;
//use crate::wgt_boost;
use crate::wgt_boosts::WgtBoosts;
use crate::BoostType;
use crate::DoutFile;
use std::collections::HashSet;
use std::io::BufRead;
use std::io::BufReader;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

pub fn write_file_string(fin: &Path, wgt_best: String) {
    let file = match File::create(&fin) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {:?}",
            &fin
        ),
    };

    let mut writer = BufWriter::new(file);
    writer.write(wgt_best.as_bytes()).unwrap();
}

pub fn create_best_para_buf(nsnv_acc_max: usize, lr: f64) -> String {
    /*     fn lr_to_str(lr: &Option<f64>) -> String {
        match lr {
            None => "1.0".to_string(),
            Some(x) => x.to_string(),
        }
    } */
    let wgt_best_para = vec![
        String::from("nsnv\t") + nsnv_acc_max.to_string().as_str(),
        String::from("learning_rate\t") + lr.to_string().as_str(),
    ];

    let wgt_best_para = wgt_best_para.join("\n");

    wgt_best_para
}

pub fn create_best_para_integrate_buf(
    nsnv_acc_max: usize,
    lr: f64,
    boost_type: BoostType,
) -> String {
    let wgt_best_para = vec![
        String::from("nsnv\t") + nsnv_acc_max.to_string().as_str(),
        String::from("learning_rate\t") + lr.to_string().as_str(),
        String::from("genetic_model\t") + boost_type.genetic_model().as_str(),
    ];

    let wgt_best_para = wgt_best_para.join("\n");

    wgt_best_para
}

pub fn create_best_wgt_buf(
    dout: &DoutFile,
    //dout: &Path,
    nsnv_acc_max: usize,
    learning_rate: f64,
    boost_type: Option<BoostType>,
) -> String {
    let dout_para = match boost_type {
        Some(boost_type) => dout.dout_para_integrate(learning_rate, boost_type),
        //Some(boost_type) => wgt_boost::io::get_dname_para_integrate(dout, lr, boost_type),
        None => dout.dout_para_lr(learning_rate),
        //wgt_boost::io::get_dname_para(dout, lr),
    };

    let fwgt = dout_para.get_file_wgt();
    //let fwgt = wgt_boost::io::get_fname_wgt(&dout_para);

    let wgts = WgtBoosts::new_from_file(&fwgt);

    // get iteration number
    let mut itern: Option<usize> = None;
    //let mut itern: usize;

    let mut rs_set = HashSet::new();
    for (wgti, wgt) in wgts.wgts().iter().enumerate() {
        match wgt.wgt().kind() {
            //WgtKind::Snv(snv_id, _, _) => {
            WgtKind::Snv(snv_wgt) => {
                let snv_id = snv_wgt.snv_id();
                let rs = snv_id.id().to_string();
                rs_set.insert(rs);
            }
            //WgtKind::SnvInteraction(snv_id_1, _, snv_id_2, _) => {
            WgtKind::SnvInteraction(snv_inter_wgt) => {
                let (snv_id_1, snv_id_2) = snv_inter_wgt.snv_ids();
                let rs_1 = snv_id_1.id().to_string();
                let rs_2 = snv_id_2.id().to_string();
                rs_set.insert(rs_1);
                rs_set.insert(rs_2);
            }
            WgtKind::Cov(_) => {}
        }

        //if !wgt.wgt().kind().is_snv() {
        //    continue;
        //}
        //let rs = wgt.wgt().kind().snv_index().id().to_string();
        //rs_set.insert(rs);

        // first get wgt of nsnv+1
        if rs_set.len() == nsnv_acc_max + 1 {
            // wgti -1 : last wgt index = iter is just before
            // wgti -1 +1: last wgt index to number of snvs
            itern = Some(wgti);
            //itern = Some(wgti - 1);
            break;
        }
    }

    log::debug!("rs_set len={:?}", rs_set.len());

    log::debug!("itern={:?}", itern);
    log::debug!("rs_set={:?}", rs_set);

    if itern.is_none() {
        panic!("Error in creating best_wgt_buf.");
    }

    // +1 for header
    let linen = itern.unwrap() + 1;

    let buf_wgt = textfile::read_file_to_end(&fwgt, None).unwrap();

    let lines = BufReader::new(&buf_wgt[..]).lines().take(linen);
    let wgt_best_vec_string: Vec<String> = lines.into_iter().map(|x| x.unwrap()).collect();

    let wgt_best_string = wgt_best_vec_string.join("\n");

    wgt_best_string
}
