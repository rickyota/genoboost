use super::SnvId;
use crate::textfile;

use std::collections::HashMap;
use std::path::Path;

/// Use rs to extract snvs.
pub fn make_use_snvs_buf(
    extract_snv_buf: Option<&[u8]>,
    snvs_in: &Vec<SnvId>,
) -> (usize, Vec<bool>) {
    let m_in: usize = snvs_in.len();

    if extract_snv_buf.is_none() {
        let use_snvs = vec![true; m_in];
        return (m_in, use_snvs);
    }

    //textfile::check_open_file(fin_snv.unwrap());

    let mut use_snvs = vec![false; m_in];

    // map: sid_in -> plink index
    let mut rs_in_to_index = HashMap::with_capacity(m_in);

    for (si, s) in snvs_in.iter().enumerate() {
        rs_in_to_index.insert(s.rs(), si);
        //rs_in_to_index.insert(s.sida(), si);
    }

    let snvs_use = load_snvs_use_buf(extract_snv_buf.unwrap());
    //let snvs_use = load_snvs_use(fin_snv.unwrap());
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

/// Use rs to extract snvs.
/// TODO: Format of fin_snv is the same as plink `--extract`.
/// TODO: deprecate
pub fn make_use_snvs(fin_snv: Option<&Path>, snvs_in: &Vec<SnvId>) -> (usize, Vec<bool>) {
    let m_in: usize = snvs_in.len();

    if fin_snv.is_none() {
        let use_snvs = vec![true; m_in];
        return (m_in, use_snvs);
    }

    //textfile::check_open_file(fin_snv.unwrap());

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

pub fn load_snvs_use_buf(snv_buf: &[u8]) -> Vec<SnvId> {
    //textfile::check_open_file(fin_snv);

    //let m_use = textfile::compute_num_line_text(fin_snv, None).unwrap();
    //let mut snvs: Vec<SnvId> = Vec::with_capacity(m_use);
    let mut snvs: Vec<SnvId> = vec![];

    // TODO: for .snvs with header
    // original order is chrom, rs, None, pos, A1, A2
    // rs, chrom, pos, A1, A2
    //let cols = [1usize, 0, 3, 4, 5];
    let cols = [0usize];
    let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(snv_buf, &cols, false);

    for vi in 0..vss[0].len() {
        snvs.push(SnvId::construct_snv_index_rs(vss[0][vi].clone()));
    }
    snvs
}

pub fn load_snvs_use(fin_snv: &Path) -> Vec<SnvId> {
    let buf = textfile::read_file_to_end(fin_snv, None).unwrap();

    load_snvs_use_buf(&buf[..])

    /*     //textfile::check_open_file(fin_snv);

    //let m_use = textfile::compute_num_line_text(fin_snv, None).unwrap();
    //let mut snvs: Vec<SnvId> = Vec::with_capacity(m_use);
    let mut snvs: Vec<SnvId> = vec![];

    // TODO: for .snvs with header
    // original order is chrom, rs, None, pos, A1, A2
    // rs, chrom, pos, A1, A2
    //let cols = [1usize, 0, 3, 4, 5];
    let cols = [0usize];
    let vss: Vec<Vec<String>> = textfile::load_table_cols(&fin_snv, &cols, false).unwrap();

    for vi in 0..vss[0].len() {
        snvs.push(SnvId::construct_snv_index_rs(vss[0][vi].clone()));
    }
    snvs */
}
