use super::SnvId;
use crate::{textfile, vec};

use std::collections::HashMap;
use std::path::Path;

/// Use rs to extract snvs.
/// TODO: Format of fin_snv is the same as plink `--extract`.
pub fn make_use_snvs_buf(
    extract_snv_buf: Option<&[u8]>,
    snvs_in: &Vec<SnvId>,
) -> (Vec<bool>, usize) {
    //) -> (usize, Vec<bool>) {
    let m_in: usize = snvs_in.len();

    if extract_snv_buf.is_none() {
        let use_snvs = vec![true; m_in];
        return (use_snvs, m_in);
    }

    let snvs_use = load_snvs_use_buf(extract_snv_buf.unwrap());

    make_use_snvs_buf_vec(&snvs_use, snvs_in)
}

// judge by rs
fn make_use_snvs_buf_vec(snvs_use: &Vec<SnvId>, snvs_in: &Vec<SnvId>) -> (Vec<bool>, usize) {
    let m_in: usize = snvs_in.len();

    // map: sid_in -> plink index
    // &str should be fine
    let mut rs_in_to_index = HashMap::with_capacity(m_in);
    for (si, s) in snvs_in.iter().enumerate() {
        rs_in_to_index.insert(s.id(), si);
    }

    let mut use_snvs = vec![false; m_in];
    let mut m: usize = 0;
    for snv in snvs_use.iter() {
        let rs = snv.id();

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

        // if you want to panic
        //use_snvs[sida_to_index[&sida]] = true;
    }

    assert_eq!(vec::count_true(&use_snvs), m);

    (use_snvs, m)
    //(m, use_snvs)
}

pub fn load_snvs_use(fin_snv: &Path) -> Vec<SnvId> {
    let buf = textfile::read_file_to_end(fin_snv, None).unwrap();

    load_snvs_use_buf(&buf[..])
}

/// SnvId with rs only
pub fn load_snvs_use_buf(snv_buf: &[u8]) -> Vec<SnvId> {
    let mut snvs: Vec<SnvId> = vec![];

    // TODO: for .snvs with header
    // original order is chrom, rs, None, pos, A1, A2
    // rs, chrom, pos, A1, A2
    //let cols = [1usize, 0, 3, 4, 5];
    let cols = [0usize];
    let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(snv_buf, &cols, false);

    for vi in 0..vss[0].len() {
        snvs.push(SnvId::new_id(vss[0][vi].clone()));
    }
    snvs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_use_snvs_buf_vec() {
        let snvs_in = vec![
            SnvId::new("rs1".to_owned(), "1", "123", "A".to_owned(), "C".to_owned()),
            SnvId::new("rs2".to_owned(), "2", "124", "A".to_owned(), "C".to_owned()),
            SnvId::new("rs3".to_owned(), "3", "125", "A".to_owned(), "C".to_owned()),
        ];

        let snvs_use = vec![
            SnvId::new("rs3".to_owned(), "3", "125", "A".to_owned(), "C".to_owned()),
            SnvId::new("rs1".to_owned(), "1", "123", "A".to_owned(), "C".to_owned()),
        ];

        let (use_snvs, m) = make_use_snvs_buf_vec(&snvs_use, &snvs_in);

        assert_eq!(m, 2);
        assert_eq!(use_snvs, vec![true, false, true]);
    }

    #[test]
    fn test_load_snvs_use_buf() {
        let buf: Vec<u8> = "rs1\t1\t123\tA\tC\nrs2\t2\t124\tA\tC\nrs3\t2\t125\tA\tC"
            .to_string()
            .into_bytes();

        let snvs = load_snvs_use_buf(&buf);

        let snvs_ans = vec![
            SnvId::new_id("rs1".to_string()),
            SnvId::new_id("rs2".to_string()),
            SnvId::new_id("rs3".to_string()),
        ];

        assert_eq!(snvs, snvs_ans);
    }
}
