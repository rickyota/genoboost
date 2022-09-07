//! module for io genotype-related module
// many func. in boosting::predict should be moved here.

use std::collections::HashMap;

pub fn create_m_to_m_in(use_snvs: &[bool]) -> (HashMap<usize, usize>, usize) {
    let mut m_to_m_in = HashMap::new();
    let mut mi = 0;
    for (m_in_i, v) in use_snvs.iter().enumerate() {
        if *v {
            m_to_m_in.insert(mi, m_in_i);
            //m_to_m_in.insert(m_in_i, mi);
            mi += 1;
        }
    }
    (m_to_m_in, mi)
}

/*
fn create_wgt_to_genotype_index(fin: &str, wgts: &[Wgt]) -> HashMap<usize, Option<usize>> {
    let m_in_wgts = wgts.len();
    let wgt_to_genotype_index = HashMap::with_capacity(m_in_wgts);

    let m_in: usize = plink::compute_num_snv(fin);
    let snvs_in = plink::load_snvs(fin, m_in);

    wgt_to_genotype_index
}
*/

#[cfg(test)]
mod tests {
    //use super::*;
}
