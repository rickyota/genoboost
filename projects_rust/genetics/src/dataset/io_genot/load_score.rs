use super::load;
use crate::genot::prelude::*;
use crate::Genot;
use crate::GenotFormat;
use crate::SnvId;
use crate::WgtTrait;
use crate::Wgts;
use crate::{io_genot, vec};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::Path;

// move to boosting_rust
pub fn load_genotypes_for_score<W: WgtTrait>(
    fin: &Path,
    gfmt: GenotFormat,
    wgts: &mut [W],
    //wgts: &mut [Wgt],
    n: usize,
    use_samples: Option<&[bool]>,
    use_missing: bool,
) -> Genot {
    let (use_snvs, is_reversed) = create_wgt_to_genotype_index_sametime(fin, gfmt, wgts);
    //let (wgt_to_m, use_snvs, is_flipped) = create_wgt_to_genotype_index(fin, wgts);

    let m = vec::count_true(&use_snvs);

    let mut g = load::generate_genot(fin, gfmt, m, n, &use_snvs, use_samples, use_missing);

    reverse_genotypes_rayon(&mut g.as_genot_mut(), &is_reversed);

    g
}

fn integrate_use_snvs(uses_snvs: &[Vec<bool>]) -> Vec<bool> {
    let mut use_snvs = uses_snvs[0].clone();
    for use_snv_ in uses_snvs.iter() {
        vec::or_bool_vec_mut(&mut use_snvs, use_snv_);
    }
    use_snvs
}

fn integrate_is_reversed(iss_reversed: &[Vec<bool>]) -> Vec<bool> {
    let mut is_reversed = iss_reversed[0].clone();
    for is_reversed_ in iss_reversed.iter() {
        vec::or_bool_vec_mut(&mut is_reversed, is_reversed_);
    }
    is_reversed
}

fn check_consistent_is_reversed(
    is_reversed: &[bool],
    iss_reversed: &[Vec<bool>],
    useds_snvs: &[Vec<bool>],
) {
    for (is_reversed_, use_snv_) in iss_reversed.iter().zip(useds_snvs.iter()) {
        for (snv_i, &use_snv_val) in use_snv_.iter().enumerate() {
            if use_snv_val {
                if is_reversed_[snv_i] != is_reversed[snv_i] {
                    panic!("Alleles are reversed within fwgt in dwgt.");
                }
            }
        }
    }
}

/*
fn integrate_use_snvs_reversed(
    uses_snvs: Vec<Vec<bool>>,
    iss_reversed: Vec<Vec<bool>>,
) -> (Vec<bool>, Vec<bool>) {
    let mut use_snvs = uses_snvs[0].clone();
    for use_snv_ in uses_snvs.iter() {
        vec::or_bool_vec_mut(&mut use_snvs, use_snv_);
    }

    let m = vec::count_true(&use_snvs);

    let mut is_reversed = vec![false; m];
    let (m_in_to_m, _) = vec::create_hashmap_index_from_use(&use_snvs);

    // each is_reversed in iss_ has different length.
    for (is_reversed_, use_snv_) in iss_reversed.iter().zip(uses_snvs.iter()) {
        let mut rev_i = 0;
        for (snv_i, &use_snv_val) in use_snv_.iter().enumerate() {
            if use_snv_val {
                let m_i = m_in_to_m[&snv_i];
                if is_reversed_[rev_i] {
                    is_reversed[m_i] = true;
                }
                rev_i += 1;
            }
        }
    }

    // check all iss_reversed are consistent
    for (is_reversed_, use_snv_) in iss_reversed.iter().zip(uses_snvs.iter()) {
        let mut rev_i = 0;
        for (snv_i, &use_snv_val) in use_snv_.iter().enumerate() {
            if use_snv_val {
                let m_i = m_in_to_m[&snv_i];
                assert_eq!(is_reversed_[rev_i], is_reversed[m_i]);
                rev_i += 1;
            }
        }
    }

    (use_snvs, is_reversed)
}
 */

// move to boosting_rust
pub fn load_genotypes_for_score_multiwgts(
    fin: &Path,
    gfmt: GenotFormat,
    wgts_multi: &mut [Wgts],
    //wgts: &mut [Wgt],
    n: usize,
    use_samples: Option<&[bool]>,
    use_missing: bool,
) -> Genot {
    let snvs_in = io_genot::load_snvs(fin, gfmt);
    let mut uses_snvs = Vec::new();
    for wgts in wgts_multi.iter_mut() {
        // length are m_in
        let use_snvs = create_wgt_to_genotype_index(&snvs_in, wgts.wgts());
        //let use_snvs = create_wgt_to_genotype_index(fin, gfmt, wgts.wgts());
        uses_snvs.push(use_snvs);
    }
    let use_snvs = integrate_use_snvs(&uses_snvs);

    //let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
    //let snvs_in = io_genot::load_snvs(fin, gfmt);
    // TODO: how to avoid snvid.clone()?
    let snvs_use: Vec<SnvId> = snvs_in
        .iter()
        .zip(use_snvs.iter())
        .filter(|(_, &b)| b)
        .map(|(snvid, _)| snvid.clone())
        .collect();

    let mut iss_reversed = Vec::new();
    let mut useds_snvs = Vec::new();
    for wgts in wgts_multi.iter_mut() {
        // length are m
        let (is_reversed, used_snvs) = set_wgts_index_reversed(wgts.wgts_mut(), &snvs_use);
        iss_reversed.push(is_reversed);
        useds_snvs.push(used_snvs);
    }
    let is_reversed = integrate_is_reversed(&iss_reversed);

    log::debug!("Reversed snvs: {}", vec::count_true(&is_reversed));

    check_consistent_is_reversed(&is_reversed, &iss_reversed, &useds_snvs);

    let m = vec::count_true(&use_snvs);

    let mut g = load::generate_genot(fin, gfmt, m, n, &use_snvs, use_samples, use_missing);

    reverse_genotypes_rayon(&mut g.as_genot_mut(), &is_reversed);

    g
}

/// use in score calculation
/// Judge the corresponding SNVs by chrom, pos, A1, A2.
/// Reversed wgts should have reversed genotypes
/// Flipped genotype do nth.
///
/// genotype_index should be in corresponding order to fin
///     (wgt_to_m, use_snvs, is_flipped): len: .. , m_in, m
///
/// The reason why is_flipped and set_index is here is they require snvs
///
fn create_wgt_to_genotype_index<W: WgtTrait>(
    snvs_in: &[SnvId],
    //fin: &Path,
    //gfmt: GenotFormat,
    wgts: &[W],
) -> Vec<bool> {
    //let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
    //let snvs_in = io_genot::load_snvs(fin, gfmt);
    let m_in = snvs_in.len();

    let mut use_snvs = vec![false; m_in];

    // map: sid_in -> genotype index (m)
    let mut sid_in_to_m_in = HashMap::with_capacity(m_in);

    for (si, s) in snvs_in.iter().enumerate() {
        sid_in_to_m_in.insert(s.to_sid(), si);
        //sid_in_to_m_in.insert(s.sida(), si);
    }

    // create use_snvs
    // Even though snv's alelles do not match, use_snvs is true, but no problem.
    for wgt in wgts.iter() {
        // if not snv, skip
        if wgt.kind().is_snv() {
            let sid = wgt.kind().snv_index().to_sid();
            if let Some(v) = sid_in_to_m_in.get(&sid) {
                use_snvs[*v] = true;
            }
        }
    }

    use_snvs
}

/// use in score calculation
/// Judge the corresponding SNVs by chrom, pos, A1, A2.
/// Reversed wgts should have reversed genotypes
/// Flipped genotype do nth.
///
/// genotype_index should be in corresponding order to fin
///     (wgt_to_m, use_snvs, is_flipped): len: .. , m_in, m
///
/// The reason why is_flipped and set_index is here is they require snvs
///
//fn create_wgt_to_genotype_index<W: WgtTrait>(
//    fin: &Path,
//    gfmt: GenotFormat,
//    wgts: &[W],
//) -> Vec<bool> {
//    let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
//    let snvs_in = io_genot::load_snvs(fin, gfmt);
//
//    let mut use_snvs = vec![false; m_in];
//
//    // map: sid_in -> genotype index (m)
//    let mut sid_in_to_m_in = HashMap::with_capacity(m_in);
//
//    for (si, s) in snvs_in.iter().enumerate() {
//        sid_in_to_m_in.insert(s.to_sid(), si);
//        //sid_in_to_m_in.insert(s.sida(), si);
//    }
//
//    // create use_snvs
//    // Even though snv's alelles do not match, use_snvs is true, but no problem.
//    for wgt in wgts.iter() {
//        // if not snv, skip
//        if wgt.kind().is_snv() {
//            let sid = wgt.kind().snv_index().to_sid();
//            if let Some(v) = sid_in_to_m_in.get(&sid) {
//                use_snvs[*v] = true;
//            }
//        }
//    }
//
//    /*
//    let (m_in_to_m, m) = vec::create_hashmap_index_from_use(&use_snvs);
//
//    // needs to reverse genotype
//    let mut is_reversed = vec![false; m];
//
//    // TODO
//    //  Not here, want to make load_dataset(&wgts) not &mut
//    // boosting.rs:wgt_sida_to_dataset_m()
//    // create wgt_to_m
//    // is_reversed will be overwritten many times but should be all right if all of the same in wgt
//    for wgt in wgts.iter_mut() {
//        if wgt.kind().is_snv() {
//            let sid = wgt.kind().snv_index().to_sid();
//            //log::debug!("sida {}", sida);
//            //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
//            if let Some(v) = sid_in_to_m_in.get(&sid) {
//                let mi = m_in_to_m[v];
//
//                // clone is necessary
//                let snv_in = snvs_in[*v].clone();
//                // clone is necessary
//                let snv_wgt = wgt.kind().snv_index().clone();
//
//                // check alleles match
//                // otherwise alleles do not match
//                if &snv_wgt == &snv_in {
//                    wgt.set_snv_index(Some(mi));
//
//                    // rev or not
//                    // here since putting inside of if raises error
//                    let is_rev = snv_wgt.is_rev(&snv_in);
//
//                    if is_rev {
//                        log::debug!(
//                            "Alleles are reversed in wgt and fplink {:?}, {:?}",
//                            &snv_wgt,
//                            &snv_in
//                        );
//                    }
//
//                    is_reversed[mi] = is_rev;
//                    //is_reversed[mi] = false;
//                } else {
//                    log::info!(
//                        "Ignore SNV: alleles do not match in wgt and fplink {:?}, {:?}",
//                        snv_wgt,
//                        &snv_in
//                    );
//                }
//            } else {
//                wgt.set_snv_index(None);
//            }
//        }
//
//        // do nth for cov
//        // not found
//        //wgt_to_m.insert(wgti, None);
//    }
//    */
//
//    use_snvs
//}

fn set_wgts_index_reversed<W: WgtTrait>(wgts: &mut [W], snvs: &[SnvId]) -> (Vec<bool>, Vec<bool>) {
    let m = snvs.len();

    // map: sid_in -> genotype index (m)
    let mut sid_to_m = HashMap::with_capacity(m);

    for (si, s) in snvs.iter().enumerate() {
        sid_to_m.insert(s.to_sid(), si);
    }

    // needs to reverse genotype
    let mut is_reversed = vec![false; m];

    // for check_reversed
    let mut used_snvs = vec![false; m];

    // TODO
    //  Not here, want to make load_dataset(&wgts) not &mut
    // boosting.rs:wgt_sida_to_dataset_m()
    // create wgt_to_m
    // is_reversed will be overwritten many times but should be all right if all of the same in wgt
    for wgt in wgts.iter_mut() {
        if wgt.kind().is_snv() {
            let sid = wgt.kind().snv_index().to_sid();
            //log::debug!("sida {}", sida);
            //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
            if let Some(&mi) = sid_to_m.get(&sid) {
                //let mi = m_in_to_m[v];

                used_snvs[mi] = true;

                // clone is necessary
                let snv_in = snvs[mi].clone();
                // clone is necessary
                let snv_wgt = wgt.kind().snv_index().clone();

                // check alleles match
                // otherwise alleles do not match
                if &snv_wgt == &snv_in {
                    wgt.set_snv_index(Some(mi));

                    // rev or not
                    // here since putting inside of if raises error
                    let is_rev = snv_wgt.is_rev(&snv_in);

                    if is_rev {
                        log::debug!(
                            "Alleles are reversed in wgt and fplink {:?}, {:?}",
                            &snv_wgt,
                            &snv_in
                        );
                    }

                    is_reversed[mi] = is_rev;
                    //is_reversed[mi] = false;
                } else {
                    log::info!(
                        "Ignore SNV: alleles do not match in wgt and fplink {:?}, {:?}",
                        snv_wgt,
                        &snv_in
                    );
                }
            } else {
                wgt.set_snv_index(None);
            }
        }

        // do nth for cov
        // not found
        //wgt_to_m.insert(wgti, None);
    }

    (is_reversed, used_snvs)
}

/// old use above
///
/// use in score calculation
/// Judge the corresponding SNVs by chrom, pos, A1, A2.
/// Reversed wgts should have reversed genotypes
/// Flipped genotype do nth.
///
/// genotype_index should be in corresponding order to fin
///     (wgt_to_m, use_snvs, is_flipped): len: .. , m_in, m
///
/// The reason why is_flipped and set_index is here is they require snvs
///
fn create_wgt_to_genotype_index_sametime<W: WgtTrait>(
    fin: &Path,
    gfmt: GenotFormat,
    wgts: &mut [W],
) -> (Vec<bool>, Vec<bool>) {
    let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
    let snvs_in = io_genot::load_snvs(fin, gfmt);

    let mut use_snvs = vec![false; m_in];

    // map: sid_in -> genotype index (m)
    let mut sid_in_to_m_in = HashMap::with_capacity(m_in);

    for (si, s) in snvs_in.iter().enumerate() {
        sid_in_to_m_in.insert(s.to_sid(), si);
        //sid_in_to_m_in.insert(s.sida(), si);
    }

    // create use_snvs
    // Even though snv's alelles do not match, use_snvs is true, but no problem.
    for wgt in wgts.iter() {
        // if not snv, skip
        if wgt.kind().is_snv() {
            let sid = wgt.kind().snv_index().to_sid();
            if let Some(v) = sid_in_to_m_in.get(&sid) {
                use_snvs[*v] = true;
            }
        }
    }

    let (m_in_to_m, m) = vec::create_hashmap_index_from_use(&use_snvs);

    // needs to reverse genotype
    let mut is_reversed = vec![false; m];

    // TODO
    //  Not here, want to make load_dataset(&wgts) not &mut
    // boosting.rs:wgt_sida_to_dataset_m()
    // create wgt_to_m
    // is_reversed will be overwritten many times but should be all right if all of the same in wgt
    for wgt in wgts.iter_mut() {
        if wgt.kind().is_snv() {
            let sid = wgt.kind().snv_index().to_sid();
            //log::debug!("sida {}", sida);
            //log::debug!("m_in {:?}", sida_in_to_m_in.get(sida));
            if let Some(v) = sid_in_to_m_in.get(&sid) {
                let mi = m_in_to_m[v];

                // clone is necessary
                let snv_in = snvs_in[*v].clone();
                // clone is necessary
                let snv_wgt = wgt.kind().snv_index().clone();

                // check alleles match
                // otherwise alleles do not match
                if &snv_wgt == &snv_in {
                    wgt.set_snv_index(Some(mi));

                    // rev or not
                    // here since putting inside of if raises error
                    let is_rev = snv_wgt.is_rev(&snv_in);

                    if is_rev {
                        log::debug!(
                            "Alleles are reversed in wgt and fplink {:?}, {:?}",
                            &snv_wgt,
                            &snv_in
                        );
                    }

                    is_reversed[mi] = is_rev;
                    //is_reversed[mi] = false;
                } else {
                    log::info!(
                        "Ignore SNV: alleles do not match in wgt and fplink {:?}, {:?}",
                        snv_wgt,
                        &snv_in
                    );
                }
            } else {
                wgt.set_snv_index(None);
            }
        }

        // do nth for cov
        // not found
        //wgt_to_m.insert(wgti, None);
    }

    (use_snvs, is_reversed)
}

/// rayon ver.
fn reverse_genotypes_rayon(g: &mut GenotMut, is_reversed: &[bool]) {
    g.iter_snv_mut()
        .zip(is_reversed.iter())
        .par_bridge()
        .filter(|(_, b)| **b)
        .for_each(|(mut gsnv, _)| reverse_genotypes_snv(&mut gsnv));
}

// {0:2, 1:1, 2:0, 3:3}
const REVERSE_GENOTYPE_AR: [u8; 4] = [2, 1, 0, 3];

#[inline]
fn reverse_genotypes_snv(gsnv: &mut GenotSnvMut) {
    for ni in 0..gsnv.n() {
        let v = gsnv.get_val_unchecked(ni) as usize;
        gsnv.set_unchecked(REVERSE_GENOTYPE_AR[v], ni);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let vec = vec![1, 2, 3, 0];
        let mut g = GenotSnv::new(&vec);
        reverse_genotypes_snv(&mut g.as_genot_snv_mut_snv());
        assert_eq!(g.vals(), vec![1, 0, 3, 2]);
    }

    #[test]
    fn test_check_consistent_is_reversed() {
        let useds_snvs = vec![
            vec![true, true, false, false, false],
            vec![false, false, true, true, true],
            vec![false, true, true, false, false],
        ];

        // 0,2 are flipped
        let iss_reversed = vec![
            vec![true, false, false, false, false],
            vec![false, false, true, false, false],
            vec![false, false, true, false, false],
        ];

        let is_reversed = vec![true, false, true, false, false];

        check_consistent_is_reversed(&is_reversed, &iss_reversed, &useds_snvs);
    }

    #[should_panic]
    #[test]
    fn test_check_consistent_is_reversed_panic() {
        let useds_snvs = vec![
            vec![true, true, false, false, false],
            vec![false, false, true, true, true],
            vec![false, true, true, false, false],
        ];

        // 0,2th snv are flipped
        // but 1st wgt says 2nd snv is not flipped
        let iss_reversed = vec![
            vec![true, false, false, false, false],
            vec![false, false, false, false, false],
            vec![false, false, true, false, false],
        ];

        let is_reversed = vec![true, false, true, false, false];

        check_consistent_is_reversed(&is_reversed, &iss_reversed, &useds_snvs);
    }

    /*
    #[test]
    fn test_integrate_use_snvs_reversed() {
        let uses_snvs = vec![
            vec![false, true, false, false, false],
            vec![false, false, true, false, false],
            vec![false, true, true, false, false],
        ];

        let iss_reversed = vec![vec![false], vec![true], vec![false, true]];

        let use_snvs_ans = vec![false, true, true, false, false];
        let is_reversed_ans = vec![false, true];

        let (use_snvs, is_reversed) = integrate_use_snvs_reversed(uses_snvs, iss_reversed);

        assert_eq!(use_snvs, use_snvs_ans);
        assert_eq!(is_reversed, is_reversed_ans);
    }

    #[should_panic]
    #[test]
    fn test_integrate_use_snvs_reversed_panic() {
        let uses_snvs = vec![
            vec![false, true, false, false, false],
            vec![false, false, true, false, false],
            vec![false, true, true, false, false],
        ];

        let iss_reversed = vec![
            vec![false],
            vec![true],
            // this is not consistent so should panic
            vec![false, false],
        ];

        let use_snvs_ans = vec![false, true, true, false, false];
        let is_reversed_ans = vec![false, true];

        let (use_snvs, is_reversed) = integrate_use_snvs_reversed(uses_snvs, iss_reversed);

        assert_eq!(use_snvs, use_snvs_ans);
        assert_eq!(is_reversed, is_reversed_ans);
    }
    */
}
