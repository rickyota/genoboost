use super::load;
use crate::genot::prelude::*;
use crate::Genot;
use crate::WgtTrait;
use crate::{plink, vec};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::Path;

// move to boosting_rust
pub fn load_genotypes_for_score<W: WgtTrait>(
    fin: &Path,
    wgts: &mut [W],
    //wgts: &mut [Wgt],
    n: usize,
    use_samples: Option<&[bool]>,
    use_missing: bool,
) -> Genot {
    let (use_snvs, is_flipped) = create_wgt_to_genotype_index(fin, wgts);
    //let (wgt_to_m, use_snvs, is_flipped) = create_wgt_to_genotype_index(fin, wgts);

    let m = vec::count_true(&use_snvs);

    let mut g = load::generate_genot(fin, m, n, &use_snvs, use_samples, use_missing);

    // FIXME: flip if wgt a1,a2
    flip_genotypes_rayon(&mut g.as_genot_mut(), &is_flipped);

    g
}

// FIXME: to check flipped A1 and A2
/// use in score calculation
/// genotype_index should be in corresponding order to fin
///     (wgt_to_m, use_snvs, is_flipped): len: .. , m_in, m
fn create_wgt_to_genotype_index<W: WgtTrait>(fin: &Path, wgts: &mut [W]) -> (Vec<bool>, Vec<bool>) {
    let m_in: usize = plink::compute_num_snv(fin).unwrap();
    let snvs_in = plink::load_snvs(fin, m_in);

    let mut use_snvs = vec![false; m_in];

    // map: sid_in -> genotype index (m)
    let mut sida_in_to_m_in = HashMap::with_capacity(m_in);

    for (si, s) in snvs_in.iter().enumerate() {
        sida_in_to_m_in.insert(s.sida(), si);
    }

    // create use_snvs
    for wgt in wgts.iter() {
        // if not snv, skip
        if wgt.kind().is_snv() {
            let sida = wgt.kind().snv_index().sida();
            if let Some(v) = sida_in_to_m_in.get(sida) {
                use_snvs[*v] = true;
            }
        }
    }

    let (m_in_to_m, m) = vec::create_hashmap_index_from_use(&use_snvs);

    let mut is_flipped = Vec::with_capacity(m);
    for _ in 0..m {
        is_flipped.push(false);
    }

    // TODO
    //  Not here, want to make load_dataset(&wgts) not &mut
    // boosting.rs:wgt_sida_to_dataset_m()
    // create wgt_to_m
    for wgt in wgts.iter_mut() {
        if wgt.kind().is_snv() {
            let sida = wgt.kind().snv_index().sida();
            //println!("sida {}", sida);
            //println!("m_in {:?}", sida_in_to_m_in.get(sida));
            if let Some(v) = sida_in_to_m_in.get(sida) {
                let mi = m_in_to_m[v];
                wgt.set_snv_index(Some(mi));

                // TODO: when flipped
                is_flipped[mi] = false;
            } else {
                wgt.set_snv_index(None);
            }
        }

        // do nth for cov
        // not found
        //wgt_to_m.insert(wgti, None);
    }

    (use_snvs, is_flipped)
}

/// rayon ver.
fn flip_genotypes_rayon(g: &mut GenotMut, is_flipped: &[bool]) {
    //let len_n = len_n(n);
    g.iter_snv_mut()
        .zip(is_flipped.iter())
        .par_bridge()
        .filter(|(_, b)| **b)
        .for_each(|(mut gsnv, _)| flip_genotypes_snv(&mut gsnv));
    /*     genotypes
    .par_chunks_mut(len_n)
    .enumerate()
    .filter(|&(mi, _)| is_flipped[mi])
    .for_each(|(_, genot)| flip_genotypes_snv(genot)); */
}

// {0:2, 1:1, 2:0, 3:3}
const FLIP_GENOTYPE_AR: [u8; 4] = [2, 1, 0, 3];

#[inline]
fn flip_genotypes_snv(gsnv: &mut GenotSnvMut) {
    for ni in 0..gsnv.n() {
        let v = gsnv.get_val_unchecked(ni) as usize;
        gsnv.set_unchecked(FLIP_GENOTYPE_AR[v], ni);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let vec = vec![1, 2, 3, 0];
        let mut g = GenotSnv::new(vec);
        flip_genotypes_snv(&mut g.as_genot_snv_mut_snv());
        assert_eq!(g.vals(), vec![1, 0, 3, 2]);
    }
}
