pub mod plink;
#[cfg(feature = "plink2")]
pub mod plink2;

use crate::genot::prelude::*;
use crate::GenotFormat;
use std::path::Path;

// TODO: unused now: for small number of snvs called from python
// TODO: mi for multiple files
// TODO: remove n
/// mi is in fin file
pub fn generate_genot_snv(
    fin: &Path,
    gfmt: GenotFormat,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    //use_missing: bool,
    fill_missing: bool,
) -> GenotSnv {
    match gfmt {
        GenotFormat::Plink1 => {
            plink::generate_genot_snv_plink(fin, gfmt, mi, n, use_samples, fill_missing)
        }
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
            call_generate_genot_snv_plink2(fin, gfmt, mi, n, use_samples, fill_missing)
            //if cfg!(feature = "plink2") {
            //    plink2::load_snv_plink2(fin, gfmt, mi, n, use_samples, use_missing)
            //} else {
            //    panic!("Cannot use plink2 in this program feature. Use --feature plink2");
            //}
        }
    }
}

#[cfg(feature = "plink2")]
fn call_generate_genot_snv_plink2(
    fin: &Path,
    gfmt: GenotFormat,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    fill_missing: bool,
) -> GenotSnv {
    plink2::generate_genot_snv_plink2(fin, gfmt, mi, n, use_samples, fill_missing)
}

#[cfg(not(feature = "plink2"))]
fn call_generate_genot_snv_plink2(
    _: &Path,
    _: GenotFormat,
    _: usize,
    _: usize,
    _: Option<&[bool]>,
    _: bool,
) -> GenotSnv {
    panic!("Cannot use plink2 in this program feature. Use --feature plink2");
}

/// Generate Vector of predictions.
pub fn generate_genot(
    fin: &Path,
    gfmt: GenotFormat,
    m: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    //use_missing: bool,
    fill_missing: bool,
) -> Genot {
    match gfmt {
        GenotFormat::Plink1 => {
            plink::generate_genot_plink(fin, gfmt, m, n, use_snvs, use_samples, fill_missing)
        }
        GenotFormat::Plink2 | GenotFormat::Plink2Vzs => {
            call_generate_genot_plink2(fin, gfmt, m, n, use_snvs, use_samples, fill_missing)
        }
    }
}

/* // error: plink2 is not loaded
fn call_generate_genot_plink2(
    fin: &Path,
    gfmt: GenotFormat,
    m: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    use_missing: bool,
) -> Genot {
    if cfg!(feature = "plink2") {
        plink2::generate_genot_plink2(fin, gfmt, m, n, use_snvs, use_samples, use_missing)
    } else {
        panic!("Cannot use plink2 in this program feature. Use --feature plink2");
    }
} */

#[cfg(feature = "plink2")]
fn call_generate_genot_plink2(
    fin: &Path,
    gfmt: GenotFormat,
    m: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    fill_missing: bool,
) -> Genot {
    plink2::generate_genot_plink2(fin, gfmt, m, n, use_snvs, use_samples, fill_missing)
}

// TODO: cleaner
#[cfg(not(feature = "plink2"))]
fn call_generate_genot_plink2(
    _: &Path,
    _: GenotFormat,
    _: usize,
    _: usize,
    _: Option<&[bool]>,
    _: Option<&[bool]>,
    _: bool,
) -> Genot {
    panic!("Cannot use plink2 in this program feature. Use --feature plink2");
}

// make missing to mode
// in pred
// TODO: mv to BaseGenotSnvMut
//fn fill_missing_snv(pred: &mut GenotSnvMut) {
//    pred.fill_missing_mode()
//    //// count 0,1,2
//    //let mut counts_allele = vec![0usize; 4];
//    ////let mut counts_allele = Vec::with_capacity(4);
//    ////for _ in 0..=3 {
//    ////    counts_allele.push(0);
//    ////}
//
//    //let n = pred.n();
//    //for ni in 0..n {
//    //    counts_allele[pred.get_val_unchecked(ni) as usize] += 1;
//    //}
//
//    //let mut mode: usize = 4;
//    //let mut mode_counts = 0;
//    //for i in 0..=2 {
//    //    if counts_allele[i] > mode_counts {
//    //        mode_counts = counts_allele[i];
//    //        mode = i;
//    //    }
//    //}
//    //let mode = mode as u8;
//    //assert_ne!(mode, 4);
//
//    //for ni in 0..n {
//    //    if pred.get_val_unchecked(ni) == 3 {
//    //        pred.set(mode, ni);
//    //    }
//    //}
//    //// check all are non-missing?
//    //// -> performance...
//}

// make missing to mode
// in pred
// TODO: mv to BaseGenotSnvMut
//pub fn fill_missing_maf(pred: &mut GenotSnvMut, maf: f64) {
//    pred.fill_missing_mode_maf(maf);
//
//    //let mode: u8 = if maf < 1.0 / 3.0f64 {
//    //    0
//    //} else if maf > 2.0 / 3.0f64 {
//    //    2
//    //} else {
//    //    1
//    //};
//
//    //let n = pred.n();
//
//    //for ni in 0..n {
//    //    if pred.get_val_unchecked(ni) == 3 {
//    //        pred.set(mode, ni);
//    //    }
//    //}
//    // check all are non-missing?
//    // -> performance...
//}

// not used?
/// x: vector of minor allele (0,1,2,3)
/// convert 3 (missing) to mode
//fn missing_to_mode(x: &mut [u8]) {
//    // count 0,1,2,3
//    let mut counts_allele = vec![0usize; 4];
//    //let mut counts_allele = Vec::with_capacity(4);
//    //for _ in 0..=3 {
//    //    counts_allele.push(0);
//    //}
//
//    for x_v in x.iter() {
//        counts_allele[*x_v as usize] += 1;
//    }
//
//    let mut mode: u8 = 4;
//    let mut mode_counts = 0;
//    for i in 0..=2 {
//        if counts_allele[i] > mode_counts {
//            mode_counts = counts_allele[i];
//            mode = i as u8;
//        }
//    }
//
//    assert_ne!(mode, 4);
//
//    for x_v in x.iter_mut() {
//        if *x_v == 3 {
//            *x_v = mode;
//        }
//    }
//}

#[cfg(test)]
mod tests {
    use crate::GenotFormat;

    use super::*;
    use crate::samples;
    use crate::{io_genot, sample, snv};
    use std::path::PathBuf;

    /// check if _file() and _file_loadpart() are the same

    #[allow(dead_code)]
    fn setup_test() -> (PathBuf, Vec<bool>, usize, usize, Vec<bool>, Vec<bool>) {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let gfmt = GenotFormat::Plink1;
        let fin_snv = None;
        let fin_sample = None;

        let m_in: usize = io_genot::compute_num_snv(&fin, gfmt).unwrap();
        log::debug!("{}", m_in);
        let n_in: usize = io_genot::compute_num_sample(&fin, gfmt).unwrap();
        println!("{}", n_in);
        // load snvs
        let snvs_in = io_genot::load_snvs(&fin, gfmt);
        let (m, use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (n, use_samples) = sample::make_use_samples(fin_sample, &fin, gfmt);
        //let ys: Vec<bool> = io_genot::load_ys(&fin, gfmt, None, None, &use_samples);
        let sample_id_to_n = samples::create_sample_id_to_n(&fin, gfmt, Some(&use_samples));
        let ys: Vec<bool> = io_genot::load_ys_buf(&fin, gfmt, None, None, &sample_id_to_n);

        (fin, ys, m, n, use_snvs, use_samples)
    }

    #[allow(dead_code)]
    fn setup_test2() -> (PathBuf, Vec<bool>, usize, usize, Vec<bool>, Vec<bool>) {
        //let fin = String::from("../../test/data/toy1/genot");
        let fin = PathBuf::from("../../test/data/1kg_n10000/1kg_n10000");
        let gfmt = GenotFormat::Plink1;
        let fin_snv = None;
        let fin_sample = None;

        let m_in: usize = io_genot::compute_num_snv(&fin, gfmt).unwrap();
        println!("{}", m_in);
        let n_in: usize = io_genot::compute_num_sample(&fin, gfmt).unwrap();
        println!("{}", n_in);
        // load snvs
        let snvs_in = io_genot::load_snvs(&fin, gfmt);
        let (m, use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (n, use_samples) = sample::make_use_samples(fin_sample, &fin, gfmt);
        let sample_id_to_n = samples::create_sample_id_to_n(&fin, gfmt, Some(&use_samples));
        let ys: Vec<bool> = io_genot::load_ys_buf(&fin, gfmt, None, None, &sample_id_to_n);

        (fin, ys, m, n, use_snvs, use_samples)
    }

    // mv to genot_struct.rs
    //#[test]
    //fn test_missing_to_mode() {
    //    let mut x = vec![0, 1, 1, 2, 1, 3];
    //    let x_exp = vec![0, 1, 1, 2, 1, 1];

    //    missing_to_mode(&mut x);

    //    for (x_v, x_exp_v) in x.iter().zip(x_exp.iter()) {
    //        assert_eq!(*x_v, *x_exp_v);
    //    }

    //    // when missing is the mode
    //    let mut x = vec![0, 1, 1, 3, 3, 3];
    //    let x_exp = vec![0, 1, 1, 1, 1, 1];

    //    missing_to_mode(&mut x);

    //    for (x_v, x_exp_v) in x.iter().zip(x_exp.iter()) {
    //        assert_eq!(*x_v, *x_exp_v);
    //    }
    //}
}
