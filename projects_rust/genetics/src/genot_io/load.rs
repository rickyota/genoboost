pub mod plink;
#[cfg(feature = "plink2")]
pub mod plink2;

use std::collections::HashMap;

use crate::genot::prelude::*;
use crate::{genot_io, vec, GenotFile};

/*
// This is not fast for plink2
//  unused now: for small number of snvs called from python
//  mi for multiple files
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
*/

/// Generate Vector of predictions.
pub fn generate_genot_simple(
    fin_genot: &GenotFile,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    fill_missing_mode: bool,
    m: Option<usize>,
    n: Option<usize>,
    mem: Option<usize>,
) -> Genot {
    let m = match m {
        Some(x) => x,
        None => match use_snvs {
            Some(x) => vec::count_true(x),
            None => genot_io::compute_num_snv(&fin_genot).unwrap(),
        },
    };

    let n = match n {
        Some(x) => x,
        None => match use_samples {
            Some(x) => vec::count_true(x),
            None => genot_io::compute_num_sample(fin_genot).unwrap(),
        },
    };

    generate_genot(
        fin_genot,
        m,
        n,
        use_snvs,
        use_samples,
        fill_missing_mode,
        mem,
        None,
    )
}

pub fn generate_genot(
    fin_genot: &GenotFile,
    m: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    fill_missing_mode: bool,
    mem: Option<usize>,
    group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
) -> Genot {
    match fin_genot {
        GenotFile::Plink1(_) => plink::generate_genot_plink(
            fin_genot,
            m,
            n,
            use_snvs,
            use_samples,
            fill_missing_mode,
            mem,
            group_to_m_in,
        ),
        GenotFile::Plink2(_) | GenotFile::Plink2Vzs(_) => {
            if group_to_m_in.is_some() {
                unimplemented!("group snv is not implemented for plink2")
            }
            call_generate_genot_plink2(
                fin_genot,
                m,
                n,
                use_snvs,
                use_samples,
                fill_missing_mode,
                mem,
            )
        }
    }
}

/* // error: plink2 is not loaded without feature=plink2
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
    fin_genot: &GenotFile,
    m: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    fill_missing_mode: bool,
    mem: Option<usize>,
) -> Genot {
    plink2::generate_genot_plink2(
        fin_genot,
        m,
        n,
        use_snvs,
        use_samples,
        fill_missing_mode,
        mem,
    )
}

#[cfg(not(feature = "plink2"))]
fn call_generate_genot_plink2(
    _: &GenotFile,
    _: usize,
    _: usize,
    _: Option<&[bool]>,
    _: Option<&[bool]>,
    _: bool,
    _: Option<usize>,
) -> Genot {
    panic!("Cannot use plink2 in this program feature. Use --feature plink2");
}
