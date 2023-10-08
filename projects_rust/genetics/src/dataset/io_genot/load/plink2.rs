use crate::genot::prelude::*;
use crate::genot_index;
use crate::{alloc, GenotFormat};
use crate::{io_genot, vec, Chrom};
use pgenlib;
use rayon::prelude::*;
use std::ffi::CString;
use std::os::unix::prelude::OsStrExt;
use std::path::Path;
use std::time::Instant;

fn fgenot_pgenlib(fin: &Path, gfmt: GenotFormat, chrom: Option<&Chrom>) -> CString {
    let fin_genot = io_genot::fname_plinks_genot(fin, gfmt, chrom);
    let fin_str = CString::new(fin_genot.as_os_str().as_bytes()).unwrap();
    fin_str
}

fn convert_use_samples_pgenlib(use_samples: Option<&[bool]>, n: usize) -> Vec<i32> {
    let use_samples_idx: Vec<i32> = match use_samples {
        None => Vec::from_iter(0..(n as i32)),
        Some(usev) => usev
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(i, _)| i.try_into().unwrap())
            .collect(),
    };

    use_samples_idx
}

// TODO: chrom
pub fn generate_genot_snv_plink2(
    fin: &Path,
    gfmt: GenotFormat,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    //use_missing: bool,
    fill_missing: bool,
) -> GenotSnv {
    let n_in = match use_samples {
        None => n,
        Some(usev) => usev.len(),
    };

    let use_samples_idx = convert_use_samples_pgenlib(use_samples, n);
    assert_eq!(use_samples_idx.len(), n);

    let fin_genot = fgenot_pgenlib(fin, gfmt, None);

    let genot_v = load_genot_snv_buf(fin_genot, mi, n_in, &use_samples_idx, n);
    //println!("genot_v {:?}", genot_v);

    let mut g_snv = GenotSnv::new_empty(n);
    assign_pred_from_genot_i8(&mut g_snv.as_genot_snv_mut_snv(), &genot_v);

    //if !use_missing {
    if fill_missing {
        //super::fill_missing_snv(&mut g_snv.as_genot_snv_mut_snv());
        g_snv.as_genot_snv_mut_snv().fill_missing_mode();
    }

    g_snv
}

fn load_genot_snv_buf(
    fin_genot: CString,
    mi: usize,
    n_in: usize,
    use_samples_idx: &[i32],
    //mut use_samples_idx: Vec<i32>,
    n: usize,
    //) -> Vec<f64> {
) -> Vec<i8> {
    let m_start = mi;
    let m_end = mi + 1;
    let use_snvs = vec![true; 1];

    let mut genot_v = vec![0i8; n];

    load_genot_snvs_extract_buf(
        fin_genot,
        m_start,
        m_end,
        &use_snvs,
        n_in,
        use_samples_idx,
        n,
        &mut genot_v,
    );
    genot_v

    /*     let nthr = rayon::current_num_threads();
    let mut genot_v = vec![0.0f64; n];
    unsafe {
        let _ = pgenlib::pgenreader_load_snv(
            genot_v.as_mut_ptr(),
            fin_genot.as_ptr(),
            mi.try_into().unwrap(),
            n_in.try_into().unwrap(),
            use_samples_idx.as_mut_ptr(),
            n.try_into().unwrap(),
            nthr.try_into().unwrap(),
        );
    }
    genot_v */
}

/* // {0:0, 1:1, 2:2, -3:3}
pub fn assign_pred_from_genot(pred: &mut GenotSnvMut, buf_mi: &[f64]) {
    for (ni, dosage) in buf_mi.iter().enumerate() {
        if *dosage < 0.0 {
            // missing
            pred.set_unchecked(3, ni);
        } else {
            let d = *dosage as u8;
            pred.set_unchecked(d, ni);
        }
    }
} */

pub fn assign_pred_from_genot_i8(pred: &mut GenotSnvMut, buf_mi: &[i8]) {
    for (ni, dosage) in buf_mi.iter().enumerate() {
        if *dosage < 0i8 {
            // missing
            pred.set_unchecked(3, ni);
        } else {
            let d = *dosage as u8;
            pred.set_unchecked(d, ni);
        }
    }
}

// load whole is fastest
// TODO: untest for split_chrom
pub fn generate_genot_plink2(
    fin: &Path,
    gfmt: GenotFormat,
    m: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    //use_missing: bool,
    fill_missing: bool,
) -> Genot {
    log::debug!("to prepare Genot plink2 m, n: {}, {}", m, n);
    let start = Instant::now();

    let mem = alloc::get_available_memory();
    log::debug!("available mem: {:?} bytes", mem);
    //log::info!("Temporary skip panic even for insufficient memory");

    let genot_byte = Genot::byte(m, n);
    // panic if available mem < genot_byte + mem_pgen_min
    let mem_pgenlib_min_check = 1usize * 1024 * 1024 * 1024;
    //let mem_pgenlib_min = 16usize * 1024 * 1024 * 1024;
    //log::info!("Temporary skip panic even for insufficient memory");
    match mem {
        Some(x) => {
            log::debug!(
                "genot + min pgenlib vs available mem, {:.3} GB + {:.3} GB vs {:.3} GB",
                alloc::mem_gb(genot_byte),
                alloc::mem_gb(mem_pgenlib_min_check),
                alloc::mem_gb(x),
            );
            if genot_byte + mem_pgenlib_min_check > x {
                panic!("Memory insufficient on preparing Genot.")
            }
        }
        None => {
            log::debug!("Could not get available memory.");
        }
    };

    // TODO: better way
    let use_snvs_v = vec![true; m];
    let use_snvs = match use_snvs {
        Some(x) => x,
        None => &use_snvs_v,
    };
    let m_in = use_snvs.len();

    let n_in = match use_samples {
        None => n,
        Some(usev) => usev.len(),
    };

    let use_samples_idx = convert_use_samples_pgenlib(use_samples, n);
    assert_eq!(use_samples_idx.len(), n);

    let mut g = Genot::new_zeros(m, n);
    let is_split_chrom = io_genot::judge_split_chrom(fin);
    if is_split_chrom {
        let mut m_begin = 0;
        let mut m_in_begin = 0;
        for chrom_i in Chrom::variants().iter() {
            log::debug!("Loading chromosome {}", chrom_i);
            let m_in_chrom = io_genot::compute_num_snv_file_chrom(fin, gfmt, Some(chrom_i));
            if m_in_chrom.is_none() {
                continue;
            }
            let m_in_chrom = m_in_chrom.unwrap();
            let m_in_end = m_in_begin + m_in_chrom;
            // count non-zero
            let m_chrom = vec::count_true(&use_snvs[m_in_begin..m_in_end]);
            let m_end = m_begin + m_chrom;
            if m_chrom == 0 {
                m_in_begin = m_in_end;
                continue;
            }
            let fin_genot = fgenot_pgenlib(fin, gfmt, Some(chrom_i));

            assign_genot(
                &mut g.as_genot_snvs_mut(m_begin, m_end),
                fin_genot,
                m_in,
                n_in,
                &use_samples_idx,
                //use_samples_idx.clone(),
                n,
                &use_snvs[m_in_begin..m_in_end],
            );

            /*             let buf = load_genot_whole_buf(fin_genot, m_in_chrom, n_in, use_samples_idx.clone(), n);
            assign_genot_buf(
                &mut g.as_genot_snvs_mut(m_begin, m_end),
                buf,
                &use_snvs[m_in_begin..m_in_end],
            ); */
            m_begin = m_end;
            m_in_begin = m_in_end;
        }
        assert_eq!(m_in_begin, use_snvs.len(), "Sth wrong.");
    } else {
        let fin_genot = fgenot_pgenlib(fin, gfmt, None);
        assign_genot(
            &mut g.as_genot_mut(),
            fin_genot,
            m_in,
            n_in,
            &use_samples_idx,
            n,
            use_snvs,
        );
    }

    // missing
    //if !use_missing {
    if fill_missing {
        g.iter_snv_mut()
            .par_bridge()
            .for_each(|mut g_snv| g_snv.fill_missing_mode());
        //.for_each(|mut g_snv| super::fill_missing_snv(&mut g_snv));
    }

    let end = start.elapsed();
    log::info!("It took {} seconds to generate genot.", end.as_secs());

    g
}

// max size allocated in pgenlib
const BUF_SIZE_PGENLIB_LIMIT: usize = 64 * 1024 * 1024 * 1024;

fn assign_genot(
    g_chrom: &mut GenotMut,
    fin_genot: CString,
    m_in_chrom: usize,
    n_in: usize,
    use_samples_idx: &[i32],
    n: usize,
    use_snvs: &[bool],
    // TODO: buf: Option<Vec<i32>> // for chrom
) {
    //let buf_size_limit: usize = BUF_SIZE_PGENLIB_LIMIT;

    let genot_byte = Genot::byte(vec::count_true(use_snvs), n);
    let buf_size_limit: usize = match alloc::get_available_memory() {
        Some(x) => x.min(BUF_SIZE_PGENLIB_LIMIT),
        None => {
            log::debug!(
                "Could not get available memory; assume there is {} GB available memory.",
                alloc::mem_gb(BUF_SIZE_PGENLIB_LIMIT)
            );
            BUF_SIZE_PGENLIB_LIMIT
        }
    };
    log::debug!("buf_size_limit: {} GB", alloc::mem_gb(buf_size_limit));

    // 1 byte (i8) per count in pgenlib
    let byte_per_snv = n * 1;
    // f64: 8 byte per count in pgenlib
    // let byte_per_snv = n * 8;
    let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
    let buf_num_snv: usize = buf_num_snv_limit.min(m_in_chrom);
    //let buf_size: usize = buf_num_snv * byte_per_snv;
    //assert_eq!(buf_size % byte_per_snv, 0);
    //assert!(buf_size <= buf_size_limit);

    // TMP
    //let buf_num_snv = 10;

    let mut buf = vec![0i8; buf_num_snv * n];
    //let mut buf = vec![0.0f64; buf_num_snv * n];

    let mut m_in_begin_loaded = 0;
    let mut m_begin_loaded = 0;
    loop {
        log::debug!("m_in_begin_loaded: {}", m_in_begin_loaded);
        let m_in_read = buf_num_snv.min(m_in_chrom - m_in_begin_loaded);
        log::debug!("m_in_read: {}", m_in_read);

        let m_in_end_loaded = m_in_begin_loaded + m_in_read;
        let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];
        let (_, m_read) = genot_index::create_m_to_m_in(use_snvs_loaded);
        log::debug!("m_read: {}", m_read);
        log::debug!("m_in_end_loded: {}", m_in_end_loaded);

        let m_end_loaded = m_begin_loaded + m_read;

        if m_read != 0 {
            //let buf = load_genot_whole_buf(fin_genot, m_in, n_in, use_samples_idx, n);
            //load_genot_snvs_buf(
            load_genot_snvs_extract_buf(
                fin_genot.clone(),
                m_in_begin_loaded,
                m_in_end_loaded,
                use_snvs_loaded,
                n_in,
                use_samples_idx,
                n,
                &mut buf,
            );
            //println!("buf {:?}", &buf[..10]);

            let mut g_chrom_part = g_chrom.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);
            //assign_genot_buf(&mut g_chrom_part, &buf, use_snvs_loaded);
            assign_genot_extract_buf(&mut g_chrom_part, &buf);
        }

        m_begin_loaded = m_end_loaded;
        m_in_begin_loaded = m_in_end_loaded;
        assert!(m_in_begin_loaded <= m_in_chrom);
        if m_in_begin_loaded == m_in_chrom {
            break;
        }
    }
    assert_eq!(m_in_begin_loaded, m_in_chrom);
}

/* fn assign_genot(
    g_chrom: &mut GenotMut,
    fin_genot: CString,
    m_in_chrom: usize,
    n_in: usize,
    use_samples_idx: Vec<i32>,
    //use_samples_idx: Vec<i32>,
    n: usize,
    use_snvs: &[bool],
    // TODO: buf: Option<Vec<i32>> // for chrom
) {
    let buf_size_limit: usize = BUF_SIZE_PGENLIB_LIMIT;
    // 8 byte per count in pgenlib
    let byte_per_snv = n * 8;
    let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
    let buf_num_snv: usize = buf_num_snv_limit.min(m_in_chrom);
    //let buf_size: usize = buf_num_snv * byte_per_snv;
    //assert_eq!(buf_size % byte_per_snv, 0);
    //assert!(buf_size <= buf_size_limit);

    // TMP
    //let buf_num_snv = 10;

    let mut buf = vec![0i8; buf_num_snv * n];
    //let mut buf = vec![0.0f64; buf_num_snv * n];

    let mut m_in_begin_loaded = 0;
    let mut m_begin_loaded = 0;
    loop {
        log::debug!("m_in_begin_loaded: {}", m_in_begin_loaded);
        let m_in_read = buf_num_snv.min(m_in_chrom - m_in_begin_loaded);
        log::debug!("m_in_read: {}", m_in_read);

        let m_in_end_loaded = m_in_begin_loaded + m_in_read;
        let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];
        let (_, m_read) = genot_index::create_m_to_m_in(use_snvs_loaded);
        log::debug!("m_read: {}", m_read);
        log::debug!("m_in_end_loded: {}", m_in_end_loaded);

        let m_end_loaded = m_begin_loaded + m_read;

        if m_read != 0 {
            //let buf = load_genot_whole_buf(fin_genot, m_in, n_in, use_samples_idx, n);
            //load_genot_snvs_extract_buf(
            load_genot_snvs_buf(
                fin_genot.clone(),
                m_in_begin_loaded,
                m_in_end_loaded,
                n_in,
                //use_samples_idx.clone(),
                //&mut use_samples_idx,
                &use_samples_idx,
                n,
                &mut buf,
            );
            //println!("buf {:?}", &buf[..10]);

            let mut g_chrom_part = g_chrom.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);
            assign_genot_buf(&mut g_chrom_part, &buf, use_snvs_loaded);
        }

        m_begin_loaded = m_end_loaded;
        m_in_begin_loaded = m_in_end_loaded;
        assert!(m_in_begin_loaded <= m_in_chrom);
        if m_in_begin_loaded == m_in_chrom {
            break;
        }
    }
    assert_eq!(m_in_begin_loaded, m_in_chrom);
} */

fn load_genot_snvs_extract_buf(
    fin_genot: CString,
    m_start: usize,
    m_end: usize,
    use_snvs: &[bool],
    n_in: usize,
    use_samples_idx: &[i32],
    n: usize,
    buf: &mut Vec<i8>,
    //buf: &mut Vec<f64>,
) {
    // too large nthr might leads mem error.
    let max_threads_pgenlib = 32;
    let nthr = rayon::current_num_threads().min(max_threads_pgenlib);
    log::debug!("nthr in rust {}", nthr);

    let m_read = vec::count_true(use_snvs);

    buf.resize(m_read * n, 0i8);

    unsafe {
        let _ = pgenlib::pgenreader_load_snvs_extract(
            buf.as_mut_ptr(),
            fin_genot.as_ptr(),
            m_start.try_into().unwrap(),
            m_end.try_into().unwrap(),
            use_snvs.as_ptr(),
            n_in.try_into().unwrap(),
            use_samples_idx.as_ptr(),
            n.try_into().unwrap(),
            nthr.try_into().unwrap(),
        );
    }
    //println!("buf_tmp {:?}", &buf_tmp[..10]);
    //println!("genot_v {:?}", genot_v);
    //bufv
}

#[allow(dead_code)]
fn load_genot_snvs_buf(
    fin_genot: CString,
    m_start: usize,
    m_end: usize,
    n_in: usize,
    use_samples_idx: &Vec<i32>,
    n: usize,
    buf: &mut Vec<i8>,
) {
    let m_in = m_end - m_start;
    let use_snvs = vec![true; m_in];

    load_genot_snvs_extract_buf(
        fin_genot,
        m_start,
        m_end,
        &use_snvs,
        n_in,
        use_samples_idx,
        n,
        //&mut buf,
        buf,
    );
    //genot_v

    /*
    let nthr = rayon::current_num_threads();
    //let nthr = 4;
    log::debug!("nthr in rust {}", nthr);

    //log::debug!("buf size {}",(m_end - m_start) * n);
    buf.resize((m_end - m_start) * n, 0.0f64);

    // TODO: convert f64->i8 in pgenlib
    unsafe {
        let _ = pgenlib::pgenreader_load_snvs(
            buf.as_mut_ptr(),
            fin_genot.as_ptr(),
            m_start.try_into().unwrap(),
            m_end.try_into().unwrap(),
            n_in.try_into().unwrap(),
            use_samples_idx.as_ptr(),
            //use_samples_idx.as_mut_ptr(),
            n.try_into().unwrap(),
            nthr.try_into().unwrap(),
        );
    }
    //println!("genot_v {:?}", genot_v);
    //bufv */
}

/// won't work for 1M SNVs x 300K samples for f64
/// not tried for i8
#[allow(dead_code)]
fn load_genot_whole_buf(
    fin_genot: CString,
    m_in: usize,
    n_in: usize,
    use_samples_idx: &[i32],
    n: usize,
) -> Vec<i8> {
    let m_start = 0;
    let m_end = m_in;
    let use_snvs = vec![true; m_in];

    let mut genot_v = vec![0i8; m_in * n];

    load_genot_snvs_extract_buf(
        fin_genot,
        m_start,
        m_end,
        &use_snvs,
        n_in,
        use_samples_idx,
        n,
        &mut genot_v,
    );
    genot_v

    /*     //let mut genot_v = vec![0.0f64; m_in * n];
    unsafe {
        let _ = pgenlib::pgenreader_load_whole(
            genot_v.as_mut_ptr(),
            fin_genot.as_ptr(),
            m_in.try_into().unwrap(),
            n_in.try_into().unwrap(),
            use_samples_idx.as_mut_ptr(),
            n.try_into().unwrap(),
            nthr.try_into().unwrap(),
        );
    }
    //println!("genot_v {:?}", genot_v);
    genot_v */
}

fn assign_genot_extract_buf(g: &mut GenotMut, buf: &[i8]) {
    //let (m_to_m_in, m_read) = genot_index::create_m_to_m_in(use_snvs);

    //assert
    //assert_eq!(g.m(), m_read);
    let n = g.n();

    g.iter_snv_mut()
        .enumerate()
        .par_bridge()
        .for_each(|(mi, mut g_snv)| {
            //let m_in_i = m_to_m_in[&mi];
            //let buf_mi = &buf[m_in_i * n..(m_in_i + 1) * n];
            let buf_mi = &buf[mi * n..(mi + 1) * n];

            //println!("buf_mi {:?}", &buf_mi[..10]);
            assign_pred_from_genot_i8(&mut g_snv, &buf_mi);
        });
}

//fn assign_genot_buf(g: &mut GenotMut, buf: &[f64], use_snvs: &[bool]) {
#[allow(dead_code)]
fn assign_genot_buf(g: &mut GenotMut, buf: &[i8], use_snvs: &[bool]) {
    let (m_to_m_in, m_read) = genot_index::create_m_to_m_in(use_snvs);

    //assert
    assert_eq!(g.m(), m_read);
    let n = g.n();

    g.iter_snv_mut()
        .enumerate()
        .par_bridge()
        .for_each(|(mi, mut g_snv)| {
            let m_in_i = m_to_m_in[&mi];
            let buf_mi = &buf[m_in_i * n..(m_in_i + 1) * n];

            //println!("buf_mi {:?}", &buf_mi[..10]);
            assign_pred_from_genot_i8(&mut g_snv, &buf_mi);
        });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{sample, snv};
    use std::path::PathBuf;

    fn setup_test3() -> (
        PathBuf,
        GenotFormat,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2Vzs;
        let fin_snv = None;
        let fin_sample = None;

        let m_in: usize = io_genot::compute_num_snv(&fin, gfmt).unwrap();
        log::debug!("{}", m_in);
        let n_in: usize = io_genot::compute_num_sample(&fin, gfmt).unwrap();
        log::debug!("{}", n_in);
        // load snvs
        let snvs_in = io_genot::load_snvs(&fin, gfmt);
        let (m, use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
        // should do this in snv::test
        assert_eq!(m, 3);
        assert_eq!(use_snvs.len(), 3);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (n, use_samples) = sample::make_use_samples(fin_sample, &fin, gfmt);
        assert_eq!(n, 10);
        assert_eq!(use_samples.len(), 10);
        let ys = vec![
            true, true, true, true, true, false, false, false, false, false,
        ];

        (fin, gfmt, ys, m, n, use_snvs, use_samples)
    }

    fn setup_test3_part() -> (
        PathBuf,
        GenotFormat,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2Vzs;

        let m_in: usize = io_genot::compute_num_snv(&fin, gfmt).unwrap();
        log::debug!("{}", m_in);
        let n_in: usize = io_genot::compute_num_sample(&fin, gfmt).unwrap();
        log::debug!("{}", n_in);
        let m = 2;
        let use_snvs = vec![true, false, true];
        let n = 5;
        let use_samples = vec![
            false, true, false, true, false, true, false, true, false, true,
        ];
        let ys = vec![
            true, true, true, true, true, false, false, false, false, false,
        ];

        (fin, gfmt, ys, m, n, use_snvs, use_samples)
    }

    fn setup_test3ref() -> (
        PathBuf,
        GenotFormat,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy3/genot.ref");
        let gfmt = GenotFormat::Plink2Vzs;
        let fin_snv = None;
        let fin_sample = None;

        let m_in: usize = io_genot::compute_num_snv(&fin, gfmt).unwrap();
        log::debug!("{}", m_in);
        let n_in: usize = io_genot::compute_num_sample(&fin, gfmt).unwrap();
        log::debug!("{}", n_in);
        // load snvs
        let snvs_in = io_genot::load_snvs(&fin, gfmt);
        let (m, use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
        // should do this in snv::test
        assert_eq!(m, 3);
        assert_eq!(use_snvs.len(), 3);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (n, use_samples) = sample::make_use_samples(fin_sample, &fin, gfmt);
        assert_eq!(n, 10);
        assert_eq!(use_samples.len(), 10);
        let ys = vec![
            true, true, true, true, true, false, false, false, false, false,
        ];

        (fin, gfmt, ys, m, n, use_snvs, use_samples)
    }

    fn setup_test3_plink2() -> (
        PathBuf,
        GenotFormat,
        Vec<bool>,
        usize,
        usize,
        Vec<bool>,
        Vec<bool>,
    ) {
        let fin = PathBuf::from("../../test/data/toy3/genot");
        let gfmt = GenotFormat::Plink2;
        let fin_snv = None;
        let fin_sample = None;

        let m_in: usize = io_genot::compute_num_snv(&fin, gfmt).unwrap();
        log::debug!("{}", m_in);
        let n_in: usize = io_genot::compute_num_sample(&fin, gfmt).unwrap();
        log::debug!("{}", n_in);
        // load snvs
        let snvs_in = io_genot::load_snvs(&fin, gfmt);
        let (m, use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
        assert_eq!(m, 3);
        assert_eq!(use_snvs.len(), 3);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (n, use_samples) = sample::make_use_samples(fin_sample, &fin, gfmt);
        assert_eq!(n, 10);
        assert_eq!(use_samples.len(), 10);
        //let ys: Vec<bool> = io_genot::load_ys(&fin, gfmt, None, None, &use_samples);
        let ys = vec![
            true, true, true, true, true, false, false, false, false, false,
        ];

        (fin, gfmt, ys, m, n, use_snvs, use_samples)
    }

    #[test]
    fn test_assign_pred_from_bed() {
        let mut g = GenotSnv::new_empty(6);
        // [2, 0, 3, 0, 1, 0]
        //let pbuf = vec![2.0f64, 0.0, 3.0, 0.0, 1.0, 0.0];
        let pbuf = vec![2i8, 0, 3, 0, 1, 0];

        assign_pred_from_genot_i8(&mut g.as_genot_snv_mut_snv(), &pbuf);
        assert_eq!(g.vals(), vec![2u8, 0, 3, 0, 1, 0]);
    }

    #[test]
    fn test_generate_genot_snv_plink2() {
        let (fin, gfmt, _, _, n, _, use_samples) = setup_test3_plink2();
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        //let g = generate_genot_plink(&fin, gfmt, m, n, &use_snvs, Some(&use_samples), true);
        let mi = 2;
        let g = generate_genot_snv_plink2(&fin, gfmt, mi, n, Some(&use_samples), false);
        assert_eq!(g.vals(), vec![2, 0, 1, 0, 1, 2, 0, 1, 0, 3]);
    }

    #[test]
    fn test_generate_genot_snv_plink2vzs() {
        let (fin, gfmt, _, _, n, _, use_samples) = setup_test3();
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        //let g = generate_genot_plink(&fin, gfmt, m, n, &use_snvs, Some(&use_samples), true);
        let mi = 2;
        let g = generate_genot_snv_plink2(&fin, gfmt, mi, n, Some(&use_samples), false);
        assert_eq!(g.vals(), vec![2, 0, 1, 0, 1, 2, 0, 1, 0, 3]);
    }

    #[test]
    fn test_generate_genot_plink2vzs_ref() {
        let (fin, gfmt, _, _, n, _, use_samples) = setup_test3ref();
        let mi = 2;
        let g = generate_genot_snv_plink2(&fin, gfmt, mi, n, Some(&use_samples), false);
        assert_eq!(g.vals(), vec![0, 2, 1, 2, 1, 0, 2, 1, 2, 3]);
    }

    #[test]
    fn test_generate_genot_snv_plink2vzs_part() {
        let (fin, gfmt, _, _, n, _, use_samples) = setup_test3_part();
        let mi = 2;
        let g = generate_genot_snv_plink2(&fin, gfmt, mi, n, Some(&use_samples),false);
        assert_eq!(g.vals(), vec![0, 0, 2, 1, 3]);
    }

    #[test]
    fn test_generate_genot_whole_plink2vzs_part() {
        let (fin, gfmt, _, m, n, use_snvs, use_samples) = setup_test3_part();
        let g = generate_genot_plink2(&fin, gfmt, m, n, Some(&use_snvs), Some(&use_samples), false);
        let mut giter = g.iter_snv();
        // [2,1,0,0,0,2,1,0,0,3]
        assert_eq!(giter.next().unwrap().vals(), vec![1, 0, 2, 0, 3]);
        assert_eq!(giter.next().unwrap().vals(), vec![0, 0, 2, 1, 3]);
        assert_eq!(giter.next(), None);
    }

    #[test]
    fn test_generate_genot_whole_plink2vzs() {
        let (fin, gfmt, _, m, n, use_snvs, use_samples) = setup_test3();
        let g = generate_genot_plink2(&fin, gfmt, m, n, Some(&use_snvs), Some(&use_samples), false);
        let mut giter = g.iter_snv();
        // [2,1,0,0,0,2,1,0,0,3]
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![2, 1, 0, 0, 0, 2, 1, 0, 0, 3]
        );
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![2, 0, 1, 1, 0, 1, 0, 2, 0, 3]
        );
        assert_eq!(
            giter.next().unwrap().vals(),
            vec![2, 0, 1, 0, 1, 2, 0, 1, 0, 3]
        );
        assert_eq!(giter.next(), None);
    }
}
