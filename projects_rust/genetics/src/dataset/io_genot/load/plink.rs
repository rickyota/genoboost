use crate::genot::prelude::*;
use crate::{alloc, GenotFormat};
use crate::{genot_index, textfile};
use crate::{io_genot, vec, Chrom};
//use cmatrix::BaseCMatrix;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufRead;
use std::io::{prelude::*, BufReader, SeekFrom};
use std::path::Path;

pub fn generate_genot_snv_plink(
    fin: &Path,
    gfmt: GenotFormat,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    //use_missing: bool,
    fill_missing: bool,
) -> GenotSnv {
    let reader = BufReader::new(File::open(io_genot::fname_plinks_genot(fin, gfmt, None)).unwrap());

    let mut g_snv = GenotSnv::new_empty(n);

    load_snv_read(
        &mut g_snv.as_genot_snv_mut_snv(),
        reader,
        mi,
        use_samples,
        n,
    );

    //if !use_missing {
    if fill_missing {
        //super::fill_missing_snv(&mut g_snv.as_genot_snv_mut_snv());
        g_snv.as_genot_snv_mut_snv().fill_missing_mode()
    }

    g_snv
}

fn load_snv_read<R: BufRead + Seek>(
    g_snv: &mut GenotSnvMut,
    mut reader: R,
    mi: usize,
    use_samples: Option<&[bool]>,
    n: usize,
) {
    let n_in = if let Some(v) = use_samples {
        v.len()
    } else {
        n
    };
    assert!(n <= n_in);

    // load to buf
    let byte_per_snv = io_genot::bed_per_snv_size(n_in);
    let buf_size = byte_per_snv;
    let mut buf: Vec<B8_2> = vec![0; buf_size];
    let buf_begin = 3 + mi * byte_per_snv;
    reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
    let loaded_byte = reader.read(&mut buf).unwrap();
    assert_eq!(loaded_byte, byte_per_snv);

    assign_pred_from_bed(g_snv, &buf, use_samples);
}

const BUF_SIZE_BED_LIMIT: usize = 64 * 1024 * 1024 * 1024;
fn load_buf_size_max(fin: &Path, gfmt: GenotFormat) -> usize {
    let mut buf_size_max = 0usize;
    for chrom_i in Chrom::variants().iter() {
        let file_size =
            textfile::file_size(&io_genot::fname_plinks_genot(&fin, gfmt, Some(chrom_i)));
        if let Some(x) = file_size {
            buf_size_max = buf_size_max.max(x);
        }
    }
    buf_size_max
}

/// Sequentially load part of bed and convert into predictions since bed could be too large to load on mem.
/// predictions: m*2*n
// load whole is fastest
// TODO: use_missing -> !fill_missing
pub fn generate_genot_plink(
    fin: &Path,
    gfmt: GenotFormat,
    m: usize,
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    //use_missing: bool,
    fill_missing: bool,
) -> Genot {
    log::debug!("to prepare Genot m, n: {}, {}", m, n);

    let mem = alloc::get_available_memory();
    log::debug!("available mem: {:?} bytes", mem);
    let genot_byte = Genot::byte(m, n);
    //log::info!("Temporary skip panic even for insufficient memory");
    match mem {
        Some(x) => {
            log::debug!(
                "genot vs available mem, {:.3} bytes vs {:.3} bytes",
                genot_byte,
                x
            );
            if genot_byte > x {
                panic!("Memory insufficient on preparing Genot.")
            }
        }
        None => {
            log::info!("Could not get available memory.");
        }
    };

    // for use_snvs=None
    let use_snvs_v = vec![true; m];
    let use_snvs = match use_snvs {
        Some(x) => x,
        None => &use_snvs_v,
    };
    //let m_in = use_snvs.len();

    let mut g = Genot::new_zeros(m, n);
    //log::debug!("done preparing Genot");

    // Two patterns to assign predictions
    // bed file is split into chrom or not
    let is_split_chrom = io_genot::judge_split_chrom(fin);

    if is_split_chrom {
        let mut m_begin = 0;
        let mut m_in_begin = 0;
        //let buf_size_limit: usize = BUF_SIZE_BED_LIMIT;
        // TODO: somehow same value as above; why not decreased by alloc in genot?

        // FIXME: fix mem for pg a*
        //log::info!("Temporary fix mem to 64GB.");
        //let mem = Some(64usize * 1024 * 1024 * 1024);
        let mem = alloc::get_available_memory();
        log::debug!("available mem: {:?} bytes", mem);

        let buf_size_limit = mem.map_or_else(|| BUF_SIZE_BED_LIMIT, |x| x.min(BUF_SIZE_BED_LIMIT));
        log::debug!("buf_size_limit: {:?} bytes", buf_size_limit);

        let buf_size_max = load_buf_size_max(fin, gfmt).min(buf_size_limit);
        log::debug!("buf_size_max {}", buf_size_max);
        let mut buf: Vec<B8_2> = vec![0; buf_size_max];
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

            // TODO
            //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

            // using  File.read() did not speed up
            //let fin_chrom = plink::fname_chrom(fin, Some(chrom_i));
            // 512 KB
            let reader_cap = 512usize * 1024;
            let reader = BufReader::with_capacity(
                reader_cap,
                File::open(io_genot::fname_plinks_genot(&fin, gfmt, Some(chrom_i))).unwrap(),
            );
            //BufReader::new(File::open(plink::fname_bed(&fin_chrom, None)).unwrap());

            // reuse buf
            buf = assign_genot(
                &mut g.as_genot_snvs_mut(m_begin, m_end),
                reader,
                &use_snvs[m_in_begin..m_in_end],
                use_samples,
                Some(buf),
            );
            m_begin = m_end;
            m_in_begin = m_in_end;
        }
        assert_eq!(m_in_begin, use_snvs.len(), "Sth wrong.");
    } else {
        let reader =
            BufReader::new(File::open(io_genot::fname_plinks_genot(&fin, gfmt, None)).unwrap());
        assign_genot(&mut g.as_genot_mut(), reader, use_snvs, use_samples, None);
    }

    // missing
    //if !use_missing {
    if fill_missing {
        g.iter_snv_mut()
            .par_bridge()
            .for_each(|mut g_snv| g_snv.fill_missing_mode());
        //.for_each(|mut g_snv| super::fill_missing_snv(&mut g_snv));
    }

    g
}

/*
// directly use assign_predictions() applicable to both
fn assign_predictions_toggle<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
) {
    let file_size = reader.seek(SeekFrom::End(0)).unwrap() as usize;
    log::debug!("seek {}", file_size);

    log::debug!("Use whole strategy to load.");
    assign_predictions_whole(g_chrom, reader, use_snvs, use_samples);
    //log::debug!("Use chunk strategy to load.");
    //assign_predictions_chunk(g_chrom, reader, use_snvs, use_samples);

    //// 64 GB
    //let thres_size = 64 * 1024 * 1024 * 1024;
    //// whole have to allocate memory same size as file
    //if file_size > thres_size {
    //    log::debug!("Use chunk strategy to load.");
    //    assign_predictions_chunk(g_chrom, reader, use_snvs, use_samples);
    //} else {
    //    log::debug!("Use whole strategy to load.");
    //    assign_predictions_whole(g_chrom, reader, use_snvs, use_samples);
    //}
} */

// TODO: use SIMD
// can I implement buf: Option<&mut Vec<u8>> ?
// now unknown way to declare Vec<u8> with longer lifetime when buf=None
// What happens if no var is loaded in each loop? can I shorten time?
/// available for both whole and chunk
fn assign_genot<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
    buf: Option<Vec<B8_2>>,
) -> Vec<B8_2> {
    let m_in_chrom = use_snvs.len();

    let n_in = if let Some(v) = use_samples {
        v.len()
    } else {
        g_chrom.n()
    };

    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    // FIXME: available mem?
    // check if 32 GB? 16GB? remains or not
    // reading smaller than buf_size is not error but reason is unknown
    // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
    // -> to deal with this, every time use Seek to change position
    let buf_size_limit: usize = BUF_SIZE_BED_LIMIT;
    //let buf_size_limit: usize = 64 * 1024 * 1024 * 1024;
    //let buf_size_limit: usize = 4 * 1024 * 1024 * 1024;
    let byte_per_snv = io_genot::bed_per_snv_size(n_in);
    let but_num_snv_limit: usize = buf_size_limit / byte_per_snv;
    let buf_num_snv: usize = but_num_snv_limit.min(m_in_chrom);
    let buf_size: usize = buf_num_snv * byte_per_snv;
    assert_eq!(buf_size % byte_per_snv, 0);
    assert!(buf_size <= buf_size_limit);

    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    // TODO: check buf_size < available

    let mut buf: Vec<B8_2> = match buf {
        Some(v) => v,
        None => vec![0; buf_size],
    };
    //let mut buf: Vec<B8_2> = vec![0; buf_size];

    let mut m_in_begin_loaded = 0;
    let mut m_begin_loaded = 0;
    // read buf length in one loop
    loop {
        //log::debug!(
        //    "m_in_begin, m_begin {},{}",
        //    m_in_begin_loaded, m_begin_loaded
        //);

        // next start is m_in_begin_load
        let buf_next_begin = 3 + m_in_begin_loaded * byte_per_snv;
        //log::debug!("buf_next_begin {}", buf_next_begin);
        reader.seek(SeekFrom::Start(buf_next_begin as u64)).unwrap();
        // https://stackoverflow.com/questions/37079342/what-is-the-most-efficient-way-to-read-a-large-file-in-chunks-without-loading-th

        let m_in_read = buf_num_snv.min(m_in_chrom - m_in_begin_loaded);
        log::debug!("m_in_read: {}", m_in_read);

        let m_in_end_loaded = m_in_begin_loaded + m_in_read;
        let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];

        let (m_to_m_in, m_read) = genot_index::create_m_to_m_in(use_snvs_loaded);

        let m_end_loaded = m_begin_loaded + m_read;

        if m_read != 0 {
            // use fill() might be fast?
            // This might be smaller than buf
            let buf_size_ = m_in_read * byte_per_snv;
            buf.resize(buf_size_, 0);
            //unsafe {
            //    buf.set_len(buf_size_);
            //}
            log::debug!("buf_size_: {}", buf_size_);
            assert_eq!(buf_size_ % byte_per_snv, 0);
            assert_eq!(buf.len(), buf_size_);

            // read_exact??
            reader.read_exact(&mut buf).unwrap();
            assert_eq!(buf.len(), buf_size_);
            log::debug!("loaded: {}", buf.len());

            // x plan 1 create vec of pred_use and  pred_use.iter()...
            //  -> unsafe when pred_use is duplicated
            // x plan 2 predictions.par_chunks_mut().filter(use)
            // o plan 3 predictions_chrom[m_buf..m_buf_end].par_chunks_mut()

            //log::debug!("m_in_end, m_end {},{}", m_in_end_loaded, m_end_loaded);

            let mut g_chrom_part = g_chrom.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);

            g_chrom_part
                .iter_snv_mut()
                .enumerate()
                .par_bridge()
                .for_each(|(mi_loaded, mut g_snv)| {
                    let m_in_read_i = m_to_m_in[&mi_loaded];
                    let buf_mi = &buf[m_in_read_i * byte_per_snv..(m_in_read_i + 1) * byte_per_snv];

                    assign_pred_from_bed(&mut g_snv, &buf_mi, use_samples);
                });
        }

        m_begin_loaded = m_end_loaded;
        m_in_begin_loaded = m_in_end_loaded;
        assert!(m_in_begin_loaded <= m_in_chrom);
        if m_in_begin_loaded == m_in_chrom {
            break;
        }
    }
    assert_eq!(m_in_begin_loaded, m_in_chrom);
    buf
}

/* // TODO: use SIMD
fn assign_predictions_whole<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
) {
    //log::debug!("g_chrom len: {}", g_chrom.len());

    let n_in = if let Some(v) = use_samples {
        v.len()
    } else {
        g_chrom.n()
    };

    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    // seek is fast
    let file_size = reader.seek(SeekFrom::End(0)).unwrap() as usize;
    log::debug!("seek {}", file_size);
    let buf_size = file_size - 3;
    log::debug!("buf_size {}", buf_size);
    if buf_size == 0 {
        panic!("Buffer size is 0.");
    }

    let byte_per_snv = plink::bed_per_snv_size(n_in);
    assert_eq!(buf_size % byte_per_snv, 0);

    /*     log::debug!(
        "Loading buffer size vs available mem, {} bytes vs {} bytes",
        buf_size,
        alloc::get_available_memory()
    );
    if buf_size > alloc::get_available_memory() {
        panic!("Memory insufficient on loading bed file.");
    } */

    let mut buf: Vec<B8_2> = vec![0; buf_size];

    let current_pos = reader.seek(SeekFrom::Start(3)).unwrap();
    assert_eq!(current_pos, 3);
    reader.read_exact(&mut buf).unwrap();
    log::debug!("loaded: {}", buf.len());
    // do not use; somehow elongate buf size
    //let loaded_byte = reader.read_to_end(&mut buf).unwrap();
    //log::debug!("loaded: {} {}", loaded_byte, buf_size);
    //println!("loaded: {} {}", loaded_byte, buf_size);
    //println!("buflen {:?}", buf.len());
    //println!("buf {:?}", buf);

    //assert_eq!(loaded_byte, buf.len(), "loaded byte and buf len");
    //assert_eq!(loaded_byte, buf_size, "reader did not read whole buf size.");
    //assert_eq!(loaded_byte % byte_per_snv, 0);

    let (m_to_m_in, _) = genot_index::create_m_to_m_in(use_snvs);

    g_chrom
        .iter_snv_mut()
        .enumerate()
        .par_bridge()
        .for_each(|(mi_loaded, mut g_snv)| {
            let m_in_read_i = m_to_m_in[&mi_loaded];
            let buf_mi = &buf[m_in_read_i * byte_per_snv..(m_in_read_i + 1) * byte_per_snv];

            assign_pred_from_bed(&mut g_snv, &buf_mi, use_samples);
        });
} */

/* // TODO: use SIMD
// slower using .read() than .read_exact()
// use when file size is too large
fn assign_predictions_chunk_old<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
) {
    let m_in_chrom = use_snvs.len();

    let n_in = if let Some(v) = use_samples {
        v.len()
    } else {
        g_chrom.n()
    };

    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    // check if 32 GB? 16GB? remains or not
    // reading smaller than buf_size is not error but reason is unknown
    // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
    // -> to deal with this, every time use Seek to change position
    //let buf_size_limit: usize = 64 * 1024 * 1024 * 1024;
    let buf_size_limit: usize = 4 * 1024 * 1024 * 1024;
    let byte_per_snv = plink::bed_per_snv_size(n_in);
    let num_snv_per_buf: usize = buf_size_limit / byte_per_snv;
    let buf_num_snv: usize = num_snv_per_buf.min(m_in_chrom);
    let buf_size: usize = buf_num_snv * byte_per_snv;
    assert_eq!(buf_size % byte_per_snv, 0);
    assert!(buf_size <= buf_size_limit);

    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    // TODO: check buf_size < available
    let mut buf: Vec<B8_2> = vec![0; buf_size];

    let mut m_in_begin_loaded = 0;
    let mut m_begin_loaded = 0;
    // read buf length in one loop
    loop {
        //log::debug!(
        //    "m_in_begin, m_begin {},{}",
        //    m_in_begin_loaded, m_begin_loaded
        //);

        // next start is m_in_begin_load
        let buf_next_begin = 3 + m_in_begin_loaded * byte_per_snv;
        //log::debug!("buf_next_begin {}", buf_next_begin);
        reader.seek(SeekFrom::Start(buf_next_begin as u64)).unwrap();
        // https://stackoverflow.com/questions/37079342/what-is-the-most-efficient-way-to-read-a-large-file-in-chunks-without-loading-th
        // use fill() might be fast?
        // This might be smaller than buf
        // read_exact??
        //let loaded_byte = reader.read_exact(&mut buf).unwrap();
        let loaded_byte = reader.read(&mut buf).unwrap();
        //if loaded_byte != buf_size, last buffer or could not load whole buffer
        log::debug!("loaded: {} {}", loaded_byte, buf_size);
        //}
        if loaded_byte == 0 {
            // all buf read
            break;
        } else {
            //log::debug!("buf {:?}", buf.len());
            // might fail when only part was loaded
            //assert_eq!(loaded_byte % byte_per_snv, 0);
            let m_in_read = loaded_byte / byte_per_snv;
            // if m_in_read==0, loaded_byte < byte_per_snv, which means loaded_byte is less than one SNV.
            assert!(m_in_read > 0, "loaded_byte is less than byte of one SNV.");

            // x plan 1 create vec of pred_use and  pred_use.iter()...
            //  -> unsafe when pred_use is duplicated
            // x plan 2 predictions.par_chunks_mut().filter(use)
            // o plan 3 predictions_chrom[m_buf..m_buf_end].par_chunks_mut()

            let m_in_end_loaded = m_in_begin_loaded + m_in_read;
            let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];

            let (m_to_m_in, m_read) = genot_index::create_m_to_m_in(use_snvs_loaded);

            // here for when m_read==0
            // or move below and do not continue for m_read==0 since unnecessary
            //m_in_begin_loaded = m_in_end_loaded;

            //if m_read == 0 {
            //    continue;
            //}

            let m_end_loaded = m_begin_loaded + m_read;

            //log::debug!("m_in_end, m_end {},{}", m_in_end_loaded, m_end_loaded);

            let mut g_chrom_part = g_chrom.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);

            g_chrom_part
                .iter_snv_mut()
                .enumerate()
                .par_bridge()
                .for_each(|(mi_loaded, mut g_snv)| {
                    let m_in_read_i = m_to_m_in[&mi_loaded];
                    let buf_mi = &buf[m_in_read_i * byte_per_snv..(m_in_read_i + 1) * byte_per_snv];

                    assign_pred_from_bed(&mut g_snv, &buf_mi, use_samples);
                });

            m_begin_loaded = m_end_loaded;
            // bug! when m_read==0, will not renewed
            m_in_begin_loaded = m_in_end_loaded;
        }
    }
    assert_eq!(m_in_chrom, m_in_begin_loaded);
} */

/* // rayon impossible since cannot read inside
pub fn assign_predictions_file_loadpart(
    g_chrom: &mut GenotMut,
    fin_chrom: &str,
    use_snvs: &[bool],
    use_samples: &[bool],
) {
    let m_in_chrom = use_snvs.len();
    // log::debug!("use_samples: {:?}", use_samples);

    // assume m and m_in in m_to_m_in are sorted
    // FIXME: use_snvs should be only chrom
    //let m_to_m_in = create_m_chrom_to_m_in_chrom(use_snvs);
    //let m_in_to_m = create_m_in_to_m_chrom(use_snvs, m_begin);

    let n_in = use_samples.len();

    // for debug
    //log::debug!("ys[0]:{:#010b}", ys[0]);

    plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    let (m_to_m_in, m_for_read) = genot_index::create_m_to_m_in(use_snvs);

    let byte_per_snv = plink::bed_per_snv_size(n_in);
    let buf_size = byte_per_snv;

    let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    let mut buf: Vec<B8_2> = vec![0; buf_size];

    assert_eq!(g_chrom.m(), m_for_read, "Sth wrong");

    // rayon impossible since cannot seek
    g_chrom
        .iter_snv_mut()
        .enumerate()
        .for_each(|(m_read_i, mut g_snv)| {
            let buf_begin = 3 + m_to_m_in[&m_read_i] * byte_per_snv;
            reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
            let loaded_byte = reader.read(&mut buf).unwrap();
            assert_eq!(loaded_byte, byte_per_snv);

            //let g_snv = g_chrom.to_genot_twin_snv_mut(m_read_i);
            //let pred = predictions_snv_s_mut(predictions_chrom, m_read_i, n);
            assign_pred_from_bed(&mut g_snv, &buf, use_samples);
        });
} */

pub fn assign_pred_from_bed(pred: &mut GenotSnvMut, buf_mi: &[B8_2], use_samples: Option<&[bool]>) {
    if let Some(use_sample) = use_samples {
        let mut ni = 0;
        for (n_in_i, v) in use_sample.iter().enumerate() {
            if *v {
                // bedcode
                let bcode = buf_to_ped_code(buf_mi, n_in_i);
                pred.set_bed_code_init_unchecked(bcode, ni);
                ni += 1;
            }
        }
    } else {
        for ni in 0..pred.n() {
            // bedcode
            let bcode = buf_to_ped_code(buf_mi, ni);
            pred.set_bed_code_init_unchecked(bcode, ni);
        }
    }
}

#[inline]
pub fn buf_to_ped_code(buf: &[B8_2], ni: usize) -> u8 {
    //buf[ni // 4], ni % 3
    byte_to_ped_code(buf[ni >> 2], ni & 3)
}

#[inline]
pub fn buf_to_count(buf: &[B8_2], ni: usize) -> u8 {
    //buf[ni // 4], ni % 3
    byte_to_count(buf[ni >> 2], ni & 3)
}

// {00: 2, 01: 3, 10: 1, 11: 0}
const CODE_TO_COUNT_AR: [u8; 4] = [2, 3, 1, 0];

// plink BED code to minor allele counts
// {00: 2, 01: 3, 10: 1, 11: 0}
#[inline]
fn code_to_count(v: u8) -> u8 {
    CODE_TO_COUNT_AR[v as usize]
}

#[inline]
fn byte_to_ped_code(v: B8_2, i: usize) -> u8 {
    (v & (0x03 << (i << 1))) >> (i << 1)
}

// extract ith code (= 2i~2i+1 bits from lower)
// ex. byte_to_count(0b11011000, 2) = code_to_count(0b10) = 1
#[inline]
fn byte_to_count(v: B8_2, i: usize) -> u8 {
    //#define count_pl(c, i) (code2count((((c) & (0x03 << (i << 1))) >> (i << 1))))
    code_to_count(byte_to_ped_code(v, i))
    //code_to_count((v & (0x03 << (i << 1))) >> (i << 1))
}

// TODO: this might be clear if use .flat_map()
// no rayon here
//fn load_x(x: &mut [u8], buf_mi: &[B8_2], use_samples: &[bool]) {
//    let mut ni = 0;
//    for (n_in_i, v) in use_samples.iter().enumerate() {
//        if *v {
//            x[ni] = buf_to_count(buf_mi, n_in_i);
//            ni += 1;
//        }
//    }
//    // ng: x could be larger thant n
//    //assert_eq!(ni, x.len());
//
//    super::missing_to_mode(x);
//}

#[cfg(test)]
mod tests {
    use crate::GenotFormat;

    use super::*;
    use crate::samples;
    use crate::{sample, snv};
    use std::{io::Cursor, path::PathBuf};

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

    #[test]
    fn test_generate_predictions_whole() {
        let (fin, _, m, n, use_snvs, use_samples) = setup_test();
        let gfmt = GenotFormat::Plink1;
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        let g = generate_genot_plink(&fin, gfmt, m, n, Some(&use_snvs), Some(&use_samples), false);
        //let g = generate_genot_plink(&fin, gfmt, m, n, Some(&use_snvs), Some(&use_samples), true);
        let mut giter = g.iter_snv();
        assert_eq!(giter.next().unwrap().vals(), vec![2, 0, 1, 0, 0]);
        assert_eq!(giter.next().unwrap().vals(), vec![1, 0, 2, 1, 0]);
        assert_eq!(giter.next().unwrap().vals(), vec![0, 2, 0, 1, 1]);
    }

    /*     #[test]
    fn test_generate_predictions_whole2() {
        let (fin, _, m, n, use_snvs, use_samples) = setup_test2();
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        generate_predictions(&fin, m, n, &use_snvs, Some(&use_samples), true);
    } */

    #[test]
    fn test_assign_pred_from_bed() {
        let mut g = GenotSnv::new_empty(4);
        // [2, 0, 3, 0, 1, 0]
        let pbuf = vec![0b11_01_11_00, 0b00_00_11_10];

        let use_samples = vec![true, false, true, true, true, false];
        assign_pred_from_bed(&mut g.as_genot_snv_mut_snv(), &pbuf, Some(&use_samples));
        assert_eq!(g.vals(), vec![2, 3, 0, 1]);
    }

    #[test]
    fn test_assign_predictions() {
        let mut g = Genot::new_zeros(2, 3);

        let reader: Vec<u8> = vec![0x6c, 0x1b, 0x01, 0b00_01_11_00, 0b00_10_00_11];
        let mut cur = Cursor::new(reader.as_slice());
        assign_genot(
            &mut g.as_genot_mut(),
            &mut cur,
            //reader.as_slice(),
            &[true, true],
            None,
            None,
        );

        assert_eq!(g.vals_snv(0), vec![2, 0, 3]);
        assert_eq!(g.vals_snv(1), vec![0, 2, 1]);
    }

    #[test]
    fn test_load_snv_read() {
        let mut g = GenotSnv::new_empty(3);

        let reader: Vec<u8> = vec![0x6c, 0x1b, 0x01, 0b00_01_11_00, 0b00_10_00_11];
        let mut cur = Cursor::new(reader.as_slice());
        load_snv_read(
            &mut g.as_genot_snv_mut_snv(),
            &mut cur,
            //reader.as_slice(),
            1,
            Some(&[true, true, true]),
            3,
        );

        assert_eq!(g.vals(), vec![0, 2, 1]);
    }

    #[test]
    fn test_code_to_count() {
        let vs_input = vec![0b00, 0b01, 0b10, 0b11];
        let expects = vec![2, 3, 1, 0];
        for (v, exp) in vs_input.iter().zip(expects.iter()) {
            assert_eq!(code_to_count(*v), *exp);
        }
    }

    #[test]
    #[should_panic]
    fn test_code_to_count_panic() {
        let v = 4;
        code_to_count(v);
    }

    #[test]
    fn test_byte_to_count() {
        let v = 0b00_01_10_11;
        let expects = vec![0, 1, 3, 2];

        for (i, exp) in expects.iter().enumerate() {
            assert_eq!(byte_to_count(v, i), *exp);
        }
    }

    #[test]
    fn test_buf_to_count() {
        let buf = vec![0b00_01_10_11, 0b00_00_11_01];
        // read the lowest of the first byte
        assert_eq!(buf_to_count(&buf, 0), 0);
        // read the higest of the first byte
        assert_eq!(buf_to_count(&buf, 3), 2);
        // read the lowest of the second byte
        assert_eq!(buf_to_count(&buf, 4), 3);
        // read the highest of the second byte
        assert_eq!(buf_to_count(&buf, 7), 2);
    }
}
