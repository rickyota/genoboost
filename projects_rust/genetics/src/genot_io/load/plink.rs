use crate::genot::prelude::*;
use crate::{alloc, GenotFile};
use crate::{genot_io, vec, Chrom};
use crate::{genotype, textfile};

use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{prelude::*, BufRead, BufReader, Read, SeekFrom};
use std::path::Path;
use std::time::Instant;
//use std::io::{prelude::*, BufReader, SeekFrom};
//use std::io::{prelude::*, Read, SeekFrom}; // for seek

/*
// save
// load single snv
pub fn generate_genot_snv_plink(
    fin: &Path,
    gfmt: GenotFormat,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    fill_missing: bool,
) -> GenotSnv {
    let reader = BufReader::new(File::open(genot_io::fname_plinks_genot(fin, gfmt, None)).unwrap());

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
    let byte_per_snv = genot_io::bed_per_snv_size(n_in);
    let buf_size = byte_per_snv;
    let mut buf: Vec<B8_2> = vec![0; buf_size];
    let buf_begin = 3 + mi * byte_per_snv;
    reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
    let loaded_byte = reader.read(&mut buf).unwrap();
    assert_eq!(loaded_byte, byte_per_snv);

    assign_pred_from_bed(g_snv, &buf, use_samples);
}
*/

const BUF_SIZE_BED_MAX: usize = 64 * 1024 * 1024 * 1024;
const BUF_SIZE_BED_MIN: usize = 1 * 1024 * 1024 * 1024;

//fn load_buf_size_max(fin: &Path, gfmt: GenotFormat) -> usize {
fn load_buf_size_max(fin_genot: &GenotFile) -> usize {
    let mut buf_size_max = 0usize;
    for chrom_i in Chrom::variants().iter() {
        let file_size = textfile::file_size(&fin_genot.genotype_file(Some(chrom_i)));
        if let Some(x) = file_size {
            buf_size_max = buf_size_max.max(x);
        }
    }
    buf_size_max
}

fn mem_buf(m: usize, n: usize, mem: Option<usize>) -> usize {
    //let mem = alloc::get_available_memory();
    let mem_avail = match mem {
        Some(x) => Some(x),
        None => alloc::get_available_memory(),
    };
    log::debug!("available mem: {:?} bytes", mem_avail);
    //log::debug!("available mem: {:?} bytes", mem);
    let genot_byte = Genot::byte(m, n);
    //log::info!("Temporary skip panic even for insufficient memory");

    let mem_buf = match mem_avail {
        Some(x) => {
            log::debug!(
                "genot + min pgenlib vs available mem, {:.3} GB + {:.3} GB vs {:.3} GB",
                alloc::mem_gb(genot_byte),
                alloc::mem_gb(BUF_SIZE_BED_MIN),
                alloc::mem_gb(x),
            );
            if genot_byte + BUF_SIZE_BED_MIN > x {
                panic!("Memory insufficient on preparing Genot.")
            }
            x - genot_byte
        }
        None => {
            log::debug!(
                "Could not get available memory; assume there is {} GB available memory.",
                alloc::mem_gb(BUF_SIZE_BED_MAX)
            );
            BUF_SIZE_BED_MAX - genot_byte
        }
    };

    mem_buf
}

fn buf_read(fin_genot: &GenotFile, chrom_i: Option<&Chrom>) -> BufReader<File> {
    // using  File.read() did not speed up
    //let fin_chrom = plink::fname_chrom(fin, Some(chrom_i));
    // 512 KB
    //let reader_cap = 512usize * 1024;
    // 1 MB
    let reader_cap = 1usize * 1024 * 1024;
    let reader = BufReader::with_capacity(
        reader_cap,
        File::open(fin_genot.genotype_file(chrom_i)).unwrap(),
    );
    reader
}

fn buf_size_limit(mem: Option<usize>) -> usize {
    let buf_size_limit = match mem {
        Some(x) => x.min(BUF_SIZE_BED_MAX),
        None => {
            log::debug!(
                "Could not get available memory; assume there is {:.3} GB available memory.",
                alloc::mem_gb(BUF_SIZE_BED_MAX)
            );
            BUF_SIZE_BED_MAX
        }
    };
    log::debug!("buf_size_limit: {} GB", alloc::mem_gb(buf_size_limit));
    buf_size_limit
}

fn buf_size_max(mem: Option<usize>, fin_genot: &GenotFile) -> usize {
    let buf_size_limit = buf_size_limit(mem);
    //let buf_size_limit = match mem {
    //    Some(x) => x.min(BUF_SIZE_BED_MAX),
    //    None => {
    //        log::debug!(
    //            "Could not get available memory; assume there is {:.3} GB available memory.",
    //            alloc::mem_gb(BUF_SIZE_BED_MAX)
    //        );
    //        BUF_SIZE_BED_MAX
    //    }
    //};

    log::debug!("buf_size_limit: {:.3} GB", alloc::mem_gb(buf_size_limit));

    let buf_size_max = load_buf_size_max(fin_genot).min(buf_size_limit);

    buf_size_max
}

pub fn group_to_m_in_chrom(
    group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
    m_in_begin: usize,
    m_in_end: usize,
) -> Option<HashMap<usize, Vec<usize>>> {
    if group_to_m_in.is_none() {
        return None;
    }

    let mut group_to_m_in_chrom = HashMap::new();
    let group_to_m_in = group_to_m_in.unwrap();

    for (group_i, m_in_set) in group_to_m_in.iter() {
        let m_in_set_chrom: Vec<usize> = m_in_set
            .iter()
            .filter(|&mi| m_in_begin <= *mi && *mi < m_in_end)
            .map(|&mi| mi - m_in_begin)
            .collect();
        if !m_in_set_chrom.is_empty() {
            group_to_m_in_chrom.insert(*group_i, m_in_set_chrom);
        }
    }

    Some(group_to_m_in_chrom)
}

/// Sequentially load part of bed and convert into predictions since bed could be too large to load on mem.
/// predictions: m*2*n
// load whole is fastest
pub fn generate_genot_plink(
    fin_genot: &GenotFile,
    m: usize, // = m_snv
    n: usize,
    use_snvs: Option<&[bool]>,
    use_samples: Option<&[bool]>,
    //use_missing: bool,
    fill_missing_mode: bool,
    mem: Option<usize>,
    group_to_m_in: Option<HashMap<usize, Vec<usize>>>,
    //m_set: Option<usize>,
) -> Genot {
    log::debug!("to prepare Genot m, n: {}, {}", m, n);

    let mem_buf = mem_buf(m, n, mem);

    // for use_snvs=None
    let use_snvs_v = vec![true; m];
    let use_snvs = match use_snvs {
        Some(x) => x,
        None => &use_snvs_v,
    };
    //let m_in = use_snvs.len();

    let mut g = Genot::new_zeros(m, n);
    //log::debug!("done preparing Genot");

    let (mut g_snvs, mut g_group) = if group_to_m_in.is_some() {
        let (g_snvs, g_group) = g.split_genot(m);
        (g_snvs, Some(g_group))
    } else {
        (g.as_genot_mut(), None)
    };

    // Two patterns to assign predictions
    // bed file is split into chrom or not
    let is_split_chrom = fin_genot.judge_split_chrom();

    if is_split_chrom {
        let mut m_begin = 0;
        let mut m_in_begin = 0;

        // FIXME: fix mem for pg a*
        //log::info!("Temporary fix mem to 64GB.");
        //let mem = Some(64usize * 1024 * 1024 * 1024);
        //todo!("");
        //let mem = alloc::get_available_memory();
        //log::debug!("available mem: {:?} bytes", mem);

        let buf_size_max = buf_size_max(mem, fin_genot);

        log::debug!("buf_size_max {}", buf_size_max);
        let mut buf: Vec<B8_2> = vec![0; buf_size_max];
        for chrom_i in Chrom::variants().iter() {
            log::debug!("Loading chromosome {}", chrom_i);
            let m_in_chrom = genot_io::compute_num_snv_file_chrom(fin_genot, Some(chrom_i));
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

            let group_to_m_in_chrom =
                group_to_m_in_chrom(group_to_m_in.as_ref(), m_in_begin, m_in_end);

            // TODO
            //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

            let reader = buf_read(fin_genot, Some(chrom_i));

            // reuse buf
            buf = assign_genot_toggle(
                &mut g_snvs.as_genot_snvs_mut(m_begin, m_end),
                reader,
                &use_snvs[m_in_begin..m_in_end],
                use_samples,
                Some(buf),
                Some(mem_buf),
                g_group.as_mut(),
                group_to_m_in_chrom.as_ref(),
            );
            m_begin = m_end;
            m_in_begin = m_in_end;
        }
        assert_eq!(m_in_begin, use_snvs.len(), "Sth wrong.");
    } else {
        let reader = buf_read(fin_genot, None);

        assign_genot_toggle(
            &mut g_snvs.as_genot_mut(),
            reader,
            use_snvs,
            use_samples,
            None,
            Some(mem_buf),
            g_group.as_mut(),
            group_to_m_in.as_ref(),
        );
    }

    // missing
    //if !use_missing {
    if fill_missing_mode {
        // print missing rate for group
        if let Some(g_group) = &mut g_group {
            let group_missing_count: Vec<usize> = g_group
                .iter_snv()
                .par_bridge()
                .map(|g_snv| g_snv.count_missing())
                .collect();

            let (hist_missing, hist_bins) = vec::histogram(&group_missing_count, 10);
            log::info!("Histgram of missing count for group");
            hist_bins
                .iter()
                .zip(hist_missing.iter())
                .for_each(|(bin, count)| {
                    log::info!("<{}: {}", bin, count);
                });
            log::info!(
                ">={}: {}",
                hist_bins.last().unwrap(),
                hist_missing.last().unwrap()
            );
        }

        // both snvs and group
        g.iter_snv_mut()
            .par_bridge()
            .for_each(|mut g_snv| g_snv.fill_missing_mode());

        //g_snvs
        //    .iter_snv_mut()
        //    .par_bridge()
        //    .for_each(|mut g_snv| g_snv.fill_missing_mode());
    }

    g
}


fn assign_genot_toggle<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
    buf: Option<Vec<B8_2>>,
    mem: Option<usize>,
    // mut is necessary
    // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
    mut g_group: Option<&mut GenotMut>,
    group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
) -> Vec<B8_2> {


    let m = vec::count_true(&use_snvs);

    //if (m as f64)/(use_snvs.len() as f64)<0.1{

    // error on ukbrap
    //    log::debug!("use load part()");
    //    assign_genot_loadpart(
    //        g_chrom,
    //        reader,
    //        use_snvs,
    //        use_samples,
    //        buf,
    //        mem,
    //        g_group,
    //        group_to_m_in,
    //        )
    //}else{
    //    log::debug!("use load whole()");
    assign_genot(
            g_chrom,
            reader,
            use_snvs,
            use_samples,
            buf,
            mem,
            g_group,
            group_to_m_in,
            )
    //}
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
    mem: Option<usize>,
    // mut is necessary
    // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
    mut g_group: Option<&mut GenotMut>,
    group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
) -> Vec<B8_2> {
    let m_in_chrom = use_snvs.len();

    let n_in = match use_samples {
        Some(v) => v.len(),
        None => g_chrom.n(),
    };

    // FIXME: available mem?
    // check if 32 GB? 16GB? remains or not
    // reading smaller than buf_size is not error but reason is unknown
    // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
    // -> to deal with this, every time use Seek to change position
    let buf_size_limit = buf_size_limit(mem);

    let byte_per_snv = bed_per_snv_size(n_in);
    let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
    let buf_num_snv: usize = buf_num_snv_limit.min(m_in_chrom);
    let buf_size: usize = buf_num_snv * byte_per_snv;
    assert_eq!(buf_size % byte_per_snv, 0);
    assert!(buf_size <= buf_size_limit);

    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    // TODO: check buf_size < available

    let mut buf: Vec<B8_2> = match buf {
        Some(v) => v,
        None => vec![0; buf_size],
    };

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

        let (m_to_m_in, m_read) = genotype::create_m_to_m_in(use_snvs_loaded);

        let m_end_loaded = m_begin_loaded + m_read;

        if m_read != 0 {
            // use fill() might be fast?
            // This might be smaller than buf
            let buf_size_ = m_in_read * byte_per_snv;
            buf.resize(buf_size_, 0);
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

                    assign_gsnv_from_bed(&mut g_snv, &buf_mi, use_samples);
                });

            // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
            if let Some(g_group) = &mut g_group {
                let group_to_m_in = group_to_m_in.unwrap();

                let n = match use_samples {
                    Some(v) => vec::count_true(v),
                    None => n_in,
                };

                g_group
                    .iter_snv_mut()
                    .enumerate()
                    .par_bridge()
                    .for_each_with(
                        GenotSnv::new_empty(n),
                        |g_snv_tmp, (mi_loaded, mut g_snv)| {
                            let group_ids_mi = group_to_m_in.get(&mi_loaded);

                            if let Some(group_ids_mi) = group_ids_mi {
                                group_ids_mi.iter().for_each(|&group_id| {
                                    let m_in_read_i = m_to_m_in[&group_id];
                                    let buf_mi = &buf[m_in_read_i * byte_per_snv
                                        ..(m_in_read_i + 1) * byte_per_snv];

                                    // initialize g_snv_tmp
                                    g_snv_tmp.fill_0();

                                    assign_gsnv_from_bed(
                                        &mut g_snv_tmp.as_genot_snv_mut_snv(),
                                        &buf_mi,
                                        use_samples,
                                    );

                                    g_snv.or_binary(&g_snv_tmp.as_genot_snv());
                                });
                            }
                            // otherwise, do nothing
                        },
                    );
            }

            // Error
            //if let Some(x) = g_group {
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

// // TODO: use SIMD
//// slower using .read() than .read_exact()
//// use when file size is too large
//fn assign_predictions_chunk_old2<R: BufRead + Seek>(
//    g_chrom: &mut GenotMut,
//    mut reader: R,
//    use_snvs: &[bool],
//    use_samples: Option<&[bool]>,
//    //buf: Option<Vec<B8_2>>,
//    mem: Option<usize>,
//    // mut is necessary
//    // https://stackoverflow.com/questions/74763962/when-is-mut-required-when-passing-an-optionmut-t
//    mut g_group: Option<&mut GenotMut>,
//    group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
//) {
//    let m_in_chrom = use_snvs.len();
//
//    let n_in = match use_samples {
//        Some(v) => v.len(),
//        None => g_chrom.n(),
//    };
//
//    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();
//
//    // check if 32 GB? 16GB? remains or not
//    // reading smaller than buf_size is not error but reason is unknown
//    // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
//    // -> to deal with this, every time use Seek to change position
//    //let buf_size_limit: usize = 64 * 1024 * 1024 * 1024;
//    let buf_size_limit: usize = 4 * 1024 * 1024 * 1024;
//
//    let byte_per_snv = bed_per_snv_size(n_in);
//    let buf_num_snv_limit: usize = buf_size_limit / byte_per_snv;
//    let buf_num_snv: usize = buf_num_snv_limit.min(m_in_chrom);
//    let buf_size: usize = buf_num_snv * byte_per_snv;
//    assert_eq!(buf_size % byte_per_snv, 0);
//    assert!(buf_size <= buf_size_limit);
//
//    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
//    // TODO: check buf_size < available
//    let mut buf: Vec<B8_2> = vec![0; buf_size];
//
//    let mut m_in_begin_loaded = 0;
//    let mut m_begin_loaded = 0;
//    // read buf length in one loop
//    loop {
//        //log::debug!(
//        //    "m_in_begin, m_begin {},{}",
//        //    m_in_begin_loaded, m_begin_loaded
//        //);
//
//        // next start is m_in_begin_load
//        let buf_next_begin = 3 + m_in_begin_loaded * byte_per_snv;
//        //log::debug!("buf_next_begin {}", buf_next_begin);
//        reader.seek(SeekFrom::Start(buf_next_begin as u64)).unwrap();
//        // https://stackoverflow.com/questions/37079342/what-is-the-most-efficient-way-to-read-a-large-file-in-chunks-without-loading-th
//        // use fill() might be fast?
//        // This might be smaller than buf
//        // read_exact??
//        //let loaded_byte = reader.read_exact(&mut buf).unwrap();
//        let loaded_byte = reader.read(&mut buf).unwrap();
//        //if loaded_byte != buf_size, last buffer or could not load whole buffer
//        log::debug!("loaded: {} {}", loaded_byte, buf_size);
//        //}
//        if loaded_byte == 0 {
//            // all buf read
//            break;
//        } else {
//            //log::debug!("buf {:?}", buf.len());
//            // might fail when only part was loaded
//            //assert_eq!(loaded_byte % byte_per_snv, 0);
//            let m_in_read = loaded_byte / byte_per_snv;
//            // if m_in_read==0, loaded_byte < byte_per_snv, which means loaded_byte is less than one SNV.
//            assert!(m_in_read > 0, "loaded_byte is less than byte of one SNV.");
//
//            // x plan 1 create vec of pred_use and  pred_use.iter()...
//            //  -> unsafe when pred_use is duplicated
//            // x plan 2 predictions.par_chunks_mut().filter(use)
//            // o plan 3 predictions_chrom[m_buf..m_buf_end].par_chunks_mut()
//
//            let m_in_end_loaded = m_in_begin_loaded + m_in_read;
//            let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];
//
//            let (m_to_m_in, m_read) = genot_index::create_m_to_m_in(use_snvs_loaded);
//
//            // here for when m_read==0
//            // or move below and do not continue for m_read==0 since unnecessary
//            //m_in_begin_loaded = m_in_end_loaded;
//
//            //if m_read == 0 {
//            //    continue;
//            //}
//
//            let m_end_loaded = m_begin_loaded + m_read;
//
//            //log::debug!("m_in_end, m_end {},{}", m_in_end_loaded, m_end_loaded);
//
//            let mut g_chrom_part = g_chrom.as_genot_snvs_mut(m_begin_loaded, m_end_loaded);
//
//            g_chrom_part
//                .iter_snv_mut()
//                .enumerate()
//                .par_bridge()
//                .for_each(|(mi_loaded, mut g_snv)| {
//                    let m_in_read_i = m_to_m_in[&mi_loaded];
//                    let buf_mi = &buf[m_in_read_i * byte_per_snv..(m_in_read_i + 1) * byte_per_snv];
//
//                    assign_pred_from_bed(&mut g_snv, &buf_mi, use_samples);
//                });
//
//            m_begin_loaded = m_end_loaded;
//            // bug! when m_read==0, will not renewed
//            m_in_begin_loaded = m_in_end_loaded;
//        }
//    }
//    assert_eq!(m_in_chrom, m_in_begin_loaded);
//}




// error on rap
// rayon impossible since cannot read inside
fn assign_genot_loadpart<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
    // no use inside
    buf_in: Option<Vec<B8_2>>,
    _mem: Option<usize>,
    mut g_group: Option<&mut GenotMut>,
    group_to_m_in: Option<&HashMap<usize, Vec<usize>>>,
)->Vec<B8_2> {

    if g_group.is_some(){
        unimplemented!("ny");
    }
    if group_to_m_in.is_some(){
        unimplemented!("ny");
    }


    //let m_in_chrom = use_snvs.len();

    // assume m and m_in in m_to_m_in are sorted
    // FIXME: use_snvs should be only chrom
    //let m_to_m_in = create_m_chrom_to_m_in_chrom(use_snvs);
    //let m_in_to_m = create_m_in_to_m_chrom(use_snvs, m_begin);

    let n_in = match use_samples {
        Some(v) => v.len(),
        None => g_chrom.n(),
    };

    // for debug
    //log::debug!("ys[0]:{:#010b}", ys[0]);

    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    let (m_to_m_in, m_read) = genotype::create_m_to_m_in(use_snvs);
    //let (m_to_m_in, m_for_read) = genot_index::create_m_to_m_in(use_snvs);

    let byte_per_snv = bed_per_snv_size(n_in);
    //let byte_per_snv = plink::bed_per_snv_size(n_in);
    let buf_size = byte_per_snv;

    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    let mut buf: Vec<B8_2> = vec![0; buf_size];
    assert_eq!(buf.len(), buf_size);

    assert_eq!(g_chrom.m(), m_read, "Sth wrong");

    // rayon impossible since cannot seek
    g_chrom
        .iter_snv_mut()
        .enumerate()
        .for_each(|(m_read_i, mut g_snv)| {

            log::debug!("m_read_i: {}", m_read_i);

            // for debug
            log::debug!("m_read_in: {}", m_to_m_in[&m_read_i]);
            let buf_begin = 3 + m_to_m_in[&m_read_i] * byte_per_snv;
            reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
            reader.read_exact(&mut buf).unwrap();
            assert_eq!(buf.len(), buf_size);
            //let loaded_byte = reader.read(&mut buf).unwrap();
            //assert_eq!(loaded_byte, byte_per_snv);

            //let g_snv = g_chrom.to_genot_twin_snv_mut(m_read_i);
            //let pred = predictions_snv_s_mut(predictions_chrom, m_read_i, n);
            assign_gsnv_from_bed(&mut g_snv, &buf, use_samples);
        });


    // return buf
    // TODO: all right when calling assign_genot() after loadpart
    match buf_in {
        Some(v) => v,
        None => buf,
    }
    
}



//// rayon impossible since cannot read inside
//pub fn assign_genot_file_loadpart(
//    g_chrom: &mut GenotMut,
//    fin_chrom: &str,
//    use_snvs: &[bool],
//    use_samples: &[bool],
//) {
//    let m_in_chrom = use_snvs.len();
//    // log::debug!("use_samples: {:?}", use_samples);
//
//    // assume m and m_in in m_to_m_in are sorted
//    // FIXME: use_snvs should be only chrom
//    //let m_to_m_in = create_m_chrom_to_m_in_chrom(use_snvs);
//    //let m_in_to_m = create_m_in_to_m_chrom(use_snvs, m_begin);
//
//    let n_in = use_samples.len();
//
//    // for debug
//    //log::debug!("ys[0]:{:#010b}", ys[0]);
//
//    plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();
//
//    let (m_to_m_in, m_for_read) = genot_index::create_m_to_m_in(use_snvs);
//
//    let byte_per_snv = plink::bed_per_snv_size(n_in);
//    let buf_size = byte_per_snv;
//
//    let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
//    let mut buf: Vec<B8_2> = vec![0; buf_size];
//
//    assert_eq!(g_chrom.m(), m_for_read, "Sth wrong");
//
//    // rayon impossible since cannot seek
//    g_chrom
//        .iter_snv_mut()
//        .enumerate()
//        .for_each(|(m_read_i, mut g_snv)| {
//            let buf_begin = 3 + m_to_m_in[&m_read_i] * byte_per_snv;
//            reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
//            let loaded_byte = reader.read(&mut buf).unwrap();
//            assert_eq!(loaded_byte, byte_per_snv);
//
//            //let g_snv = g_chrom.to_genot_twin_snv_mut(m_read_i);
//            //let pred = predictions_snv_s_mut(predictions_chrom, m_read_i, n);
//            assign_pred_from_bed(&mut g_snv, &buf, use_samples);
//        });
//}
//
pub fn bed_per_snv_size(n: usize) -> usize {
    (n + 3) / 4
}

pub fn calculate_bed_size_genotype(m: usize, n: usize) -> usize {
    m * bed_per_snv_size(n)
    //m * ((n + 3) / 4)
}

pub fn calculate_bed_size(m: usize, n: usize) -> usize {
    3 + calculate_bed_size_genotype(m, n)
    //3 + m * ((n + 3) / 4)
}

// move to mod plink
//  TODO: should unwrap here?
/// return bed_size if valid, error otherwise
//pub fn check_valid_bed(fin: &str, n: usize, m: usize) -> Result<usize, String> {
pub fn check_valid_bed(
    //fin: &Path,
    //gfmt: GenotFormat,
    fin_genot: &GenotFile,
    chrom: Option<&Chrom>,
    m: usize,
    n: usize,
) -> Result<usize, Box<dyn Error>> {
    let fin_bed = fin_genot.genotype_file(chrom);
    //let fin_bed = fname_plinks_genot(fin, gfmt, chrom);
    // check if open
    let mut reader = File::open(fin_bed)?;

    // check if size is correct
    let f_end: usize = reader.seek(SeekFrom::End(0)).unwrap() as usize;
    log::debug!("file end {}", f_end);
    let bed_size = calculate_bed_size(m, n);
    if f_end != bed_size {
        return Err(format!(
            "File size of .bed is wrong: {} vs correct {}.",
            f_end, bed_size
        )
        .into());
    }

    // check if the first 3 bytes are correct.
    reader.seek(SeekFrom::Start(0)).unwrap();
    let mut buf: Vec<u8> = Vec::with_capacity(n);
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    //log::debug!("{:?}", buf);
    if buf != vec![0x6cu8, 0x1b, 0x01] {
        return Err("Magic number of .bed file is wrong.".into());
    }
    Ok(bed_size)
}

/// g_snv must be all false
/// use g_snv.fill_0() to make sure
pub fn assign_gsnv_from_bed(
    g_snv: &mut GenotSnvMut,
    buf_mi: &[B8_2],
    use_samples: Option<&[bool]>,
) {
    if let Some(use_sample) = use_samples {
        let mut ni = 0;
        for (n_in_i, v) in use_sample.iter().enumerate() {
            if *v {
                // bedcode
                let bcode = buf_to_ped_code(buf_mi, n_in_i);
                g_snv.set_bed_code_init_unchecked(bcode, ni);
                ni += 1;
            }
        }
    } else {
        for ni in 0..g_snv.n() {
            // bedcode
            let bcode = buf_to_ped_code(buf_mi, ni);
            g_snv.set_bed_code_init_unchecked(bcode, ni);
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

// as fast as load_byte2 but complecated
#[allow(dead_code)]
fn load_byte1(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    //let bed_size = calculate_bed_size(n, m);
    let v_size = calculate_bed_size_genotype(m, n);
    log::debug!("{}", v_size);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }
    let mut reader = BufReader::new(File::open(fin).unwrap());
    let n: usize = 100_000_000;
    //let n: usize = 1_000_000;
    let mut buf: Vec<u8> = Vec::with_capacity(n);

    // first for 3
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    log::debug!("{:?}", buf);

    unsafe {
        buf.set_len(n);
    }
    let mut loop_times = 0;
    // if we know the size is 9, then
    for i in 0..((v_size - 1) / n) {
        reader.read_exact(&mut buf).unwrap();
        // https://stackoverflow.com/questions/28219231/how-to-idiomatically-copy-a-slice
        v[n * i..n * (i + 1)].copy_from_slice(&buf);
        loop_times += 1;
    }
    log::debug!("loop times{}", loop_times);
    // for remaining
    let n_remain: usize = v_size - (v_size - 1) / n * n;
    unsafe {
        buf.set_len(n_remain);
    }
    reader.read_exact(&mut buf).unwrap();
    v[(v_size - 1) / n * n..v_size].copy_from_slice(&buf);
    v
}

// fast and simple
/// Use this function to load bytes.
#[allow(dead_code)]
fn load_byte2(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    // you can use "seek"
    // use this!!
    //let bed_size = calculate_bed_size(n, m);
    let v_size = calculate_bed_size_genotype(m, n);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }
    log::debug!("v_size: {}", v_size);

    let mut reader = BufReader::new(File::open(fin).unwrap());
    let n: usize = 100_000_000;
    //let n: usize = 1_000_000;
    let mut buf: Vec<u8> = Vec::with_capacity(n);

    // first for 3
    unsafe {
        buf.set_len(3);
    }
    reader.read_exact(&mut buf).unwrap();
    log::debug!("{:?}", buf);

    unsafe {
        buf.set_len(n);
    }

    let mut loop_times = 0;
    let mut i_ptr: usize = 0;
    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            n => {
                log::debug!("i_ptr,n:{},{}", i_ptr, n);
                let buf = &buf[..n];
                v[i_ptr..i_ptr + n].copy_from_slice(&buf);
                i_ptr = i_ptr + n;
            }
        }
        loop_times += 1;
    }
    log::debug!("loop times{}", loop_times);
    v
}

#[allow(dead_code)]
fn load_byte3(fin: &str, n: usize, m: usize) -> Vec<u8> {
    let mut file = File::open(fin).unwrap();
    // this is ok but not sure fast or slow
    let mut buf: Vec<u8> = Vec::new();
    // meaningless, lastly cap is 32
    //let n = 9;
    //let mut buf: Vec<u8> = Vec::with_capacity(n + 2);
    // this wont work
    //unsafe {
    //    buf.set_len(n);
    //}
    log::debug!("len: {}", buf.len());
    log::debug!("cap: {}", buf.capacity());
    let _ = file.read_to_end(&mut buf).unwrap();
    log::debug!("len: {}", buf.len());
    // larger than cap
    log::debug!("cap: {}", buf.capacity());
    log::debug!("buf[0]: {}", buf[0]);
    //log::debug!("{:?}", buf);

    //let v_size = calculate_bed_genotype_size(n, m);
    //tmp
    let v_size = calculate_bed_size(n, m);
    let mut v: Vec<u8> = Vec::with_capacity(v_size);
    unsafe {
        v.set_len(v_size);
    }

    (&mut v).copy_from_slice(&buf);
    v
}

pub fn run_byte(fin: &Path, n: usize, m: usize) -> Vec<u8> {
    // fast but complecated
    log::debug!("way 1");
    let istart = Instant::now();
    let v = load_byte1(fin, n, m);
    log::debug!("load_byte1: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);

    // as fast as way 1
    // also easy to write
    log::debug!("way 2");
    let istart = Instant::now();
    let v = load_byte2(fin, n, m);
    log::debug!("load_byte2: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);

    /*
    // seems slow -> could be because file size was larger than mem size?
    log::debug!("way 3");
    let istart = Instant::now();
    let v = load_byte3(fin, n, m);
    log::debug!("load_byte3: {:?}", Instant::now().duration_since(istart));
    log::debug!("{}", v[0]);
    */

    v
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::samples;
    use crate::{sample, snv};
    use std::{io::Cursor, path::PathBuf};

    fn setup_test() -> (GenotFile, Vec<bool>, usize, usize, Vec<bool>, Vec<bool>) {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        //let gfmt = GenotFormat::Plink1;
        let fin_genot = GenotFile::Plink1(fin);
        //let fin_snv = None;
        //let fin_sample = None;

        let m_in: usize = genot_io::compute_num_snv(&fin_genot).unwrap();
        log::debug!("{}", m_in);
        //let n_in: usize = genot_io::compute_num_sample(&fin_genot).unwrap();
        //println!("{}", n_in);
        // load snvs
        let snvs_in = genot_io::load_snvs(&fin_genot);
        let (use_snvs, m) = snv::make_use_snvs_buf(None, &snvs_in);
        //let (m, use_snvs) = snv::make_use_snvs_buf(None, &snvs_in);
        //let (m, use_snvs) = snv::make_use_snvs(fin_snv, &snvs_in);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (use_samples, n) = sample::make_use_samples_buf(None, &fin_genot);
        //let (n, use_samples) = sample::make_use_samples(fin_sample, &fin_genot);
        //let ys: Vec<bool> = io_genot::load_ys(&fin, gfmt, None, None, &use_samples);
        let sample_id_to_n = samples::create_sample_id_to_n(&fin_genot, Some(&use_samples));
        let ys: Vec<bool> = genot_io::load_ys_buf(&fin_genot, None, None, &sample_id_to_n).unwrap();

        (fin_genot, ys, m, n, use_snvs, use_samples)
    }

    #[test]
    fn test_group_to_m_in_chrom() {
        let group_to_m_in: HashMap<usize, Vec<usize>> = HashMap::from_iter(vec![
            (0, vec![0, 1]),
            (1, vec![2, 3]),
            (2, vec![4, 5]),
            (3, vec![6, 7]),
            (4, vec![8, 9]),
            (5, vec![10, 5]),
        ]);

        let m_in_begin = 3;
        let m_in_end = 7;

        let group_to_m_in_chrom = group_to_m_in_chrom(Some(&group_to_m_in), m_in_begin, m_in_end);

        // extracted in the range
        // (1, [3])
        // (2, [4, 5])
        // (3, [6])
        // (5, [5])
        //
        // after adjusting for m_in_begin
        // (1, [0])
        // (2, [1, 2])
        // (3, [3])
        // (5, [2])
        //
        let group_to_m_in_chrom_ans: HashMap<usize, Vec<usize>> = HashMap::from_iter(vec![
            (1, vec![0]),
            (2, vec![1, 2]),
            (3, vec![3]),
            (5, vec![2]),
        ]);

        assert_eq!(group_to_m_in_chrom.unwrap(), group_to_m_in_chrom_ans);
    }

    #[test]
    fn test_generate_genot_plink_whole() {
        let (fin_genot, _, m, n, use_snvs, use_samples) = setup_test();
        //let gfmt = GenotFormat::Plink1;
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        let g = generate_genot_plink(
            &fin_genot,
            //gfmt,
            m,
            n,
            Some(&use_snvs),
            Some(&use_samples),
            false,
            None,
            None,
        );
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
    fn test_assign_gsnv_from_bed() {
        let mut g = GenotSnv::new_empty(4);
        // [2, 0, 3, 0, 1, 0]
        let pbuf = vec![0b11_01_11_00, 0b00_00_11_10];

        let use_samples = vec![true, false, true, true, true, false];
        assign_gsnv_from_bed(&mut g.as_genot_snv_mut_snv(), &pbuf, Some(&use_samples));
        assert_eq!(g.vals(), vec![2, 3, 0, 1]);
    }

    #[test]
    fn test_assign_genot() {
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
            None,
            None,
            None,
        );

        assert_eq!(g.vals_snv(0), vec![2, 0, 3]);
        assert_eq!(g.vals_snv(1), vec![0, 2, 1]);
    }

    #[test]
    fn test_assign_genot_group() {
        let mut g_snvs = Genot::new_zeros(3, 5);
        let mut g_group = Genot::new_zeros(2, 5);

        let group_m_to_m_in: HashMap<usize, Vec<usize>> =
            HashMap::from_iter(vec![(0, vec![0, 1]), (1, vec![0, 1, 2])]);
        //HashMap::from_iter(vec![(1, vec![0, 1, 2])]);

        // snv0: [0, 0, 0, 1, 3]
        // snv1: [0, 0, 1, 0, 0]
        // snv2: [0, 1, 2, 3, 1]
        // set0: snv0,1
        // set1:  snv0,1,2
        // set0: [0, 0, 1, 1, 3]
        // set1:  [0, 1, 1, 3, 3]
        let reader: Vec<u8> = vec![
            0x6c,
            0x1b,
            0x01,
            0b10_11_11_11,
            0b00_00_00_01,
            0b11_10_11_11,
            0b00_00_00_11,
            0b01_00_10_11,
            0b00_00_00_10,
        ];
        let mut cur = Cursor::new(reader.as_slice());
        assign_genot(
            &mut g_snvs.as_genot_mut(),
            &mut cur,
            &[true, true, true],
            None,
            None,
            None,
            Some(&mut g_group.as_genot_mut()),
            Some(&group_m_to_m_in),
        );

        assert_eq!(g_snvs.vals_snv(0), vec![0, 0, 0, 1, 3]);
        assert_eq!(g_snvs.vals_snv(1), vec![0, 0, 1, 0, 0]);
        assert_eq!(g_snvs.vals_snv(2), vec![0, 1, 2, 3, 1]);
        assert_eq!(g_group.vals_snv(0), vec![0, 0, 1, 1, 3]);
        assert_eq!(g_group.vals_snv(1), vec![0, 1, 1, 3, 3]);
    }

    /*
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
    */

    #[test]
    fn test_get_bed_size() {
        let m: usize = 3;
        let n: usize = 5;
        let bed_size: usize = calculate_bed_size(m, n);
        assert_eq!(bed_size, 9);

        let m: usize = 3;
        let n: usize = 4;
        let bed_size: usize = calculate_bed_size(m, n);
        assert_eq!(bed_size, 6);
    }

    #[test]
    fn test_check_valid_bed() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        //let gfmt = GenotFormat::Plink1;
        let fin_genot = GenotFile::Plink1(fin);

        let m: usize = genot_io::compute_num_snv(&fin_genot).unwrap();
        let n: usize = genot_io::compute_num_sample(&fin_genot).unwrap();
        log::debug!("m,n: {},{}", m, n);
        let bed_size = check_valid_bed(&fin_genot, None, m, n).unwrap();
        assert_eq!(bed_size, calculate_bed_size(m, n));
    }

    #[test]
    #[should_panic]
    fn test_check_valid_bed_panic_nofile() {
        // file does not exist
        let fin = PathBuf::from("./toy");
        let fin_genot = GenotFile::Plink1(fin);

        check_valid_bed(&fin_genot, None, 3, 2).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_check_valid_bed_panic() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let fin_genot = GenotFile::Plink1(fin);

        let m: usize = genot_io::compute_num_snv(&fin_genot).unwrap();
        let n: usize = genot_io::compute_num_sample(&fin_genot).unwrap();
        // size is wrong
        check_valid_bed(&fin_genot, None, m, n - 1).unwrap();
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

    #[test]
    fn test_load_byte2() {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        //let gfmt = GenotFormat::Plink1;
        let fin_genot = GenotFile::Plink1(fin);

        let fin_fam = fin_genot.sample_file(None);
        //let fin_fam = genot_io::fname_plinks_sample(&fin, gfmt, None);
        //let fin_fam = fin.clone() + ".fam";
        let n: usize = textfile::compute_num_line_text(&fin_fam, None).unwrap();

        let fin_bim = fin_genot.snv_file(None);
        //let fin_bim = genot_io::fname_plinks_snv(&fin, gfmt, None);
        //let fin_bim = fin.clone() + ".bim";
        let m: usize = textfile::compute_num_line_text(&fin_bim, None).unwrap();

        let fin_bed = fin_genot.genotype_file(None);
        //let fin_bed = genot_io::fname_plinks_genot(&fin, gfmt, None);
        //let fin_bed = fin.clone() + ".bed";
        load_byte2(&fin_bed, n, m);
    }
}
