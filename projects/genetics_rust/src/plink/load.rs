use crate::genot::prelude::*;
use crate::genot_index;
use crate::{plink, vec, Chrom};
use rayon::prelude::*;
use std::fs::File;
use std::io::BufRead;
use std::io::{prelude::*, BufReader, SeekFrom};
use std::path::Path;

// TODO: mi for multiple files
// TODO: remove n
/// mi is in fin file
pub fn load_snv(
    fin: &Path,
    mi: usize,
    n: usize,
    use_samples: Option<&[bool]>,
    use_missing: bool,
) -> GenotSnv {
    let reader = BufReader::new(File::open(plink::fname_bed(fin, None)).unwrap());

    let mut g_snv = GenotSnv::new_empty(n);

    load_snv_read(
        &mut g_snv.as_genot_snv_mut_snv(),
        reader,
        mi,
        use_samples,
        n,
    );

    if !use_missing {
        fill_missing(&mut g_snv.as_genot_snv_mut_snv());
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
    let byte_per_snv = plink::bed_per_snv_size(n_in);
    let buf_size = byte_per_snv;
    let mut buf: Vec<B8_2> = vec![0; buf_size];
    let buf_begin = 3 + mi * byte_per_snv;
    reader.seek(SeekFrom::Start(buf_begin as u64)).unwrap();
    let loaded_byte = reader.read(&mut buf).unwrap();
    assert_eq!(loaded_byte, byte_per_snv);

    assign_pred_from_bed(g_snv, &buf, use_samples);
}

/// Generate Vector of predictions.
/// Sequentially load part of bed and convert into predictions since bed could be too large to load on mem.
/// predictions: m*2*n
// load whole is fastest
pub fn generate_genot(
    fin: &Path,
    //ys: &[B8],
    m: usize,
    n: usize,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
    use_missing: bool,
) -> Genot {
    println!("to prepare Genot m, n: {}, {}", m, n);
    let mut g = Genot::new_empty(m, n);
    //println!("done preparing Genot");

    // FIXME: hang due to this part...
    // program stops and do not go next
    //println!(
    //    "len vs total mem, {} bytes vs {} bytes",
    //    g.len(),
    //    alloc::get_total_memory()
    //);
    //if g.len() > alloc::get_total_memory() {
    //    panic!("Memory insufficient.")
    //}

    // Two patterns to assign predictions
    // bed file is split into chrom or not

    let is_split_chrom = plink::judge_split_chrom(fin);

    if is_split_chrom {
        let mut m_begin = 0;
        let mut m_in_begin = 0;
        for chrom_i in Chrom::variants().iter() {
            println!("Loading chromosome {}", chrom_i);
            let m_in_chrom = plink::compute_num_snv_chrom(fin, Some(chrom_i));
            if let None = m_in_chrom {
                continue;
            }

            let m_in_chrom = m_in_chrom.unwrap();
            let m_in_end = m_in_begin + m_in_chrom;
            // count non-zero
            let m_chrom = vec::count_true(&use_snvs[m_in_begin..m_in_end]);
            let m_end = m_begin + m_chrom;

            // TODO
            //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();
            // TODO: no use fname_chrom()
            let fin_chrom = plink::fname_chrom(fin, Some(chrom_i));
            let reader = BufReader::new(File::open(plink::fname_bed(&fin_chrom, None)).unwrap());

            assign_predictions(
                &mut g.as_genot_snvs_mut(m_begin, m_end),
                reader,
                &use_snvs[m_in_begin..m_in_end],
                use_samples,
            );
            m_begin = m_end;
            m_in_begin = m_in_end;
        }
        assert_eq!(m_in_begin, use_snvs.len(), "Sth wrong.");
    } else {
        let reader = BufReader::new(File::open(plink::fname_bed(&fin, None)).unwrap());
        assign_predictions(&mut g.as_genot_mut(), reader, use_snvs, use_samples);
    }

    // missing
    if !use_missing {
        g.iter_snv_mut()
            .par_bridge()
            .for_each(|mut g_snv| fill_missing(&mut g_snv));
    }

    g
}

/*
// wrongly assumes reader.read(&buf) always load buf size
fn assign_predictions_old<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    //mut reader: BufReader<File>,
    mut reader: R,
    //fin_chrom: &str,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
) {
    //println!("g_chrom len: {}", g_chrom.len());

    let m_in_chrom = use_snvs.len();
    //println!("m_in_chrom: {}", m_in_chrom);
    //println!("use_samples: {:?}", use_samples);

    // assume m and m_in in m_to_m_in are sorted
    // FIXME: use_snvs should be only chrom
    //let m_to_m_in = create_m_chrom_to_m_in_chrom(use_snvs);
    //let m_in_to_m = create_m_in_to_m_chrom(use_snvs, m_begin);

    let n_in = if let Some(v) = use_samples {
        v.len()
    } else {
        g_chrom.n()
    };
    //println!("n_in: {}", n_in);
    //let n_in = use_samples.len();

    //println!("m_in len {}", plink::bed_per_snv_size(n_in) * m_in_chrom);

    // for debug
    //println!("ys[0]:{:#010b}", ys[0]);

    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    // check if 32 GB? 16GB? remains or not

    // TODO: somehow only load part of buf_size
    // If violate `assert_eq!(loaded_byte % byte_per_snv, 0);`, probably because only loading part of buf.
    // -> make buf_size_limit smaller will solve the problem
    // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
    // reading smaller than buf_size is not error but reason is unknown
    // -> should to deal with this?
    // -> YOU SHOULD
    // TODO: to deal with this, every time use Seek to change position
    // only 2,147,479,552 byte was read
    //let buf_size_limit: usize = 4 * 1024 * 1024 * 1024; //error
    //let buf_size_limit: usize = 3 * 1024 * 1024 * 1024; //error
    let buf_size_limit: usize = 2 * 1024 * 1024 * 1024; // ok
                                                        //let buf_size_limit: usize = 1 * 1024 * 1024 * 1024; //ok
                                                        //let buf_size_limit: usize = 16 * 1024 * 1024;
                                                        //let buf_size_limit: usize = 32 * 1024 * 1024;
    let byte_per_snv = plink::bed_per_snv_size(n_in);
    //let byte_per_snv: usize = (n_in + 3) / 4;
    let num_snv_per_buf: usize = buf_size_limit / byte_per_snv;
    //let buf_num_snv: usize = num_snv_per_buf;
    let buf_num_snv: usize = num_snv_per_buf.min(m_in_chrom);
    //println!("buf_num_snv {}", buf_num_snv);
    let buf_size: usize = buf_num_snv * byte_per_snv;
    //let buf_size: usize = num_snv_per_buf * byte_per_snv;
    assert_eq!(buf_size % byte_per_snv, 0);
    assert!(buf_size <= buf_size_limit);
    //println!("buf_size {}", buf_size);

    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    let mut buf: Vec<B8_2> = vec![0; buf_size];

    // skip the magic numbers
    // How to skip?
    // https://stackoverflow.com/questions/42243355/how-to-advance-through-data-from-the-stdioread-trait-when-seek-isnt-impleme

    //let loaded_byte = reader.read(&mut vec![0; 3]).unwrap();
    // Error: somehow consume raise error when loading large file.
    //reader.consume(3);
    // This requires Seek but easy to use
    reader.seek(SeekFrom::Start(3)).unwrap();

    let mut m_in_begin_loaded = 0;
    let mut m_begin_loaded = 0;
    // read buf length in one loop
    loop {
        //println!(
        //    "m_in_begin, m_begin {},{}",
        //    m_in_begin_loaded, m_begin_loaded
        //);

        let loaded_byte = reader.read(&mut buf).unwrap();
        if loaded_byte != buf_size {
            // last buffer or error
            println!("loaded: {} {}", loaded_byte, buf_size);
        }
        if loaded_byte == 0 {
            // all buf read
            break;
        } else {
            //println!("buf {:?}", buf.len());
            assert_eq!(loaded_byte % byte_per_snv, 0);
            let m_in_read = loaded_byte / byte_per_snv;

            // x plan 1 create vec of pred_use and  pred_use.iter()...
            //  -> unsafe when pred_use is duplicated
            // x plan 2 predictions.par_chunks_mut().filter(use)
            // o plan 3 predictions_chrom[m_buf..m_buf_end].par_chunks_mut()

            let m_in_end_loaded = m_in_begin_loaded + m_in_read;
            let use_snvs_loaded = &use_snvs[m_in_begin_loaded..m_in_end_loaded];

            let (m_to_m_in, m_read) = genot_index::create_m_to_m_in(use_snvs_loaded);
            //let (m_to_m_in, m_for_read) = create_m_to_m_in(use_snvs_loaded);

            // here for when m_read==0
            m_in_begin_loaded = m_in_end_loaded;

            if m_read == 0 {
                continue;
            }

            let m_end_loaded = m_begin_loaded + m_read;

            //println!("m_in_end, m_end {},{}", m_in_end_loaded, m_end_loaded);

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
            //m_in_begin_loaded = m_in_end_loaded;
        }
    }
    assert_eq!(m_in_chrom, m_in_begin_loaded);

    //let m_in_end = m_in_begin_loaded;
    //m_in_end
}
*/

// TODO: use SIMD
// solve read() error
fn assign_predictions<R: BufRead + Seek>(
    g_chrom: &mut GenotMut,
    mut reader: R,
    use_snvs: &[bool],
    use_samples: Option<&[bool]>,
) {
    //println!("g_chrom len: {}", g_chrom.len());

    let m_in_chrom = use_snvs.len();
    //println!("m_in_chrom: {}", m_in_chrom);
    //println!("use_samples: {:?}", use_samples);

    // assume m and m_in in m_to_m_in are sorted
    // FIXME: use_snvs should be only chrom
    //let m_to_m_in = create_m_chrom_to_m_in_chrom(use_snvs);
    //let m_in_to_m = create_m_in_to_m_chrom(use_snvs, m_begin);

    let n_in = if let Some(v) = use_samples {
        v.len()
    } else {
        g_chrom.n()
    };
    //println!("n_in: {}", n_in);

    //println!("m_in len {}", plink::bed_per_snv_size(n_in) * m_in_chrom);

    // for debug
    //println!("ys[0]:{:#010b}", ys[0]);

    //plink::check_valid_bed(fin_chrom, None, m_in_chrom, n_in).unwrap();

    // check if 32 GB? 16GB? remains or not

    // reading smaller than buf_size is not error but reason is unknown
    // [here](https://doc.rust-lang.org/std/io/trait.Read.html#tymethod.read)
    // -> to deal with this, every time use Seek to change position
    //let buf_size_limit: usize = 8 * 1024 * 1024 * 1024;
    let buf_size_limit: usize = 4 * 1024 * 1024 * 1024;
    let byte_per_snv = plink::bed_per_snv_size(n_in);
    let num_snv_per_buf: usize = buf_size_limit / byte_per_snv;
    let buf_num_snv: usize = num_snv_per_buf.min(m_in_chrom);
    let buf_size: usize = buf_num_snv * byte_per_snv;
    assert_eq!(buf_size % byte_per_snv, 0);
    assert!(buf_size <= buf_size_limit);
    //println!("buf_size {}", buf_size);

    //let mut reader = BufReader::new(File::open(plink::fname_bed(fin_chrom, None)).unwrap());
    let mut buf: Vec<B8_2> = vec![0; buf_size];

    let mut m_in_begin_loaded = 0;
    let mut m_begin_loaded = 0;
    // read buf length in one loop
    loop {
        //println!(
        //    "m_in_begin, m_begin {},{}",
        //    m_in_begin_loaded, m_begin_loaded
        //);

        // next start is m_in_begin_load
        let buf_next_begin = 3 + m_in_begin_loaded * byte_per_snv;
        //println!("buf_next_begin {}", buf_next_begin);
        reader.seek(SeekFrom::Start(buf_next_begin as u64)).unwrap();
        // This might be smaller than buf
        let loaded_byte = reader.read(&mut buf).unwrap();
        //if loaded_byte != buf_size {
        //    // last buffer or could not load whole buffer
        //    println!("loaded: {} {}", loaded_byte, buf_size);
        //}
        if loaded_byte == 0 {
            // all buf read
            break;
        } else {
            //println!("buf {:?}", buf.len());
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

            //println!("m_in_end, m_end {},{}", m_in_end_loaded, m_end_loaded);

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
}

/* // rayon impossible since cannot read inside
pub fn assign_predictions_file_loadpart(
    g_chrom: &mut GenotMut,
    fin_chrom: &str,
    use_snvs: &[bool],
    use_samples: &[bool],
) {
    let m_in_chrom = use_snvs.len();
    // println!("use_samples: {:?}", use_samples);

    // assume m and m_in in m_to_m_in are sorted
    // FIXME: use_snvs should be only chrom
    //let m_to_m_in = create_m_chrom_to_m_in_chrom(use_snvs);
    //let m_in_to_m = create_m_in_to_m_chrom(use_snvs, m_begin);

    let n_in = use_samples.len();

    // for debug
    //println!("ys[0]:{:#010b}", ys[0]);

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
                let bcode = plink::buf_to_ped_code(buf_mi, n_in_i);
                pred.set_bed_code_unchecked(bcode, ni);
                ni += 1;
            }
        }
    } else {
        for ni in 0..pred.n() {
            // bedcode
            let bcode = plink::buf_to_ped_code(buf_mi, ni);
            pred.set_bed_code_unchecked(bcode, ni);
        }
    }
}

// make missing to mode
// in pred
fn fill_missing(pred: &mut GenotSnvMut) {
    // count 0,1,2
    let mut counts_allele = Vec::with_capacity(4);
    for _ in 0..=3 {
        counts_allele.push(0);
    }

    let n = pred.n();
    //let len_n = len_n(n);
    // cannot borrowed late
    //let pred_s0m = &pred[..len_n_p];
    //let pred_s1m = &pred[len_n_p..];
    for ni in 0..n {
        counts_allele[pred.get_val_unchecked(ni) as usize] += 1;
        // counts_allele[predictions_both_to_count(pred, ni, len_n) as usize] += 1;
    }

    let mut mode: usize = 4;
    let mut mode_counts = 0;
    for i in 0..=2 {
        if counts_allele[i] > mode_counts {
            mode_counts = counts_allele[i];
            mode = i;
        }
    }
    let mode = mode as u8;
    assert_ne!(mode, 4);

    for ni in 0..n {
        if pred.get_val_unchecked(ni) == 3 {
            pred.set(mode, ni);
        }
    }
    // check all are non-missing?
    // -> performance...
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{io::Cursor, path::PathBuf};

    /// check if _file() and _file_loadpart() are the same

    fn setup_test() -> (PathBuf, Vec<bool>, usize, usize, Vec<bool>, Vec<bool>) {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let fin_snv = None;
        let fin_sample = None;

        let m_in: usize = plink::compute_num_snv(&fin).unwrap();
        println!("{}", m_in);
        let n_in: usize = plink::compute_num_sample(&fin).unwrap();
        println!("{}", n_in);
        // load snvs
        let snvs_in = plink::load_snvs(&fin, m_in);
        let (m, use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (n, use_samples) = plink::make_use_samples(fin_sample, &fin, n_in);
        let ys: Vec<bool> = plink::load_ys(&fin, n, &use_samples);

        (fin, ys, m, n, use_snvs, use_samples)
    }

    #[allow(dead_code)]
    fn setup_test2() -> (PathBuf, Vec<bool>, usize, usize, Vec<bool>, Vec<bool>) {
        //let fin = String::from("../../test/data/toy1/genot");
        let fin = PathBuf::from("../../test/data/1kg_n10000/1kg_n10000");
        let fin_snv = None;
        let fin_sample = None;

        let m_in: usize = plink::compute_num_snv(&fin).unwrap();
        println!("{}", m_in);
        let n_in: usize = plink::compute_num_sample(&fin).unwrap();
        println!("{}", n_in);
        // load snvs
        let snvs_in = plink::load_snvs(&fin, m_in);
        let (m, use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);
        //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
        let (n, use_samples) = plink::make_use_samples(fin_sample, &fin, n_in);
        let ys: Vec<bool> = plink::load_ys(&fin, n, &use_samples);

        (fin, ys, m, n, use_snvs, use_samples)
    }

    #[test]
    fn test_generate_predictions_whole() {
        let (fin, _, m, n, use_snvs, use_samples) = setup_test();
        //let use_snvs = vec![true; use_snvs.len()];
        //let use_samples = vec![true; use_samples.len()];
        let g = generate_genot(&fin, m, n, &use_snvs, Some(&use_samples), true);
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
        let mut g = Genot::new_empty(2, 3);

        let reader: Vec<u8> = vec![0x6c, 0x1b, 0x01, 0b00_01_11_00, 0b00_10_00_11];
        let mut cur = Cursor::new(reader.as_slice());
        assign_predictions(
            &mut g.as_genot_mut(),
            &mut cur,
            //reader.as_slice(),
            &[true, true],
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
}
