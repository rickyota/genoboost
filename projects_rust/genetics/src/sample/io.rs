use crate::io_genot::{self, GenotFormat};
use crate::textfile::{self};
use crate::vec;

use std::collections::HashMap;
use std::path::Path;

pub fn make_use_samples_buf(
    extract_sample_buf: Option<&[u8]>,
    fin: &Path,
    gfmt: GenotFormat,
) -> (usize, Vec<bool>) {
    let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();

    if extract_sample_buf.is_none() {
        let use_samples = vec![true; n_in];
        return (n_in, use_samples);
    }

    let samples_in: Vec<String> = io_genot::load_samples_id(fin, gfmt, None);
    //log::debug!("samples in plink {}", samples_in.len());

    // map: sample_in -> index
    let mut sample_in_to_index: HashMap<String, usize> = HashMap::with_capacity(n_in);
    for si in 0..samples_in.len() {
        //log::debug!("sample_in {:?}",
        //    samples_in[si].clone(),
        //);
        sample_in_to_index.insert(
            //samples::sample_id(samples_in[0][si].clone(), &samples_in[1][si]),
            samples_in[si].clone(),
            si,
        );
        //sample_in_to_index.insert(samples_in[0][si].clone() + ":" + &samples_in[1][si], si);
    }

    //textfile::check_open_file(fin_sample.unwrap());
    let samples_use: Vec<String> = load_samples_use_buf(extract_sample_buf.unwrap());
    //let samples_use: Vec<String> = load_samples_use(fin_sample.unwrap());
    //let samples_use: Vec<(String, String)> = load_samples_use(fin_sample.unwrap());
    //let n_use = samples_use.len();

    let mut use_samples = vec![false; n_in];
    for sample in samples_use.into_iter() {
        //log::debug!("sample {:?}",sample);
        //let sample_id = samples::sample_id(sample.0, &sample.1);
        if let Some(v) = sample_in_to_index.get(&sample) {
            use_samples[*v] = true;
        }
    }
    let n = vec::count_true(&use_samples);

    return (n, use_samples);
}

pub fn make_use_samples(
    fin_sample: Option<&Path>,
    fin: &Path,
    gfmt: GenotFormat,
) -> (usize, Vec<bool>) {
    let buf = fin_sample.map(|x| textfile::read_file_to_end(x, None).unwrap());

    make_use_samples_buf(buf.as_deref(), fin, gfmt)

    /*     let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();

    if fin_sample.is_none() {
        //for _ in 0..n_in {
        //    use_samples.push(true);
        //}
        let use_samples = vec![true; n_in];
        return (n_in, use_samples);
    }

    let samples_in: Vec<String> = io_genot::load_samples_id(fin, gfmt, None);

    // map: sample_in -> index
    let mut sample_in_to_index: HashMap<String, usize> = HashMap::with_capacity(n_in);
    for si in 0..samples_in.len() {
        sample_in_to_index.insert(
            //samples::sample_id(samples_in[0][si].clone(), &samples_in[1][si]),
            samples_in[si].clone(),
            si,
        );
        //sample_in_to_index.insert(samples_in[0][si].clone() + ":" + &samples_in[1][si], si);
    }

    //textfile::check_open_file(fin_sample.unwrap());
    let samples_use: Vec<String> = load_samples_use(fin_sample.unwrap());
    //let samples_use: Vec<(String, String)> = load_samples_use(fin_sample.unwrap());
    //let n_use = samples_use.len();

    let mut use_samples = vec![false; n_in];
    //let mut n: usize = 0;
    //for nui in 0..n_use {
    for sample in samples_use.into_iter() {
        //let sample_id = samples::sample_id(sample.0, &sample.1);
        //let sample_id: String = sample.0.clone() + ":" + &sample.1;
        if let Some(v) = sample_in_to_index.get(&sample) {
            use_samples[*v] = true;
            //n += 1;
        }
    }
    let n = vec::count_true(&use_samples);

    return (n, use_samples); */
}

/*
// does not assume FID==IID
pub fn make_use_samples(
    fin_sample: Option<&Path>,
    fin: &Path,
    gfmt: GenotFormat,
) -> (usize, Vec<bool>) {
    let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();

    if fin_sample.is_none() {
        //for _ in 0..n_in {
        //    use_samples.push(true);
        //}
        let use_samples = vec![true; n_in];
        return (n_in, use_samples);
    }

    textfile::check_open_file(fin_sample.unwrap());

    let mut use_samples = vec![false; n_in];
    // initialize with false
    //for _ in 0..n_in {
    //    use_samples.push(false);
    //}

    let fin_fam = io_genot::fname_fam_exist_chrom(fin,gfmt).unwrap();
    //let fin_fam = get_fname_fam(fin, Some(1));
    // fid, iid
    let cols = [0usize, 1];
    let samples_in: Vec<Vec<String>> = textfile::load_table_cols(&fin_fam, &cols, false).unwrap();
    //log::debug!("vals: {}", samples_in[0][0]);

    /*
    #[inline]
    fn sample_unique(fid: String, iid: &str) -> String {
        fid + ":" + iid
        //fid.clone() + ":" + iid
    }
    */

    // map: sample_in -> index
    let mut sample_in_to_index: HashMap<String, usize> = HashMap::with_capacity(n_in);
    for si in 0..samples_in[0].len() {
        sample_in_to_index.insert(
            samples::sample_id(samples_in[0][si].clone(), &samples_in[1][si]),
            si,
        );
        //sample_in_to_index.insert(samples_in[0][si].clone() + ":" + &samples_in[1][si], si);
    }

    let samples_use: Vec<(String, String)> = load_samples_use(fin_sample.unwrap());
    //let n_use = samples_use.len();

    let mut n: usize = 0;
    //for nui in 0..n_use {
    for sample in samples_use.into_iter() {
        let sample_id = samples::sample_id(sample.0, &sample.1);
        //let sample_id: String = sample.0.clone() + ":" + &sample.1;
        if let Some(v) = sample_in_to_index.get(&sample_id) {
            use_samples[*v] = true;
            n += 1;
        }
    }
}

*/

/// return vec of tuple of (fid, iid)
/// TODO: or just return Vec<Vec<String>>?
pub fn load_samples_use_buf(extract_sample_buf: &[u8]) -> Vec<String> {
    let cols = [0usize];
    let vss: Vec<Vec<String>> = textfile::load_table_cols_buf(extract_sample_buf, &cols, false);

    let mut samples: Vec<String> = vec![];
    for vi in 0..vss[0].len() {
        samples.push(vss[0][vi].clone());
        //samples.push((vss[0][vi].clone(), vss[1][vi].clone()));
    }
    samples
}

/// return vec of tuple of (fid, iid)
/// TODO: or just return Vec<Vec<String>>?
/// TODO: deprecate
pub fn load_samples_use(fin_sample: &Path) -> Vec<String> {
    let buf = textfile::read_file_to_end(fin_sample, None).unwrap();

    load_samples_use_buf(&buf[..])

    /*     //textfile::check_open_file(fin_sample);

    //let n_use = textfile::compute_num_line_text(fin_sample, None).unwrap();

    //let mut samples: Vec<(String, String)> = Vec::with_capacity(n_use);
    //let mut snvs: Vec<Snv> = Vec::with_capacity(m_use);

    //let fin_bim = get_fname_bim(fin_snv, Some(1));

    // load fid and iid
    //let cols = [0usize, 1];
    //let vss: Vec<Vec<String>> = textfile::load_table_cols(&fin_sample, &cols, false).unwrap();

    let cols = [0usize];
    let vss: Vec<Vec<String>> = textfile::load_table_cols(&fin_sample, &cols, false).unwrap();

    let mut samples: Vec<String> = vec![];
    for vi in 0..vss[0].len() {
        samples.push(vss[0][vi].clone());
        //samples.push((vss[0][vi].clone(), vss[1][vi].clone()));
    }
    samples */
}

pub fn vals_align_id(
    vals: &[String],
    ids: &[String],
    sample_id_to_n: &HashMap<String, usize>,
    //n: usize,
) -> Vec<String> {
    let n = sample_id_to_n.len();
    // TODO: better way?
    let mut vals_align: Vec<String> = vec![String::from(""); n];
    //let mut vals_align: Vec<String> = vec![String::from(""); n];

    for (n_in_i, val) in vals.iter().enumerate() {
        let sample_id = ids[n_in_i].clone();
        //let sample_id = samples::sample_id(ids.0[n_in_i].clone(), &ids.1[n_in_i]);

        if let Some(ni) = sample_id_to_n.get(&sample_id) {
            vals_align[*ni] = val.clone();
        }
        // else, sample id in .cov is not used.
    }

    // panic if any value is not assigned
    if vals_align.iter().any(|v| *v == "") {
        for n_in_i in 0..vals.len() {
            let sample_id = ids[n_in_i].clone();
            //let sample_id = samples::sample_id(ids.0[n_in_i].clone(), &ids.1[n_in_i]);

            if let Some(ni) = sample_id_to_n.get(&sample_id) {
                if vals_align[*ni] == "" {
                    panic!("Some sample is not found in fin_phe: {}.", sample_id);
                }
            }
        }
    }

    vals_align
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vals_align_id() {
        let vals: Vec<String> = ["1", "2", "3"].iter().map(|x| x.to_string()).collect();
        let ids: Vec<String> = ["id1", "id2", "id3"]
            .iter()
            .map(|x| x.to_string())
            .collect();
        let id_to_n: HashMap<String, usize> =
            HashMap::from([("id3".to_string(), 0), ("id1".to_string(), 1)]);

        let vals = vals_align_id(&vals, &ids, &id_to_n);
        assert_eq!(vals, vec!["3", "1"]);
    }
}
