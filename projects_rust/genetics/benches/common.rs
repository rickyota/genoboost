use std::path::Path;

use genetics::{alloc, Dataset, GenotFormat};

type B8 = u8;

#[allow(dead_code)]
pub fn predictions_len(m: usize, n: usize) -> usize {
    let len_n = n / 8 + 5;
    let len_pred = m * 2 * len_n;
    len_pred
}

#[allow(dead_code)]
pub fn setup_test_empty(n: usize, m: usize) -> Vec<B8> {
    let len_pred = predictions_len(m, n);

    let mut predictions: Vec<B8> = alloc::with_capacity_align_u8(len_pred);
    //let mut predictions: Vec<B8> = Vec::with_capacity(len_pred);
    for _ in 0..len_pred {
        predictions.push(0x00);
    }
    predictions
}

#[allow(dead_code)]
fn setup_vars(
    fin: &Path,
    fin_snv: Option<&Path>,
    fin_sample: Option<&Path>,
) -> (Dataset, Vec<f64>, Vec<f64>) {
    //let m_in: usize = plink::compute_num_snv(&fin).unwrap();
    //println!("m_in: {}", m_in);
    //let n_in: usize = plink::compute_num_sample(&fin).unwrap();
    //println!("n_in: {}", n_in);
    // load snvs
    //let snvs_in: Vec<SnvId> = plink::load_snvs(&fin, m_in);
    //println!("snvs_in {}", &snvs_in.len());
    //let (m, use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);
    //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
    //println!("use_snvs {}", &use_snvs.len());
    //let (n, use_samples) = plink::make_use_samples(fin_sample, &fin, n_in);
    //println!("m,n: {},{}", m, n);
    //let ys: Vec<B8> = plink::load_ys_b8(&fin, n, &use_samples);

    //(fin, ys, m, n, use_snvs, use_samples)

    //let predictions = genot_bi::generate_predictions(&fin, m, n, &use_snvs, &use_samples);

    let snv_buf = fin_snv.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());
    let sample_buf = fin_sample.map(|x| genetics::textfile::read_file_to_end(x, None).unwrap());

    let dataset: Dataset = Dataset::new(
        fin,
        GenotFormat::Plink1,
        None,
        None,
        "",
        snv_buf.as_deref(),
        //fin_snv,
        sample_buf.as_deref(),
        //fin_sample,
        false,
        None,
        None,
    );
    let n = dataset.samples().samples_n();
    let m = dataset.snvs().snvs_n();

    let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
    for _ in 0..n {
        ps.push(1.0f64 / (n as f64));
    }
    for _ in n..n + 32 {
        ps.push(0.0);
    }

    let mut losss = Vec::with_capacity(2 * m);
    for _ in 0..2 * m {
        losss.push(0.0f64);
    }

    (dataset, ps, losss)
}

/*
pub fn setup_test_n100k() -> (DatasetTwin, Vec<f64>, Vec<f64>) {
    let fin = String::from("../../test/data/1kg_n100000/1kg_n100000");
    //let fin = String::from("../../test/data/1kg_n100000/genot");
    // too large
    //let fin = String::from("../../test/boosting/data/1kg_n1000000/1kg_n1000000");
    //let fin = String::from("../../test/boosting/data/toy1/toy1");
    //let fin_snv = String::from("../../data/boosting/1000g/1kg_phase1_all/use_snvs_n100000.snvs");
    //let fin_snv = Some(fin_snv.as_str());
    // all snvs are too much
    let fin_snv = None;
    let fin_sample = None;

    setup_vars(&fin, fin_snv, fin_sample)

    /*
    let m_in: usize = plink::compute_num_snv(&fin).unwrap();
    println!("m_in: {}", m_in);
    let n_in: usize = plink::compute_num_sample(&fin).unwrap();
    println!("n_in: {}", n_in);
    // load snvs
    let snvs_in: Vec<Snv> = plink::load_snvs(&fin, m_in);
    println!("snvs_in {}", &snvs_in.len());
    let (m, use_snvs) = plink::make_use_snvs(fin_snv, &snvs_in);
    //let (m,use_snvs: Vec<bool>) = plink::make_use_snv(fin, snvs_in);
    println!("use_snvs {}", &use_snvs.len());
    let (n, use_samples) = plink::make_use_samples(fin_sample, &fin, n_in);
    println!("m,n: {},{}", m, n);
    let ys: Vec<B8> = plink::load_ys_b8(&fin, n, &use_samples);

    //(fin, ys, m, n, use_snvs, use_samples)

    let predictions = predict::generate_predictions(&fin, m, n, &use_snvs, &use_samples);

    let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
    for _ in 0..n {
        ps.push(1.0f64 / (n as f64));
    }
    for _ in n..n + 32 {
        ps.push(0.0);
    }

    let mut losss = Vec::with_capacity(2 * m);
    for _ in 0..2 * m {
        losss.push(0.0f64);
    }

    (predictions, ys, m, n, ps, losss)
    */
}

pub fn setup_test_n10k() -> (DatasetTwin, Vec<f64>, Vec<f64>) {
    let fin = String::from("../../test/data/1kg_n10000/1kg_n10000");
    //let fin = String::from("../../test/boosting/data/1kg_n10000/genot");
    // too large
    //let fin = String::from("../../test/boosting/data/1kg_n1000000/1kg_n1000000");
    let fin_snv = None;
    let fin_sample = None;

    setup_vars(&fin, fin_snv, fin_sample)
}
 */
