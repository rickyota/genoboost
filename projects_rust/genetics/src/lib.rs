// change arg name from v to ar, and for v,i

pub mod alloc;
pub mod cov;
pub mod dataset;
pub mod dataset_file;
pub mod dout_file;
//pub mod genet_path;
pub mod genot;
pub mod genot_calc;
pub mod genot_io;
pub mod genotype;
pub mod param;
pub mod pgs;
pub mod regression;
pub mod sample;
pub mod score;
pub mod snv;
pub mod sum_stat;
pub mod textfile;
pub mod vec;
pub mod wgt;
mod wgts;

use std::collections::HashMap;

// you can call Snv with `use io_rust::Snv` instead of `use io_rust::snv::Snv`
// should call with `use io_rust::wgt::Model` to avoid confusion
// -> no meaning
//pub use wgt::{CovWgt, Model, SnvWgt, Wgt, WgtKind};
// only make pub of Wgt
pub use cov::CovKind;
pub use dataset::Dataset;
pub use dataset_file::DatasetFile;
pub use dout_file::{DoutFile, DoutScoreFile, WgtDoutOrFile};
pub use genot::genot_struct;
pub use genot::genot_struct::{Genot, B8, B8_2, THRESHOLD_SNV};
pub use genot::prelude::*; // for test_pgenlib
pub use genot_io::GenotFile;
pub use param::LdCriteria;
pub use sample::covs::{CovId, Covs, CovsTrait};
pub use sample::samples;
pub use sample::samples::Samples;
pub use score::{SampleScore, Sum3by3};
pub use snv::snvs;
pub use snv::Snvs;
pub use snv::{Chrom, SnvId};
use std::path::Path;
use std::path::PathBuf;
pub use sum_stat::SumStat;
pub use wgt::{Model, SnvInteractionWgt, SnvWgt, Wgt, WgtKind, WgtTrait};
pub use wgts::Wgts;

// for assert_float_absolute_eq
#[macro_use]
extern crate assert_float_eq;

//pub fn

// TODO: create and use dout_file.rs
pub fn run_score(
    dout_score: &DoutScoreFile,
    //dout_score: &Path,
    dfile: &DatasetFile,
    //fin_genot: &GenotFile,
    //fin_phe: Option<&Path>,
    //cov_name: Option<&str>,
    wgt_d_f: &WgtDoutOrFile,
    //dout_wgt: Option<&Path>,
    //fout_wgt: Option<&Path>,
    //fin_sample: Option<&Path>,
    concat: Option<&str>,
    no_concat: Option<&str>,
    is_resume: bool,
    fill_missing_in_dataset: bool,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    missing_to_mode: bool,
    missing_to_mean: bool,
    is_nonadd: bool,
    mem: Option<usize>,
) {
    //fin_genot.check_valid_open();

    //if let Some(dout_wgt) = dout_wgt {
    match wgt_d_f {
        WgtDoutOrFile::Dout(dout_wgt) => {
            let files_wgt_ori = dout_wgt.get_files_wgt();
            //let files_wgt_ori = wgt::io::get_files_wgt(dout_wgt);

            if files_wgt_ori.len() == 0 {
                panic!("No wgt files exist in dir: {:?}.", dout_wgt);
            }

            let files_wgt = if let Some(concat_para) = concat {
                if is_resume {
                    log::debug!("--is-resume but --concat is set. Ignore --is-resume.");
                }
                //files_wgt_ori
                // only use files with concat_para
                files_wgt_ori
                    .iter()
                    .filter(|fwgt| dout_file::is_fwgt_concat(fwgt, concat_para))
                    .map(|x| x.clone())
                    .collect::<Vec<PathBuf>>()
                //.filter(|fwgt| genet_path::is_fwgt_concat(fwgt, concat_para))
            } else if let Some(concat_para) = no_concat {
                if is_resume {
                    log::debug!("--is-resume but --no-concat is set. Ignore --is-resume.");
                }
                // only use files without concat_para
                files_wgt_ori
                    .iter()
                    .filter(|fwgt| !dout_file::is_fwgt_concat(fwgt, concat_para))
                    .map(|x| x.clone())
                    .collect::<Vec<PathBuf>>()
                //.filter(|fwgt| !genet_path::is_fwgt_concat(fwgt, concat_para))
            } else {
                // concat.is_none() and no_concat.is_none()
                if is_resume {
                    // exclude fwgt if .score already exist
                    let mut v = Vec::new();
                    for fwgt in files_wgt_ori {
                        let fscore = dout_score.fscore_from_fwgt(&fwgt);
                        //let fscore = genet_path::fscore_from_fwgt(dout_score, &fwgt);
                        if !textfile::exist_file(&fscore) {
                            v.push(fwgt);
                        }
                    }
                    v
                } else {
                    files_wgt_ori
                }
            };

            log::debug!("wgt files: {:?}", files_wgt);

            if files_wgt.len() == 0 {
                log::info!("All scores exist.");
                return;
            }

            let mut wgts_vec: Vec<Wgts> = files_wgt
                .iter()
                .map(|file_wgt| Wgts::new_from_file(file_wgt, is_nonadd))
                .collect();
            //.map(|file_wgt| Wgts::new_from_file(file_wgt, use_snv_pos, is_nonadd))

            let dataset = Dataset::new_datasetfile_score_genetics(
                dfile,
                //fin_genot,
                //phe_buf.as_deref(),
                //cov_name,
                //sample_buf.as_deref(),
                &mut wgts_vec,
                fill_missing_in_dataset,
                allow_nonexist_snv,
                use_snv_pos,
                mem,
            );

            //let phe_buf = fin_phe.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
            //let sample_buf = fin_sample.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());

            //// TODO: use dfile
            //let dataset = Dataset::new_score_genetics(
            //    fin_genot,
            //    phe_buf.as_deref(),
            //    cov_name,
            //    sample_buf.as_deref(),
            //    &mut wgts_vec,
            //    fill_missing_in_test,
            //    allow_nonexist_snv,
            //    use_snv_pos,
            //    mem,
            //);

            //if let Some(concat_para) = concat {
            //    for fwgt in files_wgt.iter() {
            //        if !(fwgt
            //            .file_stem()
            //            .unwrap()
            //            .to_str()
            //            .unwrap()
            //            .contains(&("_".to_string() + concat_para + "-")))
            //        {
            //            panic!("Parameter not in file name. Avoid using --concat.");
            //        }

            //        if fwgt
            //            .file_stem()
            //            .unwrap()
            //            .to_str()
            //            .unwrap()
            //            .split(&("_".to_string() + concat_para + "-"))
            //            .collect::<Vec<&str>>()
            //            .last()
            //            .unwrap()
            //            .contains("_")
            //        {
            //            panic!("File stem name should end with _(para)-*. ");
            //        }
            //    }
            //}

            let samples_id = dataset.samples().names();

            let mut score_paras: Vec<SampleScore> = vec![];
            //let mut score_paras: Vec<Vec<f64>> = vec![];

            // rayon in score::calculate_score
            if concat.is_some() {
                for (wgts, file_wgt) in wgts_vec.iter().zip(files_wgt.iter()) {
                    log::debug!("file_wgt: {:?}", file_wgt);
                    let score_v = score::calc_score_wgts(
                        &wgts,
                        &dataset,
                        allow_nonexist_snv,
                        missing_to_mode,
                        missing_to_mean,
                    );
                    score_paras.push(score_v);
                }

                //let mut score_paras: Vec<Vec<f64>> = vec![vec![]; wgts_vec.len()];
                //score_paras_for_concat
                //    .iter_mut()
                //    .zip(wgts_vec.iter())
                //    .enumerate()
                //    .for_each(|(wgt_i, (score_para, wgts))| {
                //        log::debug!("wgt_i: {:?}", wgt_i);
                //        let score_v = score::calculate_score(&wgts, &dataset);
                //        // error
                //        score_para.copy_from_slice(&score_v);
                //    });
            } else {
                wgts_vec
                    .iter()
                    .zip(files_wgt.iter())
                    .for_each(|(wgts, file_wgt)| {
                        log::debug!("file_wgt: {:?}", file_wgt);
                        let score_v = score::calc_score_wgts(
                            &wgts,
                            &dataset,
                            allow_nonexist_snv,
                            missing_to_mode,
                            missing_to_mean,
                        );
                        // --no-concat or no args
                        let fscore = dout_score.fscore_from_fwgt(file_wgt);
                        //let fout_score = genet_path::fscore_from_fwgt(dout_score, file_wgt);
                        score::write_scores(&fscore, &score_v, samples_id);
                    });
                // score_paras is empty
            }

            // Not Good
            // rayon here!
            // Speed up when several wgts are input.
            // -> TODO: much better if able to use rayon in score::calculate_score
            //if concat.is_some() {
            //    let mut score_paras_for_concat: Vec<Vec<f64>> = vec![vec![]; wgts_vec.len()];

            //    score_paras_for_concat
            //        .par_iter_mut()
            //        .enumerate()
            //        .for_each(|(wgt_i, score_para)| {
            //            log::debug!("wgt_i: {:?}", wgt_i);
            //            let score_v = score::calculate_score(&wgts_vec[wgt_i], &dataset);
            //            score_para.copy_from_slice(&score_v);
            //        });

            //    score_paras = score_paras_for_concat;
            //} else {
            //    wgts_vec
            //        .par_iter()
            //        .zip(files_wgt.par_iter())
            //        .for_each(|(wgts, file_wgt)| {
            //            log::debug!("file_wgt: {:?}", file_wgt);
            //            let score_v = score::calculate_score(&wgts, &dataset);
            //            // --no-concat or no args
            //            let fout_score = genet_path::fscore_from_fwgt(dout_score, file_wgt);
            //            score::write_scores(&fout_score, &score_v, samples_id);
            //        });
            //    // score_paras is empty
            //}

            //for (wgts, file_wgt) in wgts_vec.iter().zip(files_wgt.iter()) {
            //    log::debug!("file_wgt: {:?}", file_wgt);

            //    //score::score(&fout_score, &wgts, &dataset, &samples_id);
            //    let score_v = score::calculate_score(&wgts, &dataset);

            //    if concat.is_some() {
            //        score_paras.push(score_v);
            //    } else {
            //        // --no-concat or no args
            //        let fout_score = genet_path::fscore_from_fwgt(dout_score, file_wgt);
            //        score::write_scores(&fout_score, &score_v, samples_id)
            //    }
            //}

            if let Some(concat_para) = concat {
                // grouping scores with the same fout_concat
                // ex. ['clump_p-0.1_n-100', 'clump_p-0.1_n-200'] -> 'clump_p-0.1_n.score'
                //  ['clump_p-0.2_n-100', 'clump_p-0.2_n-200'] -> 'clump_p-0.2_n.score'

                let fouts_concat: Vec<PathBuf> = files_wgt
                    .iter()
                    .map(|fwgt| dout_score.fscore_concat(concat_para, fwgt))
                    .collect();
                //.map(|x| genet_path::fscore_concat(dout_score, concat_para, x))
                //let fout_concat = fscore_concat(dout_score, concat_para, &files_wgt[0]);
                //log::debug!("write to: {:?}", fout_concat);

                let paras: Vec<String> = files_wgt
                    .iter()
                    .map(|x| dout_file::para_from_fwgt(x, concat_para))
                    .collect();
                //.map(|x| genet_path::para_from_fwgt(x, concat_para))

                let foutd_paras: HashMap<PathBuf, Vec<String>> = fouts_concat
                    .iter()
                    .zip(paras.iter())
                    .fold(HashMap::new(), |mut paras_each, (fout, para)| {
                        paras_each
                            .entry(fout.clone())
                            .or_insert(vec![])
                            .push(para.clone());
                        paras_each
                    });

                let foutd_scores: HashMap<PathBuf, Vec<SampleScore>> = fouts_concat
                    .iter()
                    .zip(score_paras.into_iter())
                    .fold(HashMap::new(), |mut scores_each, (fout, score)| {
                        scores_each
                            .entry(fout.clone())
                            .or_insert(vec![])
                            .push(score);
                        scores_each
                    });

                //let foutd_scores: HashMap<PathBuf, Vec<Vec<f64>>> = fouts_concat
                //    .iter()
                //    .zip(score_paras.iter())
                //    .fold(HashMap::new(), |mut scores_each, (fout, score)| {
                //        scores_each
                //            .entry(fout.clone())
                //            .or_insert(vec![])
                //            .push(score.clone());
                //        scores_each
                //    });

                for (fout_concat, paras) in foutd_paras.iter() {
                    let score_para = foutd_scores.get(fout_concat).unwrap();

                    score::write_scores_paras(
                        &fout_concat,
                        &score_para,
                        concat_para,
                        &paras,
                        samples_id,
                    )
                }
                //score::write_scores_paras_nopheno(
                //    &fout_concat,
                //    &score_paras,
                //    concat_para,
                //    &paras,
                //    samples_id,
                //)
            }
        }
        //} else if let Some(file_wgt) = fout_wgt {
        WgtDoutOrFile::File(file_wgt) => {
            wgt::io::check_file_wgt_exist(&file_wgt);
            let wgts = Wgts::new_from_file(&file_wgt, is_nonadd);
            //let wgts = Wgts::new_from_file(&file_wgt, use_snv_pos, is_nonadd);
            let mut wgts_vec = vec![wgts];

            let dataset = Dataset::new_datasetfile_score_genetics(
                dfile,
                //fin_genot,
                //phe_buf.as_deref(),
                //cov_name,
                //sample_buf.as_deref(),
                &mut wgts_vec,
                fill_missing_in_dataset,
                allow_nonexist_snv,
                use_snv_pos,
                mem,
            );

            //let phe_buf = fin_phe.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
            //let sample_buf = fin_sample.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
            //let dataset = Dataset::new_score_genetics(
            //    fin_genot,
            //    phe_buf.as_deref(),
            //    cov_name,
            //    sample_buf.as_deref(),
            //    &mut wgts_vec,
            //    fill_missing_in_test,
            //    allow_nonexist_snv,
            //    use_snv_pos,
            //    mem,
            //);

            let fout_score = dout_score.fscore_from_fwgt(file_wgt);
            //let fout_score = genet_path::fscore_from_fwgt(dout_score, file_wgt);
            //let fname_score = file_wgt.file_name().unwrap().to_str().unwrap().to_owned() + ".score";
            //let fout_score = dout_score.join(fname_score);

            //score::score(&fout_score, &wgts_vec[0], &dataset, &samples_id);
            let score_v = score::calc_score_wgts(
                &wgts_vec[0],
                &dataset,
                allow_nonexist_snv,
                missing_to_mode,
                missing_to_mean,
            );
            let samples_id = dataset.samples().names();
            score::write_scores(&fout_score, &score_v, samples_id)
        }
    }
}

pub fn test_pgenlib(fin_genot: &GenotFile, _fin_sample: Option<&Path>) -> Genot {
    let snvs_in = genot_io::load_snvs(fin_genot);

    for (si, snv) in snvs_in.iter().enumerate() {
        if si < 10 {
            println!("snv {:?}", snv);
            //log::debug!("snv {:?}", snv);
        }
    }

    let m_in: usize = genot_io::compute_num_snv(fin_genot).unwrap();
    log::debug!("m_in: {}", m_in);
    let n_in: usize = genot_io::compute_num_sample(fin_genot).unwrap();
    log::debug!("n_in: {}", n_in);

    let use_snvs = if m_in < 100 {
        //let use_snvs = vec![true; m_in];
        vec![true; m_in]
    } else {
        //let mut u_ = vec![true; 100];
        //u_.extend(vec![false; m_in - 100]);
        let mut u_ = vec![false; m_in - 100];
        u_.extend(vec![true; 100]);
        u_
    };

    println!("to generate");
    let g =
        genot_io::load::generate_genot(fin_genot, m_in, n_in, Some(&use_snvs), None, false, None);
    println!("generated");

    for (si, snv) in g.iter_snv().enumerate() {
        if si < 10 {
            //log::debug!("genot {:?}", &snv.vals()[..10]);
            println!("genot {:?}", &snv.vals()[..10]);
        }
    }
    g
}
