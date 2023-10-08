// change arg name from v to ar, and for v,i

use std::collections::HashMap;

pub mod alloc;
pub mod cov;
pub mod dataset;
pub mod genot;
pub mod genot_index;
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
//pub mod io_genot;
//pub mod samples;
//pub mod snvs;
//pub mod tsv;

// you can call Snv with `use io_rust::Snv` instead of `use io_rust::snv::Snv`
// should call with `use io_rust::wgt::Model` to avoid confusion
// -> no meaning
//pub use wgt::{CovWgt, Model, SnvWgt, Wgt, WgtKind};
// only make pub of Wgt
pub use cov::{CovKind, Var};
pub use dataset::io_genot;
pub use dataset::io_genot::load;
pub use dataset::samples;
pub use dataset::samples::{CovId, Covs, Samples};
pub use dataset::snvs;
pub use dataset::snvs::Snvs;
pub use dataset::Dataset;
pub use genot::genot_struct;
pub use genot::genot_struct::{Genot, B8, B8_2, THRESHOLD_SNV};
pub use genot::prelude::*;
pub use io_genot::GenotFormat;
pub use snv::{Chrom, SnvId};
use std::path::Path;
use std::path::PathBuf;
pub use sum_stat::SumStat;
pub use wgt::{CovWgt, Model, ModelType, SnvWgt, Wgt, WgtKind, WgtTrait};
use wgts::Wgts;

#[macro_use]
extern crate assert_float_eq;

// TODO: better
fn fscore_from_fwgt(dscore: &Path, fwgt: &Path) -> PathBuf {
    let fname_score = fwgt.file_stem().unwrap().to_str().unwrap().to_owned() + ".score";
    //let fname_score = fwgt.file_name().unwrap().to_str().unwrap().to_owned() + ".score";
    let fout_score = dscore.join(fname_score);
    fout_score
}

fn fscore_concat(dscore: &Path, para: &str, fwgt: &Path) -> PathBuf {
    // fwgt to get method name
    // ex. 'clump_p-0.1_n-100' -> 'clump_p-0.1_n.score'

    let method = fwgt
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .split(&("_".to_string() + para + "-"))
        .collect::<Vec<&str>>()[0]
        .to_string();

    // [method]_n.score
    //let method = fwgt
    //    .file_name()
    //    .unwrap()
    //    .to_str()
    //    .unwrap()
    //    .split("_")
    //    .collect::<Vec<&str>>()[0]
    //    .to_string();
    let fname_score = method + "_" + para + ".score";
    let fout_score = dscore.join(fname_score);
    fout_score
}

fn para_from_fwgt(fwgt: &Path, para: &str) -> String {
    let para = fwgt
        .file_stem() // exclude .wgt here
        .unwrap()
        .to_str()
        .unwrap()
        .split(&("_".to_string() + para + "-"))
        .collect::<Vec<&str>>()[1]
        .to_string();
    para
}

fn is_fwgt_concat(fwgt: &Path, concat_para: &str) -> bool {
    let is_concat = fwgt
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .contains(&("_".to_string() + concat_para + "-"));

    if is_concat {
        // check if file stem name end with _(para)-*
        if fwgt
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
            .split(&("_".to_string() + concat_para + "-"))
            .collect::<Vec<&str>>()
            .last()
            .unwrap()
            .contains("_")
        {
            panic!("File stem name should end with _(para)-*. ");
        }
    };
    is_concat
}

pub fn run_score(
    dout_score: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    fin_phe: Option<&Path>,
    //phe_name: Option<&str>,
    cov_name: Option<&str>,
    dout_wgt: Option<&Path>,
    fout_wgt: Option<&Path>,
    //fin_cov: Option<&Path>,
    fin_sample: Option<&Path>,
    concat: Option<&str>,
    no_concat: Option<&str>,
    is_resume: bool,
    allow_nonexist_snv: bool,
    use_snv_pos: bool,
    is_nonadd: bool,
) {
    io_genot::check_valid_fin(fin, gfmt);

    //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
    //let (_, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
    //let samples_id = io_genot::load_samples_id(fin, gfmt, Some(&use_samples));

    if let Some(dout_wgt) = dout_wgt {
        let files_wgt_ori = wgt::io::get_files_wgt(dout_wgt);

        if files_wgt_ori.len() == 0 {
            panic!("No wgt files exist.");
        }

        let files_wgt = if let Some(concat_para) = concat {
            if is_resume {
                log::debug!("--is-resume but --concat is set. Ignore --is-resume.");
            }
            //files_wgt_ori
            // only use files with concat_para
            files_wgt_ori
                .iter()
                .filter(|fwgt| is_fwgt_concat(fwgt, concat_para))
                .map(|x| x.clone())
                .collect::<Vec<PathBuf>>()
        } else if let Some(concat_para) = no_concat {
            if is_resume {
                log::debug!("--is-resume but --no-concat is set. Ignore --is-resume.");
            }
            // only use files without concat_para
            files_wgt_ori
                .iter()
                .filter(|fwgt| !is_fwgt_concat(fwgt, concat_para))
                .map(|x| x.clone())
                .collect::<Vec<PathBuf>>()
        } else {
            // concat.is_none() and no_concat.is_none()
            if is_resume {
                // exclude fwgt if .score already exist
                let mut v = Vec::new();
                for fwgt in files_wgt_ori {
                    let fscore = fscore_from_fwgt(dout_score, &fwgt);
                    if !textfile::exist_file(&fscore) {
                        v.push(fwgt);
                    }
                }
                v
            } else {
                files_wgt_ori
            }
        };

        //let files_wgt = if is_resume & concat.is_none() {
        //    // exclude fwgt if .score already exist
        //    let mut v = Vec::new();
        //    for fwgt in files_wgt_ori {
        //        let fscore = fscore_from_fwgt(dout_score, &fwgt);
        //        if !textfile::exist_file(&fscore) {
        //            v.push(fwgt);
        //        }
        //    }
        //    v
        //} else {
        //    if is_resume & concat.is_some() {
        //        log::debug!("--is-resume but --concat is set. Ignore --is-resume.");
        //    }
        //    files_wgt_ori
        //};

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

        let phe_buf = fin_phe.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
        let sample_buf = fin_sample.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());

        let dataset = Dataset::new_score_genetics(
            fin,
            gfmt,
            phe_buf.as_deref(),
            //fin_phe,
            //phe_name,
            cov_name,
            sample_buf.as_deref(),
            //fin_sample,
            &mut wgts_vec,
            allow_nonexist_snv,
            use_snv_pos,
        );

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
        let mut score_paras: Vec<Vec<f64>> = vec![];
        for (wgts, file_wgt) in wgts_vec.iter().zip(files_wgt.iter()) {
            log::debug!("file_wgt: {:?}", file_wgt);

            //score::score(&fout_score, &wgts, &dataset, &samples_id);
            let score_v = score::score_nowrite(&wgts, &dataset);

            if concat.is_some() {
                score_paras.push(score_v);
            } else {
                // --no-concat or no args
                let fout_score = fscore_from_fwgt(dout_score, file_wgt);
                score::write_scores_nopheno(&fout_score, &score_v, samples_id)
            }
        }

        if let Some(concat_para) = concat {
            // grouping scores with the same fout_concat
            // ex. ['clump_p-0.1_n-100', 'clump_p-0.1_n-200'] -> 'clump_p-0.1_n.score'
            //  ['clump_p-0.2_n-100', 'clump_p-0.2_n-200'] -> 'clump_p-0.2_n.score'

            let fouts_concat: Vec<PathBuf> = files_wgt
                .iter()
                .map(|x| fscore_concat(dout_score, concat_para, x))
                .collect();
            //let fout_concat = fscore_concat(dout_score, concat_para, &files_wgt[0]);
            //log::debug!("write to: {:?}", fout_concat);

            let paras: Vec<String> = files_wgt
                .iter()
                .map(|x| para_from_fwgt(x, concat_para))
                .collect();

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

            let foutd_scores: HashMap<PathBuf, Vec<Vec<f64>>> = fouts_concat
                .iter()
                .zip(score_paras.iter())
                .fold(HashMap::new(), |mut scores_each, (fout, score)| {
                    scores_each
                        .entry(fout.clone())
                        .or_insert(vec![])
                        .push(score.clone());
                    scores_each
                });

            for (fout_concat, paras) in foutd_paras.iter() {
                let score_para = foutd_scores.get(fout_concat).unwrap();

                score::write_scores_paras_nopheno(
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
    } else if let Some(file_wgt) = fout_wgt {
        wgt::io::check_file_wgt_exist(&file_wgt);
        let wgts = Wgts::new_from_file(&file_wgt, is_nonadd);
        //let wgts = Wgts::new_from_file(&file_wgt, use_snv_pos, is_nonadd);
        let mut wgts_vec = vec![wgts];

        let phe_buf = fin_phe.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
        let sample_buf = fin_sample.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
        let dataset = Dataset::new_score_genetics(
            fin,
            gfmt,
            phe_buf.as_deref(),
            //fin_phe,
            //phe_name,
            cov_name,
            sample_buf.as_deref(),
            //fin_sample,
            &mut wgts_vec,
            allow_nonexist_snv,
            use_snv_pos,
        );
        let fout_score = fscore_from_fwgt(dout_score, file_wgt);
        //let fname_score = file_wgt.file_name().unwrap().to_str().unwrap().to_owned() + ".score";
        //let fout_score = dout_score.join(fname_score);

        //score::score(&fout_score, &wgts_vec[0], &dataset, &samples_id);
        let score_v = score::score_nowrite(&wgts_vec[0], &dataset);
        let samples_id = dataset.samples().names();
        score::write_scores_nopheno(&fout_score, &score_v, samples_id)
    } else {
        panic!("sth wrong.")
    }
}

pub fn test_pgenlib(fin: &Path, gfmt: GenotFormat, _fin_sample: Option<&Path>) -> Genot {
    let snvs_in = io_genot::load_snvs(fin, gfmt);

    for (si, snv) in snvs_in.iter().enumerate() {
        if si < 10 {
            println!("snv {:?}", snv);
            //log::debug!("snv {:?}", snv);
        }
    }

    let m_in: usize = io_genot::compute_num_snv(fin, gfmt).unwrap();
    log::debug!("m_in: {}", m_in);
    let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
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
    let g = load::generate_genot(fin, gfmt, m_in, n_in, Some(&use_snvs), None, false);
    println!("generated");

    for (si, snv) in g.iter_snv().enumerate() {
        if si < 10 {
            //log::debug!("genot {:?}", &snv.vals()[..10]);
            println!("genot {:?}", &snv.vals()[..10]);
        }
    }
    g
}
