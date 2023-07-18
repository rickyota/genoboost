// change arg name from v to ar, and for v,i

pub mod alloc;
pub mod cov;
pub mod dataset;
pub mod genot;
pub mod genot_index;
pub mod pgs;
pub mod regression;
pub mod sample;
mod score;
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
pub use dataset::samples;
pub use dataset::samples::{CovId, Covs, Samples};
pub use dataset::snvs;
pub use dataset::snvs::Snvs;
pub use dataset::Dataset;
pub use genot::genot_struct;
pub use genot::genot_struct::{Genot, B8, B8_2, THRESHOLD_SNV};
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
    let fname_score = fwgt.file_name().unwrap().to_str().unwrap().to_owned() + ".score";
    let fout_score = dscore.join(fname_score);
    fout_score
}

pub fn run_score(
    dout_score: &Path,
    fin: &Path,
    gfmt: GenotFormat,
    fin_phe: Option<&Path>,
    phe_name: Option<&str>,
    cov_name: Option<&str>,
    dout_wgt: Option<&Path>,
    fout_wgt: Option<&Path>,
    //fin_cov: Option<&Path>,
    fin_sample: Option<&Path>,
    is_resume: bool,
) {
    io_genot::check_valid_fin(fin, gfmt);

    //let n_in: usize = io_genot::compute_num_sample(fin, gfmt).unwrap();
    let (_, use_samples) = sample::make_use_samples(fin_sample, fin, gfmt);
    let samples_id = io_genot::load_samples_id(fin, gfmt, Some(&use_samples));

    if let Some(dout_wgt) = dout_wgt {
        let files_wgt_ori = wgt::io::get_files_wgt(dout_wgt);

        if files_wgt_ori.len() == 0 {
            panic!("No wgt files exist.");
        }

        let files_wgt = if is_resume {
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
        };

        log::debug!("wgt files: {:?}", files_wgt);

        if files_wgt.len() == 0 {
            log::info!("All scores exist.");
            return;
        }

        let mut wgts_vec: Vec<Wgts> = files_wgt
            .iter()
            .map(|file_wgt| Wgts::new_from_file(file_wgt))
            .collect();

        let phe_buf = fin_phe.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
        let sample_buf = fin_sample.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());

        let dataset = Dataset::new_score_genetics(
            fin,
            gfmt,
            phe_buf.as_deref(),
            //fin_phe,
            phe_name,
            cov_name,
            sample_buf.as_deref(),
            //fin_sample,
            &mut wgts_vec,
        );

        for (wgts, file_wgt) in wgts_vec.iter().zip(files_wgt.iter()) {
            log::debug!("file_wgt: {:?}", file_wgt);

            let fout_score = fscore_from_fwgt(dout_score, file_wgt);
            //let fname_score = file_wgt.file_name().unwrap().to_str().unwrap().to_owned() + ".score";
            //let fout_score = dout_score.join(fname_score);
            score::score(&fout_score, &wgts, &dataset, &samples_id);
        }
    } else if let Some(file_wgt) = fout_wgt {
        wgt::io::check_file_wgt_exist(&file_wgt);
        let wgts = Wgts::new_from_file(&file_wgt);
        let mut wgts_vec = vec![wgts];

        let phe_buf = fin_phe.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
        let sample_buf = fin_sample.map(|x| crate::textfile::read_file_to_end(x, None).unwrap());
        let dataset = Dataset::new_score_genetics(
            fin,
            gfmt,
            phe_buf.as_deref(),
            //fin_phe,
            phe_name,
            cov_name,
            sample_buf.as_deref(),
            //fin_sample,
            &mut wgts_vec,
        );
        let fout_score = fscore_from_fwgt(dout_score, file_wgt);
        //let fname_score = file_wgt.file_name().unwrap().to_str().unwrap().to_owned() + ".score";
        //let fout_score = dout_score.join(fname_score);

        score::score(&fout_score, &wgts_vec[0], &dataset, &samples_id);
    } else {
        panic!("sth wrong.")
    }
}
