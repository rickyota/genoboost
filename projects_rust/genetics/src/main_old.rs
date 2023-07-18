//! General lib and application of general functions
//! - Polygenic score
//! 
//! 
//! 

#[macro_use]
extern crate clap;
use clap::{AppSettings, Arg,  ArgMatches, SubCommand};
//use rayon;
use std::{path::PathBuf, time::Instant};
//use clap::app_from_crate; //error

fn check_usize(v: String) -> Result<(), String> {
    match v.parse::<usize>() {
        Ok(_) => Ok(()),
        Err(_) => {
            Err("The value should be able to be parsed into non-negative integer.".to_owned())
        }
    }
}

fn get_matches() -> ArgMatches<'static> {
    let matches = app_from_crate!()
        .setting(AppSettings::DeriveDisplayOrder)
        .setting(AppSettings::SubcommandRequiredElseHelp) //require subcommand or show help
        .arg(Arg::from_usage("--num_threads [NUM] 'Number of threads.'").validator(check_usize))
        .arg(Arg::from_usage("--verbose 'Verbose.'"))
        .subcommand( SubCommand::with_name("score").about("Calculate score.")
            .arg(Arg::from_usage("--dir_score <FILE> 'Directory of output file.'"))
            .arg(Arg::from_usage("--file_plink <FILE> 'Prefix of a plink file. If those files are separated into chromosomes, set '%' where chromosome number is.'"))
            .arg(Arg::from_usage("--file_phe [FILE] 'File of phenotype to use.'"))
            .arg(Arg::from_usage("--phe [VAL] 'Phonotype. must use with --file_phe.'").default_value("None")) // TODO: validator
            .arg(Arg::from_usage("--file_sample [FILE] 'File of samples to use.'"))
            // maybe later
            //.arg(Arg::from_usage("--file_cov [FILE] 'File of covariates to use.'"))
            .arg(Arg::from_usage("--dir_wgt [FILE] 'Directory of a weight file. Must use either --dir_wgt or --file_wgt.'"))
            .arg(Arg::from_usage("--file_wgt [FILE] 'Prefix of a weight file. Must use either --dir_wgt or --file_wgt.'"))
            .arg(Arg::from_usage("--resume 'Do not overwrite existing score.'"))
        )
        .get_matches();
    matches
}


fn main() {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            log::info!("Able to use SIMD.")
        } else {
            log::info!("Not able to use SIMD since avx2 is not detected.")
        }
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        log::info!("Not able to use SIMD since arch is not x86 or x86_64.")
    }

    let matches = get_matches();

    let verbose = matches.is_present("verbose");
    if verbose {
        std::env::set_var("RUST_LOG", "debug");
    } else {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::init();

    log::debug!("matches {:?}", matches);

    let start = Instant::now();

    if let Some(n_threads) = matches.value_of("num_threads") {
        let n_threads = n_threads.parse::<usize>().unwrap();
        //n_threads=n_threads.min(num_cpus::get());
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .unwrap();
    };
    // otherwise, use default thread number
    log::debug!("num_thread set: {}", rayon::current_num_threads());



    if let Some(ref matches) = matches.subcommand_matches("score") {
        let dout_score = PathBuf::from(matches.value_of("dir_score").unwrap());
        let fin = PathBuf::from(matches.value_of("file_plink").unwrap());
        let fin_phe = matches.value_of("file_phe").map(|x| PathBuf::from(x));
        let phe_name = match matches.value_of("phe").unwrap() {
            "None" => None,
            z => Some(z),
        };
        let fin_sample = matches.value_of("file_sample").map(|x| PathBuf::from(x));
        //let fin_cov = matches.value_of("file_cov").map(|x| PathBuf::from(x));
        let is_resume = matches.is_present("resume");
        let dout_wgt = matches.value_of("dir_wgt").map(|x| PathBuf::from(x));
        let fout_wgt = matches.value_of("file_wgt").map(|x| PathBuf::from(x));
        if dout_wgt.is_some() & fout_wgt.is_some() {
            panic!("Do not indicate both --dir_wgt and --file_wgt.");
        } else if dout_wgt.is_none() & fout_wgt.is_none() {
            panic!("Indicate either --dir_wgt or --file_wgt.");
        }


        genetics::run_score(
            &dout_score,
            &fin,
            fin_phe.as_deref(),
            phe_name.as_deref(),
            dout_wgt.as_deref(), // use enum?
            fout_wgt.as_deref(),
            //fin_cov.as_deref(),
            fin_sample.as_deref(),
            is_resume
        );
    }

    let end = start.elapsed();
    log::info!("It took {} seconds.", end.as_secs());
    log::info!("Done!!");

}
