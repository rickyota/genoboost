//! General lib and application of general functions
//! - Polygenic score
//!
//!
//!

use clap::{Parser, ValueEnum};
//#[macro_use]
//extern crate clap;
//use clap::{AppSettings, Arg,  ArgMatches, SubCommand};
//use rayon;
use genetics::GenotFormat;
use std::{path::PathBuf, time::Instant};
//use clap::app_from_crate; //error

//#[derive(Debug, Parser)]
//struct Cli {
//    #[command(subcommand)]
//    command: Commands,
//
//    // global=true makes you able to `-- train --verbose`
//    #[arg(long, global = true, help = "Number of threads")]
//    threads: Option<usize>,
//    #[arg(long, global = true, help = "Verbose")]
//    verbose: bool,
//}

//#[derive(Debug, Subcommand)]
//enum Commands {
//    #[command(about = "score")]
//    Score(ScoreArgs),
//}

#[derive(Parser, Debug)]
#[command(author, version, about = "test pgenlib")]
struct Args {
    #[arg(long)]
    file_genot: String,
    #[arg(long, value_enum)]
    genot_format: GenotFormatArg,
    #[arg(long)]
    file_sample: Option<String>,
    #[arg(long, help = "Number of threads")]
    threads: Option<usize>,
    #[arg(long, help = "Verbose")]
    verbose: bool,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, ValueEnum)]
enum GenotFormatArg {
    Plink,
    Plink2,
    Plink2Vzs,
}

impl GenotFormatArg {
    pub fn to_naive(self) -> GenotFormat {
        match self {
            GenotFormatArg::Plink => GenotFormat::Plink1,
            GenotFormatArg::Plink2 => GenotFormat::Plink2,
            GenotFormatArg::Plink2Vzs => GenotFormat::Plink2Vzs,
        }
    }
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

    let start = Instant::now();

    let args = Args::parse();
    println!("args: {:?}", args);

    //let cli = Cli::parse();
    //println!("cli: {:?}", cli);

    let verbose = args.verbose;
    if verbose {
        std::env::set_var("RUST_LOG", "debug");
    } else {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::init();

    if let Some(n_threads) = args.threads {
        //n_threads=n_threads.min(num_cpus::get());
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .unwrap();
    };
    // otherwise, use default thread number
    log::debug!("num_thread set: {}", rayon::current_num_threads());

    //let dout_score = PathBuf::from(args.dir_score);
    let fin = PathBuf::from(args.file_genot);
    let _genot_format = args.genot_format.to_naive();
    //let fin_phe = args.file_phe.map(|x| PathBuf::from(x));
    //let phe_name = args.phe;
    //let cov_name = args.cov;
    let fin_sample = args.file_sample.map(|x| PathBuf::from(x));

    //let dout_wgt = args.dir_wgt.map(|x| PathBuf::from(x));
    //let fout_wgt = args.file_wgt.map(|x| PathBuf::from(x));

    //let is_resume = args.resume;

    //let concat = args.concat;

    println!("start test_pgenlib");
    //let g=genetics::test_pgenlib(&fin, genot_format, fin_sample.as_deref());

    let g_plink2 = genetics::test_pgenlib(&fin, GenotFormat::Plink2Vzs, fin_sample.as_deref());
    let g_plink1 = genetics::test_pgenlib(&fin, GenotFormat::Plink1, fin_sample.as_deref());

    //assert_eq!(g_plink2, g_plink1);
    if g_plink2 == g_plink1 {
        println!("OK: plink1==plink2");
    } else {
        println!("**WARNING**: plink1!=plink2");
    }

    let end = start.elapsed();
    log::info!("It took {} seconds.", end.as_secs());
    log::info!("Done!!");
}
