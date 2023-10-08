//! General lib and application of general functions
//! - Polygenic score
//!

use clap::{ArgGroup, Args, Parser, Subcommand, ValueEnum};
//#[macro_use]
//extern crate clap;
//use clap::{AppSettings, Arg,  ArgMatches, SubCommand};
//use rayon;
use genetics::GenotFormat;
use std::{path::PathBuf, time::Instant};
//use clap::app_from_crate; //error

#[derive(Debug, Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    // global=true makes you able to `-- train --verbose`
    #[arg(long, global = true, help = "Number of threads")]
    threads: Option<usize>,
    #[arg(long, global = true, help = "Verbose")]
    verbose: bool,
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(about = "score")]
    Score(ScoreArgs),
}

#[derive(Debug, Args)]
#[command(group(ArgGroup::new("wgt_dir_or_file").required(true).args(["dir_wgt", "file_wgt"])))]
struct ScoreArgs {
    #[arg(long)]
    dir_score: String,
    #[arg(long)]
    file_genot: String,
    #[arg(long, value_enum)]
    genot_format: GenotFormatArg,
    #[arg(long)]
    dir_wgt: Option<String>,
    #[arg(long)]
    file_wgt: Option<String>,
    #[arg(long)]
    file_sample: Option<String>,
    #[arg(long)]
    file_phe: Option<String>,
    //#[arg(long)]
    //phe: Option<String>,
    #[arg(long)]
    cov: Option<String>,
    #[arg(long)]
    resume: bool,
    #[arg(
        long,
        help = "Concat parameters into one file. Only use wgt file name is *_n-[#snv].wgt for --concat 'n'. "
    )]
    concat: Option<String>,
    #[arg(
        long,
        help = "Opposite of --concat. Use wgts not calculated in --concat."
    )]
    no_concat: Option<String>,
    #[arg(
        long,
        help = "Allow snvs or alleles not in genot. The score is ignored."
    )]
    allow_nonexist_snv: bool,
    #[arg(
        long,
        help = "When matching snvs in wgt and genot, use chromosome and posotion not variant id to match."
    )]
    use_snv_pos: bool,
    #[arg(long, help = "Use column score0-score2.")]
    nonadd: bool,
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

    let cli = Cli::parse();
    println!("cli: {:?}", cli);

    let verbose = cli.verbose;
    if verbose {
        std::env::set_var("RUST_LOG", "debug");
    } else {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::init();

    if let Some(n_threads) = cli.threads {
        //n_threads=n_threads.min(num_cpus::get());
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .unwrap();
    };
    // otherwise, use default thread number
    log::debug!("num_thread set: {}", rayon::current_num_threads());

    match cli.command {
        Commands::Score(args) => {
            let dout_score = PathBuf::from(args.dir_score);
            let fin = PathBuf::from(args.file_genot);
            let genot_format = args.genot_format.to_naive();
            let fin_phe = args.file_phe.map(|x| PathBuf::from(x));
            //let phe_name = args.phe;
            let cov_name = args.cov;
            let fin_sample = args.file_sample.map(|x| PathBuf::from(x));

            let dout_wgt = args.dir_wgt.map(|x| PathBuf::from(x));
            let fout_wgt = args.file_wgt.map(|x| PathBuf::from(x));

            let is_resume = args.resume;

            let concat = args.concat;
            let no_concat = args.no_concat;

            if concat.is_some() & no_concat.is_some() {
                panic!("--concat and --no-concat cannot be used together.");
            }

            let allow_nonexist_snv = args.allow_nonexist_snv;
            let use_snv_pos = args.use_snv_pos;
            let is_nonadd = args.nonadd;

            genetics::run_score(
                &dout_score,
                &fin,
                genot_format,
                fin_phe.as_deref(),
                //phe_name.as_deref(),
                cov_name.as_deref(),
                dout_wgt.as_deref(), // use enum?
                fout_wgt.as_deref(),
                //fin_cov.as_deref(),
                fin_sample.as_deref(),
                concat.as_deref(),
                no_concat.as_deref(),
                is_resume,
                allow_nonexist_snv,
                use_snv_pos,
                is_nonadd,
            );
        }
    }

    let end = start.elapsed();
    log::info!("It took {} seconds.", end.as_secs());
    log::info!("Done!!");
}
