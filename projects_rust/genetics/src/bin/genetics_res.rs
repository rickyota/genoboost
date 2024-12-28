//! General lib and application of general functions
//! - Polygenic score
//!

use clap::{ArgGroup, Args, Parser, Subcommand, ValueEnum};
use genetics::{DatasetFile, DoutScoreFile, GenotFile, WgtDoutOrFile};
use indoc::indoc;
use std::{path::PathBuf, time::Instant};

#[derive(Debug, Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    // global=true makes you able to `-- train --verbose`
    #[arg(long, global = true, help = "Number of threads")]
    threads: Option<usize>,
    #[arg(long, global = true, help = "Verbose")]
    verbose: bool,
    //#[arg(long, global = true, help = "Memory [GB]")]
    #[arg(long, global = true, help = "Memory [MB]")]
    memory: Option<usize>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(about = "score")]
    Score(ScoreArgs),
}

#[derive(Debug, Args)]
#[command(group(ArgGroup::new("wgt_dir_or_file").required(true).args(["dir_wgt", "file_wgt"])))]
#[command(group(ArgGroup::new("fill_missing").args(["missing_to_mode", "missing_to_mean"])))]
struct ScoreArgs {
    #[arg(long)]
    dir_score: String,
    #[arg(long, help = ".zst is allowed.")]
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
    file_freq: Option<String>,
    #[arg(long)]
    resume: bool,
    #[arg(
        long,
        help = "Concat parameters into one file. Wgt file name should be *_n-[#snv].wgt for --concat 'n'. "
    )]
    concat: Option<String>,
    #[arg(
        long,
        help = "Opposite of --concat. Use wgts not calculated in --concat."
    )]
    no_concat: Option<String>,
    #[arg(
        long,
        help = "Allow snvs or alleles in weight file not in genot. The score is ignored."
    )]
    allow_nonexist_snv: bool,
    #[arg(
        long,
        help = "When matching snvs in wgt and genot, use chromosome and posotion not variant id to match."
    )]
    use_snv_pos: bool,
    #[arg(long)]
    missing_to_mode: bool,
    #[arg(long)]
    missing_to_mean: bool,
    #[arg(long, help = "Use column score0-score2.")]
    nonadd: bool,
    #[arg(
        long,
        help = indoc!{"For backward compatibility."}
    )]
    fill_missing_in_dataset: bool,
    // TODO: maf file for missing values
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, ValueEnum)]
enum GenotFormatArg {
    Plink,
    Plink2,
    Plink2Vzs,
}

impl GenotFormatArg {
    pub fn to_genot_file(self, fin: PathBuf) -> GenotFile {
        match self {
            GenotFormatArg::Plink => GenotFile::Plink1(fin),
            GenotFormatArg::Plink2 => GenotFile::Plink2(fin),
            GenotFormatArg::Plink2Vzs => GenotFile::Plink2Vzs(fin),
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

    let mem = cli.memory.map(|x| x * 1024 * 1024);
    //let mem = cli.memory.map(|x| x * 1024 * 1024 * 1024);
    log::debug!("Memory : {:?} Byte", mem);

    match cli.command {
        Commands::Score(args) => {
            let dout_score = DoutScoreFile::new(PathBuf::from(args.dir_score));
            //let dout_score = PathBuf::from(args.dir_score);

            let fin = PathBuf::from(args.file_genot);
            let fin_genot = args.genot_format.to_genot_file(fin);
            //let genot_format = args.genot_format.to_naive();
            let fin_phe = args.file_phe.map(|x| PathBuf::from(x));
            //let phe_name = args.phe;
            let cov_name = args.cov;
            let fin_sample = args.file_sample.map(|x| PathBuf::from(x));
            let fin_freq = args.file_freq.map(|x| PathBuf::from(x));

            let mut dfile =
                DatasetFile::new(fin_genot, fin_phe, None, cov_name, None, fin_sample, None);
            dfile.update_file_freq(fin_freq);
            dfile.reads();
            let dfile = dfile;
            dfile.check_valid_fin();

            let dout_wgt = args.dir_wgt.map(|x| PathBuf::from(x));
            let fout_wgt = args.file_wgt.map(|x| PathBuf::from(x));
            let wgt_d_f = WgtDoutOrFile::new_path(dout_wgt, fout_wgt);

            //let is_resume = args.resume;

            let concat = args.concat;
            let no_concat = args.no_concat;

            if concat.is_some() && no_concat.is_some() {
                panic!("--concat and --no-concat cannot be used together.");
            }

            //let allow_nonexist_snv = args.allow_nonexist_snv;
            //let use_snv_pos = args.use_snv_pos;
            let is_nonadd = args.nonadd;

            genetics::run_score(
                &dout_score,
                &dfile,
                //&fin_genot,
                //fin_phe.as_deref(),
                //cov_name.as_deref(),
                &wgt_d_f,
                //dout_wgt.as_deref(), // use enum?
                //fout_wgt.as_deref(),
                //fin_sample.as_deref(),
                concat.as_deref(),
                no_concat.as_deref(),
                args.resume,
                args.fill_missing_in_dataset,
                args.allow_nonexist_snv,
                args.use_snv_pos,
                args.missing_to_mode,
                args.missing_to_mean,
                is_nonadd,
                mem,
            );
        }
    }

    let end = start.elapsed();
    log::info!("It took {} seconds.", end.as_secs());
    log::info!("Done!!");
}
