//! Application of **Genoboost**.
//! 
//! same input as main.rs in gwasmethods should work
//!
//! Logitnomissing
//! When the denominator of s2 is 0.0, (no eps and no samples for minor homozygotes), s2 is set the same as s1 (=dominant model).
//! When the denominator of s1, s0 is 0.0, (no samples for major homozygotes or heterozygotes), s0, s1 is set to 0.0
//!
//!
// TODO: ensure the same para when resuming
// TODO: (optional) write down extract snvs from top
// TODO: how to get memory?
// samples indicate both fid and iid
// TODO: trimming sample weights: only use samples with large weights on choosing SNVs: Friedman, J., Hastie, T. and Tibshirani, R. (2000) ‘Additive logistic regression: a statistical view of boosting (With discussion and a rejoinder by the authors)’, Annals of statistics, 28(2), pp. 337–407. doi:10.1214/aos/1016218223.
// split main for genoboost_pruning, ssgenoboost
// TODO: use ref/alt not A1/A2 and fix lims2
// TODO: format of use_snvs, use_sample should be the same as plink --extract

use clap::{Args, Parser, Subcommand};
//use crate::boosting::{BoostMethod, BoostParam, IterationNumber};
use boosting::{self, BoostMethod, BoostParam};
//use rayon;
use genetics::GenotFormat;
use std::{path::PathBuf, time::Instant};

#[derive(Debug, Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    // globa=true makes you able to `-- trian --verbose`
    #[arg(long, global = true, help = "Number of threads")]
    num_threads: Option<usize>,
    #[arg(long, global = true, help = "Verbose")]
    verbose: bool,
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(about = "train")]
    Train(TrainArgs),
    #[command(about = "score")]
    Score(ScoreArgs),
}

#[derive(Debug, Args)]
//#[command(group(ArgGroup::new("how_to_input").required(true).args(["input", "input_from_stdin"])))]
struct TrainArgs {
    #[arg(long)]
    dir: String,
    #[arg(long)]
    file_plink: String,
    //#[arg(long,value_enum,default_value_t = BoostType::genoboost)]
    //boost_type: BoostType,
    #[arg(long, default_value_t = String::from("genoboost"))]
    boost_type: String,
    #[arg(long)]
    file_sample: Option<String>,
    #[arg(long)]
    file_sample_val: Option<String>,
    #[arg(long)]
    file_phe: Option<String>,
    #[arg(long)]
    phe: Option<String>,
    #[arg(long)]
    file_cov: Option<String>,
    // TODO: require file_cov
    #[arg(long)]
    cov: Option<String>,
    #[arg(long)]
    file_snv: Option<String>,
    #[arg(long)]
    iter_snv: Option<usize>,
    #[arg(long)]
    iter: Option<usize>,
    #[arg(long)]
    learning_rates: Option<Vec<f64>>,
    #[arg(long)]
    batch: Option<String>,
    #[arg(long)]
    clip_sample_weight: Option<String>,
    #[arg(long)]
    clip_sample_wls_weight: Option<String>,
    #[arg(long)]
    eps: Option<String>,
    #[arg(long)]
    eff_eps: Option<String>,
    #[arg(long)]
    use_adjloss: bool,
    #[arg(long)]
    use_const_for_loss: bool,
    #[arg(long)]
    resume: bool,
    #[arg(long)]
    write_loss: bool,
}

#[derive(Debug, Args)]
struct ScoreArgs {
    #[arg(long)]
    dir_score: String,
    #[arg(long)]
    file_plink: String,
    #[arg(long)]
    dir_wgt: Option<String>,
    #[arg(long)]
    file_wgt: Option<String>,
    #[arg(long)]
    file_sample: Option<String>,
    #[arg(long)]
    file_phe: Option<String>,
    #[arg(long)]
    phe: Option<String>,
    #[arg(long)]
    file_cov: Option<String>,
    #[arg(long)]
    iters: Vec<usize>,
    #[arg(long)]
    learning_rates: Option<Vec<f64>>,
    #[arg(long)]
    use_iter: bool,
}

//#[derive(Debug, Clone, ValueEnum)]
//enum BoostType {
//    genoboost,
//    ada,
//    constada,
//    freemodelmissing,
//    logit,
//    logitnomissing,
//    logitadd,
//}

// TODO: do not allow string "None"

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

    if let Some(n_threads) = cli.num_threads {
        //n_threads=n_threads.min(num_cpus::get());
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .unwrap();
    };
    // otherwise, use default thread number
    log::debug!("num_thread set: {}", rayon::current_num_threads());

    match cli.command {
        Commands::Train(args) => {
            let dout = PathBuf::from(args.dir);
            let fin = PathBuf::from(args.file_plink);
            let fin_phe = args.file_phe.map(|x| PathBuf::from(x));
            // TODO: does not allow "None"
            let phe_name = match args.phe {
                None => None,
                // &*String -> str
                Some(y) => match &*y {
                    "None" => None,
                    z => Some(z.to_string()),
                },
            };
            let fin_sample = args.file_sample.map(|x| PathBuf::from(x));
            let fin_cov = args.file_cov.map(|x| PathBuf::from(x));
            let fin_sample_val = args.file_sample_val.map(|x| PathBuf::from(x));
            let is_monitor = fin_sample_val.is_some();
            let is_resume = args.resume;
            let fin_snv = args.file_snv.map(|x| PathBuf::from(x));

            let learning_rates: Vec<Option<f64>> = match args.learning_rates {
                None => vec![None],
                Some(lrs) => lrs.iter().map(|x| Some(*x)).collect(),
            };
            let cov_way = match args.cov {
                None => None,
                Some(y) => match &*y {
                    "None" => None,
                    z => Some(z.to_string()),
                },
            };
            let batch_way = match args.batch {
                None => None,
                Some(y) => match &*y {
                    "None" => None,
                    z => Some(z.to_string()),
                },
            };
            let eps_way = match args.eps {
                None => None,
                Some(y) => match &*y {
                    "None" | "none" => None,
                    z => Some(z.to_string()),
                },
            };
            let eff_eps = match args.eff_eps {
                None => None,
                Some(y) => match &*y {
                    "None" => None,
                    z => Some(z.to_string()),
                },
            };
            let clip_sample_weight = match args.clip_sample_weight {
                None => None,
                Some(y) => match &*y {
                    "None" => None,
                    z => Some(z.to_string()),
                },
            };
            let clip_sample_wls_weight = match args.clip_sample_wls_weight {
                None => None,
                Some(y) => match &*y {
                    "None" => None,
                    z => Some(z.to_string()),
                },
            };

            // set learning rate later
            let boost_param = BoostParam::default()
                .set_boost_type(&*args.boost_type)
                .set_loss_func("logistic")
                .set_sample_weight_clip(clip_sample_weight.as_deref())
                .set_sample_weight_wls_clip(clip_sample_wls_weight.as_deref())
                .set_eps(eps_way.as_deref())
                .set_cov_way(cov_way.as_deref())
                .set_batch_way(batch_way.as_deref())
                .set_eff_eps(eff_eps.as_deref());

            // TODO: better way in clap
            if args.iter.is_some() & args.iter_snv.is_some() {
                panic!("Do not indicate both --iter and --iter_snv.");
            } else if args.iter.is_none() & args.iter_snv.is_none() {
                panic!("Indicate either --iter or --iter_snv.");
            }

            // TODO: cleaner
            let boost_param = if args.iter.is_some() {
                boost_param.set_iteration(args.iter.unwrap())
            } else {
                boost_param.set_iteration_snv(args.iter_snv.unwrap())
            };
            boost_param.check();
            log::debug!("boost_param {:?}", boost_param);

            let boost_method = BoostMethod::Classic;

            log::info!("dout {:?}", dout);
            log::info!("file_plink {:?}", fin);
            log::info!("file_snv {:?}", fin_snv);
            log::info!("file_cov {:?}", fin_cov);
            log::info!("file_sample {:?}", fin_sample);
            log::info!("boost_param {:?}", boost_param);

            let use_adjloss = args.use_adjloss;
            let use_const_for_loss = args.use_const_for_loss;
            let is_write_loss = args.write_loss;

            // Genoboost
            crate::boosting::run_boosting(
                &dout,
                &fin,
                GenotFormat::Plink,
                fin_phe.as_deref(),
                phe_name.as_deref(),
                None,
                boost_method,
                boost_param,
                fin_snv.as_deref(),
                fin_sample.as_deref(),
                fin_cov.as_deref(),
                fin_sample_val.as_deref(),
                use_adjloss,
                use_const_for_loss,
                is_resume,
                is_write_loss,
                None, //prune_snv,
                &learning_rates,
                is_monitor,
            );
        }
        Commands::Score(args) => {
            let dout_score = PathBuf::from(args.dir_score);
            let fin = PathBuf::from(args.file_plink);
            let fin_phe = args.file_phe.map(|x| PathBuf::from(x));
            let phe_name = match args.phe {
                None => None,
                // &*String -> str
                Some(y) => match &*y {
                    "None" => None,
                    z => Some(z.to_string()),
                },
            };
            let fin_sample = args.file_sample.map(|x| PathBuf::from(x));
            let fin_cov = args.file_cov.map(|x| PathBuf::from(x));

            // TODO: better way in clap
            if args.dir_wgt.is_some() & args.file_wgt.is_some() {
                panic!("Do not indicate both --iter and --iter_snv.");
            } else if args.dir_wgt.is_none() & args.file_wgt.is_none() {
                panic!("Indicate either --iter or --iter_snv.");
            }

            let dout_wgt = args.dir_wgt.map(|x| PathBuf::from(x));
            let fout_wgt = args.file_wgt.map(|x| PathBuf::from(x));
            let mut iterations = args.iters;
            iterations.sort();
            iterations.dedup();
            log::info!("iters {:?}", iterations);
            let learning_rates: Vec<Option<f64>> = match args.learning_rates {
                None => vec![None],
                Some(lrs) => lrs.iter().map(|x| Some(*x)).collect(),
            };
            let use_iter = args.use_iter;

            crate::boosting::run_boosting_score(
                &dout_score,
                &fin,
                GenotFormat::Plink,
                fin_phe.as_deref(),
                phe_name.as_deref(),
                None,
                false,
                Some(&iterations),
                dout_wgt.as_deref(), // use enum?
                fout_wgt.as_deref(),
                fin_cov.as_deref(),
                fin_sample.as_deref(),
                //boost_param,
                &learning_rates,
                use_iter,
            );
        }
    }

    let end = start.elapsed();
    log::info!("It took {} seconds.", end.as_secs());
    log::info!("Done!!");
}
