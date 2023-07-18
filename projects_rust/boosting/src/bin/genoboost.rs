//! Application of **Genoboost** for research.
//! Input plink file to run Genoboost.
//!
//! changed inputs;
//!     - GenotFormatArg
//!     - remove fin_cov
//!
//! do not allow string "None"
//!
//!
//! Input is one of following
//! 1. plink2 + fin_phe
//! 2. plink2 (cov and phe in .psam) : fin_phe is option
//! 3. plink1 (phe) + fin_phe (cov) : phe_name is option
//! 4. plink1 + fin_phe (cov and phe)
//!
//!
//! TODO: assume FID==IID for plink1 or plink2
//!
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

use clap::{ArgGroup, Args, Parser, Subcommand, ValueEnum};
//use crate::boosting::{BoostMethod, BoostParam, IterationNumber};
use boosting::{self, BoostMethod, BoostParamsTypes};
use genetics::GenotFormat;
//use rayon;
use std::{path::PathBuf, time::Instant};

#[derive(Debug, Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    // global=true makes you able to `-- trian --verbose`
    #[arg(long, global = true, help = "Number of threads")]
    threads: Option<usize>,
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

//#[command(group(ArgGroup::new("iter_or_snv").required(true).args(["iter_snv", "iter"])))]
#[derive(Debug, Args)]
#[command(group(ArgGroup::new("iter_or_snv").args(["iter_snv", "iter"])))]
struct TrainArgs {
    #[arg(long)]
    dir: String,
    #[arg(long)]
    file_genot: String,
    #[arg(long, value_enum, default_value_t = GenotFormatArg::Plink1)]
    genot_format: GenotFormatArg,
    //#[arg(long,value_enum,default_value_t = BoostType::genoboost)]
    //boost_type: BoostType,
    //#[arg(long, default_value_t = String::from("genoboost"))]
    // only use defeault
    //#[arg(long, default_value_t = String::from("nonadd"))]
    //boost_type: String,
    #[arg(long)]
    file_sample: Option<String>,
    #[arg(long)]
    file_sample_val: Option<String>,
    // option for covs and phes are in .psam
    #[arg(long)]
    file_phe: Option<String>,
    // option for phe in plink1
    #[arg(long)]
    phe: Option<String>,
    // parse later
    #[arg(long)]
    cov: String,
    //#[arg(long)]
    //file_cov: Option<String>,
    #[arg(long)]
    covway: Option<String>,
    #[arg(long)]
    file_snv: Option<String>,
    #[arg(long)]
    iter_snv: Option<usize>,
    #[arg(long)]
    iter: Option<usize>,
    #[arg(long, value_parser, num_args = 1.., value_delimiter = ' ', default_value = "0.5 0.2 0.1 0.05")]
    learning_rates: Vec<f64>,
    //learning_rates: Option<Vec<f64>>,
    //#[arg(long)]
    //batch: Option<String>,
    //#[arg(long)]
    //clip_sample_weight: Option<String>,
    //#[arg(long)]
    //clip_sample_wls_weight: Option<String>,
    //#[arg(long)]
    //eps: Option<String>,
    //#[arg(long)]
    //effeps: Option<String>,
    #[arg(long)]
    resume: bool,
    #[arg(long)]
    write_loss: bool,
    // --integrate-only
    //#[arg(long, default_value_t = true)]
    //integrate: bool,
    // TODO: boost_types for --integrate
    /// Number of cross validation. Ignored with --file-sample-val
    #[arg(long, default_value_t = 1)]
    //cross_validation: usize,
    cross_validation: usize,
    /// For sample split in cross validation
    #[arg(long)]
    seed: Option<u64>,
    #[arg(long)]
    train_only: bool,
}

#[derive(Debug, Args)]
#[command(group(ArgGroup::new("wgt").required(true).args(["dir_wgt", "file_wgt"])))]
struct ScoreArgs {
    #[arg(long)]
    dir_score: String,
    #[arg(long)]
    file_genot: String,
    #[arg(long, value_enum, default_value_t = GenotFormatArg::Plink1)]
    genot_format: GenotFormatArg,
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
    // TODO: remove --cov and read from wgt?
    #[arg(long)]
    cov: String,
    //#[arg(long)]
    //file_cov: Option<String>,
    // if indicated, do not use para_best and calc score of all paras
    #[arg(long)]
    iters: Option<Vec<usize>>,
    #[arg(long)]
    learning_rates: Option<Vec<f64>>,
    #[arg(long)]
    use_iter: bool,
    /// calculate scores for all cross validation weights
    #[arg(long, default_value_t = 1)]
    cross_validation: usize,
    #[arg(long)]
    train_only: bool,
}

// concat with boosting_rest
#[derive(Copy, Clone, PartialEq, Eq, Debug, ValueEnum)]
enum GenotFormatArg {
    Plink1,
    Plink2,
    Plink2Vzs,
}

impl GenotFormatArg {
    pub fn to_naive(self) -> GenotFormat {
        match self {
            GenotFormatArg::Plink1 => GenotFormat::Plink1,
            GenotFormatArg::Plink2 => GenotFormat::Plink2,
            GenotFormatArg::Plink2Vzs => GenotFormat::Plink2Vzs,
        }
    }
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
        Commands::Train(args) => {
            let dout = PathBuf::from(args.dir);
            //let fin = PathBuf::from(args.file_plink);
            let fin = PathBuf::from(args.file_genot);
            let genot_format = args.genot_format.to_naive();
            let fin_phe = args.file_phe.map(|x| PathBuf::from(x));
            let phe_name = args.phe;
            let cov_name = args.cov;
            let fin_sample = args.file_sample.map(|x| PathBuf::from(x));
            //let fin_cov = args.file_cov.map(|x| PathBuf::from(x));
            let fin_sample_val = args.file_sample_val.map(|x| PathBuf::from(x));
            let is_monitor = fin_sample_val.is_some();
            let is_resume = args.resume;
            let fin_snv = args.file_snv.map(|x| PathBuf::from(x));

            let learning_rates: Vec<f64> = args.learning_rates;

            //let learning_rates: Vec<f64> = match args.learning_rates {
            //    None => vec![1.0f64],
            //    Some(lrs) => lrs.iter().map(|x| *x).collect(),
            //};
            let cov_way = Some("first".to_string());
            //let cov_way = args.covway;
            let batch_way = Some("fix-50-10-50".to_string());
            //let batch_way = args.batch;
            let eps_way: Option<String> = None;
            //let eps_way = args.eps;
            let eff_eps = Some("lims2gmodeloverkeepsignprop-100-4-4-0.8".to_string());
            //let eff_eps = args.effeps;
            let clip_sample_weight: Option<String> = None;
            //let clip_sample_weight = args.clip_sample_weight;
            let clip_sample_wls_weight: Option<String> = None;
            //let clip_sample_wls_weight = args.clip_sample_wls_weight;

            // TODO: raise error if duplicate nonadd in .set_boost_types()
            let boost_params = BoostParamsTypes::default()
                .set_boost_types(vec!["nonadd".to_string(), "add".to_string()])
                .set_loss_func("logistic")
                .set_learning_rates(learning_rates)
                .set_sample_weight_clip(clip_sample_weight.as_deref())
                .set_sample_weight_wls_clip(clip_sample_wls_weight.as_deref())
                .set_eps(eps_way.as_deref())
                .set_cov_way(cov_way.as_deref())
                .set_batch_way(batch_way.as_deref())
                .set_eff_eps(eff_eps.as_deref());

            //let boost_params =if args.integrate{
            //BoostParams::default()
            //    .set_boost_types(vec!["nonadd","add"])
            //    .set_loss_func("logistic")
            //    .set_learning_rates(learning_rates)
            //    .set_sample_weight_clip(clip_sample_weight.as_deref())
            //    .set_sample_weight_wls_clip(clip_sample_wls_weight.as_deref())
            //    .set_eps(eps_way.as_deref())
            //    .set_cov_way(cov_way.as_deref())
            //    .set_batch_way(batch_way.as_deref())
            //    .set_eff_eps(eff_eps.as_deref())
            //}else{
            //BoostParams::default()
            //    .set_boost_type(&*args.boost_type)
            //    .set_loss_func("logistic")
            //    .set_learning_rates(learning_rates)
            //    .set_sample_weight_clip(clip_sample_weight.as_deref())
            //    .set_sample_weight_wls_clip(clip_sample_wls_weight.as_deref())
            //    .set_eps(eps_way.as_deref())
            //    .set_cov_way(cov_way.as_deref())
            //    .set_batch_way(batch_way.as_deref())
            //    .set_eff_eps(eff_eps.as_deref())
            //};

            let boost_params = if args.iter.is_some() {
                boost_params.set_iteration(args.iter.unwrap())
            } else if args.iter_snv.is_some() {
                boost_params.set_iteration_snv(args.iter_snv.unwrap())
            } else {
                // else: integrate
                boost_params
            };

            //boost_params.check();
            log::debug!("boost_param {:?}", boost_params);

            let boost_method = BoostMethod::Classic;

            log::info!("dout {:?}", dout);
            log::info!("file_plink {:?}", fin);
            log::info!("file_snv {:?}", fin_snv);
            //log::info!("file_cov {:?}", fin_cov);
            log::info!("file_sample {:?}", fin_sample);
            log::info!("boost_params {:?}", boost_params);

            let use_adjloss = true;
            //let use_adjloss = args.use_adjloss;
            let use_const_for_loss = false;
            let is_write_loss = args.write_loss;

            let cross_vali: Option<usize> = if args.train_only {
                None
            } else {
                Some(args.cross_validation)
            };
            let seed = args.seed;

            // Genoboost
            crate::boosting::run_boosting_integrate_cv(
                &dout,
                &fin,
                genot_format,
                fin_phe.as_deref(),
                phe_name.as_deref(),
                &cov_name,
                boost_method,
                &boost_params,
                fin_snv.as_deref(),
                fin_sample.as_deref(),
                //fin_cov.as_deref(),
                fin_sample_val.as_deref(),
                use_adjloss,
                use_const_for_loss,
                is_resume,
                is_write_loss,
                None, //prune_snv,
                //&learning_rates,
                is_monitor,
                cross_vali,
                seed,
            );
        }
        Commands::Score(args) => {
            let dout_score = PathBuf::from(args.dir_score);
            let fin = PathBuf::from(args.file_genot);
            let genot_format = args.genot_format.to_naive();
            let fin_phe = args.file_phe.map(|x| PathBuf::from(x));
            let phe_name = args.phe;
            let cov_name = args.cov;
            let fin_sample = args.file_sample.map(|x| PathBuf::from(x));
            //let fin_cov = args.file_cov.map(|x| PathBuf::from(x));

            let dout_wgt = args.dir_wgt.map(|x| PathBuf::from(x));
            let fout_wgt = args.file_wgt.map(|x| PathBuf::from(x));
            let is_every_para = args.iters.is_some();
            let iterations = if is_every_para {
                let mut iterations = args.iters.unwrap();
                iterations.sort();
                iterations.dedup();
                log::info!("iters {:?}", iterations);
                Some(iterations)
            } else {
                None
            };
            let learning_rates: Vec<f64> = match args.learning_rates {
                None => vec![1.0f64],
                Some(lrs) => lrs.iter().map(|x| *x).collect(),
            };
            //let learning_rates: Vec<Option<f64>> = match args.learning_rates {
            //    None => vec![None],
            //    Some(lrs) => lrs.iter().map(|x| Some(*x)).collect(),
            //};
            let use_iter = args.use_iter;
            let cross_vali: Option<usize> = if args.train_only {
                None
            } else {
                Some(args.cross_validation)
            };

            crate::boosting::run_boosting_score_cv(
                &dout_score,
                &fin,
                genot_format,
                fin_phe.as_deref(),
                phe_name.as_deref(),
                Some(&cov_name),
                is_every_para,
                iterations.as_deref(),
                dout_wgt.as_deref(), // use enum?
                fout_wgt.as_deref(),
                //fin_cov.as_deref(),
                fin_sample.as_deref(),
                //boost_param,
                &learning_rates,
                use_iter,
                cross_vali,
            );
        }
    }

    let end = start.elapsed();
    log::info!("It took {} seconds.", end.as_secs());
    log::info!("Done!!");
}
