//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
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

#[macro_use]
extern crate clap;
//use crate::boosting::{BoostMethod, BoostParam, IterationNumber};
use boosting::{self, BoostMethod, BoostParam};
use clap::{AppSettings, Arg, ArgGroup, ArgMatches, SubCommand};
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

fn check_f64(v: String) -> Result<(), String> {
    match v.parse::<f64>() {
        Ok(_) => Ok(()),
        Err(_) => Err("The value should be able to be parsed into float value.".to_owned()),
    }
}

fn parse_arg_f64_none(v: &str) -> Option<f64> {
    if v == "None" {
        return None;
    }
    Some(v.parse::<f64>().unwrap())
}

// need this for clap
fn check_f64_none(v: String) -> Result<(), String> {
    if v == "None" {
        return Ok(());
    }
    match v.parse::<f64>() {
        Ok(_) => Ok(()),
        Err(_) => Err("The value should be able to be parsed into float value.".to_owned()),
    }
}

// TODO: use clap v4
// TODO: do not allow string "None"
fn get_matches() -> ArgMatches<'static> {
    let matches = app_from_crate!()
        .setting(AppSettings::DeriveDisplayOrder)
        .setting(AppSettings::SubcommandRequiredElseHelp) //require subcommand or show help
        .arg(Arg::from_usage("--num_threads [NUM] 'Number of threads.'").validator(check_usize))
        .arg(Arg::from_usage("--verbose 'Verbose.'"))
        .subcommand( SubCommand::with_name("train").about("Train model.")
            //.subcommand( SubCommand::with_name("ss").about("Only use summary statistics."))
            .arg(Arg::from_usage("--use_adjloss 'Use adjusted loss.'"))
            .arg(Arg::from_usage("--use_const_for_loss 'Use const value for adjusted loss.'"))
            .arg(Arg::from_usage("--resume 'Resume training.'"))
            .arg(Arg::from_usage("--cov [COV] 'When to use cov. first/every(integer)/none'").requires("file_cov"))
            .arg(Arg::from_usage("--batch [BATCH] 'When to use cov. every(integer)/none'").requires("file_cov"))
            .arg(Arg::from_usage("--write_loss 'Write down loss.'"))
            .arg(Arg::from_usage("--ss 'Only use summary statistics.'").requires("fin_loss"))
            .arg(Arg::from_usage("--pruning 'Use pruning method.'"))
            .group(ArgGroup::with_name("method").args(&["ss","pruning"]))
            // <FILE> : required
            // [FILE] : optional
            .arg(Arg::from_usage("--dir <FILE> 'Directory to output.'"))
            .arg(Arg::from_usage("--file_plink <FILE> 'Prefix of a plink file. If those files are separated into chromosomes, set '%' where chromosome number is. If you use subcommand 'ss', indicate reference panel.'"))
            .arg(Arg::from_usage("--iter [NUM] 'Number of iterations.'"))
            //.arg(Arg::from_usage("--iter [NUM] 'Number of iterations.'").default_value("1000").validator(check_usize))
            .arg(Arg::from_usage("--iter_snv [NUM] 'Number of SNVs for iterations.'"))
            //.arg(Arg::from_usage("--iter_snv [NUM] 'Number of SNVs for iterations.'").default_value("1000").validator(check_usize))
            .arg(
                Arg::from_usage("--boost_type [BOOSTTYPE] 'Boosting method to use.'")
                    .possible_values(&["genoboost","ada","constada","freemodelmissing","logit","logitnomissing","logitadd"])
                    .default_value("genoboost"),
            )
            .arg(Arg::from_usage("--file_sample [FILE] 'File of samples to use.'"))
            .arg(Arg::from_usage("--file_sample_val [FILE] 'File of samples to use.'"))
            .arg(Arg::from_usage("--file_phe [FILE] 'File of phenotype to use.'"))
            // TODO
            //.arg(Arg::from_usage("--file_plink_val [FILE] 'Prefix of a plink file for validation. If those files are separated into chromosomes, set '%' where chromosome number is. If you use subcommand 'ss', indicate reference panel.'"))
            //.arg(Arg::from_usage("--file_phe_val [FILE] 'File of phenotype to use.'"))
            .arg(Arg::from_usage("--phe [VAL] 'Phonotype. must use with --file_phe.'").default_value("None")) // TODO: validator
            .arg(Arg::from_usage("--file_snv [FILE] 'File of SNVs to use.'"))
            .arg(Arg::from_usage("--file_cov [FILE] 'File of covariates to use.'"))
            .arg(Arg::from_usage("--file_loss [FILE] 'File of loss function, which required only for --ss'"))
            .arg(Arg::from_usage("--learning_rates [VALS] 'Learning rates.'").default_value("None").multiple(true))
            //.arg(Arg::from_usage("--learning_rate [VAL] 'Learning rate.'").default_value("None").validator(check_f64_none))
            .arg(Arg::from_usage("--clip_sample_weight [VAL] 'Proportion and method of clipped samples.'").default_value("None")) // TODO: validator
            .arg(Arg::from_usage("--clip_sample_wls_weight [VAL] 'Proportion and method of clipped samples.'").default_value("None")) // TODO: validator
            .arg(Arg::from_usage("--eps [VAL] 'Eps.'").default_value("medcase2allcell")) // TODO: validator
            .arg(Arg::from_usage("--eff_eps [VAL] 'Eff Eps.'").default_value("None")) // TODO: validator
            .arg(Arg::from_usage("--dom_and_rec 'Choose both dominant model and recessive model at one iteration.'"))
            .arg(Arg::from_usage("--prune_snv [VAL] 'Proportion of SNVs to be pruned from the candidate SNVs.'").default_value("None").validator(check_f64_none))
            .arg(Arg::from_usage("--clump_r2 [VAL] 'Threshold of clump correlation.'").default_value("0.1").validator(check_f64))
            // use .required_if() for r2 for pruning-genoboost
        )
        .subcommand( SubCommand::with_name("score").about("Calculate score.")
            .arg(Arg::from_usage("--dir_score <FILE> 'Directory of output file.'"))
            //.arg(Arg::from_usage("--dout_score <FILE> 'Prefix of output file.'"))
            //.arg(Arg::from_usage("--fout <FILE> 'Prefix of output file.'"))
            .arg(Arg::from_usage("--file_plink <FILE> 'Prefix of a plink file. If those files are separated into chromosomes, set '%' where chromosome number is.'"))
            .arg(Arg::from_usage("--file_phe [FILE] 'File of phenotype to use.'"))
            .arg(Arg::from_usage("--phe [VAL] 'Phonotype. must use with --file_phe.'").default_value("None")) // TODO: validator
            .arg(Arg::from_usage("--file_sample [FILE] 'File of samples to use.'"))
            .arg(Arg::from_usage("--file_cov [FILE] 'File of covariates to use.'"))
            .arg(Arg::from_usage("--dir_wgt [FILE] 'Directory of a weight file. Must use either --dir_wgt or --file_wgt.'"))
            .arg(Arg::from_usage("--file_wgt [FILE] 'Prefix of a weight file. Must use either --dir_wgt or --file_wgt.'"))
            //.arg(Arg::from_usage("--fin_wgt <FILE> 'Prefix of a weight file.'"))
            // FIXME: default=None
            .arg(Arg::from_usage("--iters [NUMS] 'Numbers of iterations and number of SNVs.'").multiple(true))
            //.arg(Arg::from_usage("--iters [NUMS] 'Numbers of iterations.'").default_value("10").min_values(1))
            //.arg( Arg::from_usage("--iters [NUMS] 'Numbers of iterations.'").default_values(&["1", "10"]).min_values(1))
            .arg(Arg::from_usage("--learning_rates [VALS] 'Learning rates. Necessary when use --dir_wgt'").default_value("None").multiple(true))
            .arg(Arg::from_usage("--use_iter 'Output score on iterations. Otherwise, only output number of SNVs.'"))
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

    if let Some(ref matches) = matches.subcommand_matches("train") {
        let dout = PathBuf::from(matches.value_of("dir").unwrap());
        //let fout = matches.value_of("fout").unwrap();
        let fin = PathBuf::from(matches.value_of("file_plink").unwrap());
        let fin_phe = matches.value_of("file_phe").map(|x| PathBuf::from(x));
        let phe_name = match matches.value_of("phe").unwrap() {
            "None" => None,
            z => Some(z),
        };
        let fin_sample = matches.value_of("file_sample").map(|x| PathBuf::from(x));
        let fin_cov = matches.value_of("file_cov").map(|x| PathBuf::from(x));
        //let fin_cov = matches.value_of("fin_cov");

        let fin_sample_val = matches
            .value_of("file_sample_val")
            .map(|x| PathBuf::from(x));
        let is_monitor = fin_sample_val.is_some();

        //let boost_type = matches.value_of("boost_type").unwrap();
        let learning_rates: Vec<Option<f64>> = matches
            .values_of("learning_rates")
            .unwrap()
            .map(|s| parse_arg_f64_none(s))
            .collect();
        let prune_snv = parse_arg_f64_none(matches.value_of("prune_snv").unwrap());

        let is_write_loss = matches.is_present("write_loss");
        //let is_write_loss = true;

        // make this Option<> here? or in Boostingparam?
        // should be here, or there should be option to make Option<> as default
        let clip_sample_weight = match matches.value_of("clip_sample_weight").unwrap() {
            "None" => None,
            z => Some(z),
        };

        let clip_sample_wls_weight = match matches.value_of("clip_sample_wls_weight").unwrap() {
            "None" => None,
            z => Some(z),
        };

        let eps_way = match matches.value_of("eps").unwrap() {
            "None" | "none" => None,
            z => Some(z),
        };
        //let eps_way = matches.value_of("eps").unwrap();

        let is_dom_rec = matches.is_present("dom_and_rec");

        // FIXME: add assert for no --cov and --file_cov exist
        let cov_way = match matches.is_present("cov") {
            false => None,
            true => match matches.value_of("cov").unwrap() {
                "None" => None,
                z => Some(z),
            },
        };
        let batch_way = match matches.is_present("batch") {
            false => None,
            true => match matches.value_of("batch").unwrap() {
                "None" => None,
                z => Some(z),
            },
        };
        let eff_eps = match matches.is_present("eff_eps") {
            false => None,
            true => match matches.value_of("eff_eps").unwrap() {
                "None" => None,
                z => Some(z),
            },
        };

        // set learning rate later
        let boost_param = BoostParam::default()
            .set_boost_type(matches.value_of("boost_type").unwrap())
            .set_loss_func("logistic")
            .set_sample_weight_clip(clip_sample_weight)
            .set_sample_weight_wls_clip(clip_sample_wls_weight)
            .set_eps(eps_way)
            .set_is_dom_rec(is_dom_rec)
            .set_cov_way(cov_way)
            .set_batch_way(batch_way)
            .set_eff_eps(eff_eps);

        // TODO: better way in clap
        if matches.value_of("iter").is_some() & matches.value_of("iter_snv").is_some() {
            panic!("Do not indicate both --iter and --iter_snv.");
        } else if matches.value_of("iter").is_none() & matches.value_of("iter_snv").is_none() {
            panic!("Indicate either --iter or --iter_snv.");
        }

        // TODO: cleaner
        let boost_param = if matches.value_of("iter").is_some() {
            boost_param.set_iteration(matches.value_of("iter").unwrap().parse::<usize>().unwrap())
            //IterationNumber::Iteration(matches.value_of("iter").unwrap().parse::<usize>().unwrap())
        } else {
            boost_param.set_iteration_snv(
                matches
                    .value_of("iter_snv")
                    .unwrap()
                    .parse::<usize>()
                    .unwrap(),
            )
            //IterationNumber::Snv(
            //    matches
            //        .value_of("iter_snv")
            //        .unwrap()
            //        .parse::<usize>()
            //        .unwrap(),
            //)
        };

        boost_param.check();

        log::debug!("boost_param {:?}", boost_param);

        /*         // set learning rate later
        let boost_param = BoostParam::new_str(
            iteration,
            boost_type,
            "logistic",
            None,    //learning_rate,
            eps_way, //"medcase2allcell",
            clip_sample_weight,
            clip_sample_wls_weight,
            is_dom_rec,
            cov_way,
            batch_way,
            eff_eps,
        ); */

        let fin_snv = matches.value_of("file_snv").map(|x| PathBuf::from(x));
        let fin_wgt_cov = matches.value_of("file_wgt_cov").map(|x| PathBuf::from(x));
        let fin_loss = matches.value_of("file_loss").map(|x| PathBuf::from(x));

        let is_resume = matches.is_present("resume");

        let boost_method = if matches.is_present("pruning") {
            let clump_r2 = matches
                .value_of("clump_r2")
                .unwrap()
                .parse::<f64>()
                .unwrap();
            BoostMethod::Pruning(clump_r2)
        } else if matches.is_present("ss") {
            let clump_r2 = matches
                .value_of("clump_r2")
                .unwrap()
                .parse::<f64>()
                .unwrap();
            BoostMethod::Ss(clump_r2)
        } else {
            // classic boosting
            BoostMethod::Classic
        };

        log::info!("dout {:?}", dout);
        log::info!("file_plink {:?}", fin);
        log::info!("file_snv {:?}", fin_snv);
        log::info!("file_cov {:?}", fin_cov);
        log::info!("file_loss {:?}", fin_loss);
        log::info!("file_sample {:?}", fin_sample);
        log::info!("file_wgt_cov {:?}", fin_wgt_cov);
        log::info!("boost_param {:?}", boost_param);

        let use_adjloss = matches.is_present("use_adjloss");
        let use_const_for_loss = matches.is_present("use_const_for_loss");

        // Genoboost
        crate::boosting::run_boosting(
            &dout,
            &fin,
            fin_phe.as_deref(),
            phe_name.as_deref(),
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
            prune_snv,
            &learning_rates,
            is_monitor,
        );
    } else if let Some(ref matches) = matches.subcommand_matches("score") {
        let dout_score = PathBuf::from(matches.value_of("dir_score").unwrap());
        let fin = PathBuf::from(matches.value_of("file_plink").unwrap());
        let fin_phe = matches.value_of("file_phe").map(|x| PathBuf::from(x));
        let phe_name = match matches.value_of("phe").unwrap() {
            "None" => None,
            z => Some(z),
        };
        let fin_sample = matches.value_of("file_sample").map(|x| PathBuf::from(x));
        let fin_cov = matches.value_of("file_cov").map(|x| PathBuf::from(x));
        let dout_wgt = matches.value_of("dir_wgt").map(|x| PathBuf::from(x));
        let fout_wgt = matches.value_of("file_wgt").map(|x| PathBuf::from(x));
        //let dout_wgt = PathBuf::from(matches.value_of("dir_wgt").unwrap());
        // TODO: better way in clap
        if dout_wgt.is_some() & fout_wgt.is_some() {
            panic!("Do not indicate both --dir_wgt and --file_wgt.");
        } else if dout_wgt.is_none() & fout_wgt.is_none() {
            panic!("Indicate either --dir_wgt or --file_wgt.");
        }
        let use_iter = matches.is_present("use_iter");
        // FIXME: allow iter=None for fixed wgt
        let mut iterations: Vec<usize> = matches
            .values_of("iters")
            .unwrap()
            .map(|s| {
                s.parse::<usize>().unwrap_or_else(|| {
                    panic!("Iters should be able to be parsed to non-negative integer")
                })
            })
            .collect();
        // sort and deduplication
        iterations.sort();
        iterations.dedup();
        log::info!("iters {:?}", iterations);
        let learning_rates: Vec<Option<f64>> = matches
            .values_of("learning_rates")
            .unwrap()
            .map(|s| parse_arg_f64_none(s))
            .collect();

        // FIXME: remove later by getting boost_type from wgt config
        //let boost_param = BoostParam::new_type2();

        crate::boosting::run_boosting_score(
            &dout_score,
            &fin,
            fin_phe.as_deref(),
            phe_name.as_deref(),
            &iterations,
            dout_wgt.as_deref(), // use enum?
            fout_wgt.as_deref(),
            fin_cov.as_deref(),
            fin_sample.as_deref(),
            //boost_param,
            &learning_rates,
            use_iter,
        );
    }

    let end = start.elapsed();
    log::info!("It took {} seconds.", end.as_secs());
    log::info!("Done!!");
}
