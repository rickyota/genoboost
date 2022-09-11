//! Application of **Genoboost**.
//! Input plink file to run Genoboost.
// TODO: add examples
// TODO: exclude snv if all genotypes are the same
// TODO: usually do not consume -> if you have to clone() inside but original can be discarded, then use consume
// TODO: do not use `str: String` but `string: String`
// TODO: format of use_snvs, use_sample should be the same as plink --extract
// TODO: input several r2 at once
// TODO: ensure the same para when resuming
// TODO: (optional) write down extract snvs from top
// TODO: align A1, A2
// TODO: --verbose
// TODO: trimming sample weights: only use samples with large weights on choosing SNVs: Friedman, J., Hastie, T. and Tibshirani, R. (2000) ‘Additive logistic regression: a statistical view of boosting (With discussion and a rejoinder by the authors)’, Annals of statistics, 28(2), pp. 337–407. doi:10.1214/aos/1016218223.
// split main for genoboost_pruning, ssgenoboost

#[macro_use]
extern crate clap;
use boosting_rust::{self, BoostMethod, BoostParam};
use clap::{AppSettings, Arg, ArgGroup, ArgMatches, SubCommand};
use rayon;
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

fn get_matches() -> ArgMatches<'static> {
    let matches = app_from_crate!()
        .setting(AppSettings::DeriveDisplayOrder)
        .setting(AppSettings::SubcommandRequiredElseHelp) //require subcommand or show help
        .arg(Arg::from_usage("--num_threads [NUM] 'Number of threads.'").validator(check_usize))
        .subcommand( SubCommand::with_name("train").about("Train model.")
            //.subcommand( SubCommand::with_name("ss").about("Only use summary statistics."))
            .arg(Arg::from_usage("--resume 'Resume training.'"))
            .arg(Arg::from_usage("--write_loss 'Write down loss.'"))
            //.arg(Arg::from_usage("[write_loss] --write_loss 'Write down loss.'"))
            .arg(Arg::from_usage("--ss 'Only use summary statistics.'").requires("fin_loss"))
            .arg(Arg::from_usage("--pruning 'Use pruning method.'"))
            .group(ArgGroup::with_name("method").args(&["ss","pruning"]))
            // <FILE> : required
            // [FILE] : optional
            .arg(Arg::from_usage("--dir <FILE> 'Directory to output.'"))
            //.arg(Arg::from_usage("--fout <FILE> 'Prefix of output file.'"))
            .arg(Arg::from_usage("--file_plink <FILE> 'Prefix of a plink file. If those files are separated into chromosomes, set '%' where chromosome number is. If you use subcommand 'ss', indicate reference panel.'"))
            //.arg(Arg::from_usage("--fin [FILE] 'Prefix of a plink file. If those files are separated into chromosomes, set '*' where chromosome number is. If you use subcommand 'ss', indicate reference panel.'"))
            .arg(Arg::from_usage("--iter [NUM] 'Number of iterations.'").default_value("1000").validator(check_usize))
            .arg(
                Arg::from_usage("--boost_type [BOOSTTYPE] 'Boosting method to use.'")
                    .possible_values(&["genoboost","ada","constada","freemodelmissing"])
                    .default_value("genoboost"),
            )
            .arg(Arg::from_usage("--file_sample [FILE] 'File of samples to use.'"))
            .arg(Arg::from_usage("--file_snv [FILE] 'File of SNVs to use.'"))
            .arg(Arg::from_usage("--file_cov [FILE] 'File of covariates to use.'"))
            .arg(Arg::from_usage("--file_loss [FILE] 'File of loss function, which required only for --ss'"))
            .arg(Arg::from_usage("--learning_rates [VALS] 'Learning rates.'").default_value("None").multiple(true))
            //.arg(Arg::from_usage("--learning_rate [VAL] 'Learning rate.'").default_value("None").validator(check_f64_none))
            .arg(Arg::from_usage("--clip_sample_weight [VAL] 'Proportion and method of clipped samples.'").default_value("None")) // TODO: validator
            .arg(Arg::from_usage("--dom_and_rec 'Choose both dominant model and recessive model at one iteration.'"))
            .arg(Arg::from_usage("--prune_snv [VAL] 'Proportion of SNVs to be pruned from the candidate SNVs.'").default_value("None").validator(check_f64_none))
            .arg(Arg::from_usage("--clump_r2 [VAL] 'Threshold of clump correlation.'").default_value("0.1").validator(check_f64))
            // use .required_if() for r2 for pruning-genoboost
        )
        .subcommand( SubCommand::with_name("score").about("Calculate score.")
            .arg(Arg::from_usage("--dir_score <FILE> 'Prefix of output file.'"))
            //.arg(Arg::from_usage("--dout_score <FILE> 'Prefix of output file.'"))
            //.arg(Arg::from_usage("--fout <FILE> 'Prefix of output file.'"))
            .arg(Arg::from_usage("--file_plink <FILE> 'Prefix of a plink file. If those files are separated into chromosomes, set '*' where chromosome number is.'"))
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
            println!("Able to use SIMD.")
        } else {
            println!("Not able to use SIMD since avx2 is not detected.")
        }
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        println!("Not able to use SIMD since arch is not x86 or x86_64.")
    }

    let matches = get_matches();

    println!("matches {:?}", matches);

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
    println!("num_thread set: {}", rayon::current_num_threads());

    if let Some(ref matches) = matches.subcommand_matches("train") {
        let dout = PathBuf::from(matches.value_of("dir").unwrap());
        //let fout = matches.value_of("fout").unwrap();
        let fin = PathBuf::from(matches.value_of("file_plink").unwrap());
        let fin_sample = matches.value_of("file_sample").map(|x| PathBuf::from(x));
        let fin_cov = matches.value_of("file_cov").map(|x| PathBuf::from(x));
        //let fin_cov = matches.value_of("fin_cov");
        let iteration = matches.value_of("iter").unwrap().parse::<usize>().unwrap();
        let boost_type = matches.value_of("boost_type").unwrap();
        let learning_rates: Vec<Option<f64>> = matches
            .values_of("learning_rates")
            .unwrap()
            .map(|s| parse_arg_f64_none(s))
            .collect();
        //let learning_rate = parse_arg_f64_none(matches.value_of("learning_rate").unwrap());
        let prune_snv = parse_arg_f64_none(matches.value_of("prune_snv").unwrap());

        let is_write_loss = matches.is_present("write_loss");
        //let is_write_loss = true;

        // make this Option<> here? or in Boostingparam?
        // should be here, or there should eb option to make Option<> as default
        let clip_sample_weight = match matches.value_of("clip_sample_weight").unwrap() {
            "None" => None,
            z => Some(z),
        };

        let is_dom_rec = matches.is_present("dom_and_rec");

        // set learning rate later
        let boost_param = BoostParam::new_str(
            iteration,
            boost_type,
            "logistic",
            None,
            //learning_rate,
            "medcase2allcell",
            clip_sample_weight,
            is_dom_rec,
        );

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

        println!("dout {:?}", dout);
        println!("file_plink {:?}", fin);
        println!("file_snv {:?}", fin_snv);
        println!("file_cov {:?}", fin_cov);
        println!("file_loss {:?}", fin_loss);
        println!("file_sample {:?}", fin_sample);
        println!("file_wgt_cov {:?}", fin_wgt_cov);
        println!("boost_param {:?}", boost_param);

        // Genoboost
        boosting_rust::run_boosting(
            &dout,
            &fin,
            iteration,
            boost_method,
            boost_param,
            fin_snv.as_deref(),
            fin_sample.as_deref(),
            fin_cov.as_deref(),
            is_resume,
            is_write_loss,
            prune_snv,
            &learning_rates,
        );
    } else if let Some(ref matches) = matches.subcommand_matches("score") {
        let dout_score = PathBuf::from(matches.value_of("dir_score").unwrap());
        let fin = PathBuf::from(matches.value_of("file_plink").unwrap());
        let fin_sample = matches.value_of("file_sample").map(|x| PathBuf::from(x));
        let fin_cov = matches.value_of("file_cov").map(|x| PathBuf::from(x));
        let dout_wgt = matches.value_of("dir_wgt").map(|x| PathBuf::from(x));
        let fout_wgt = matches.value_of("file_wgt").map(|x| PathBuf::from(x));
        //let dout_wgt = PathBuf::from(matches.value_of("dir_wgt").unwrap());
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
                s.parse::<usize>()
                    .expect("Iters should be able to be parsed to non-negative integer")
            })
            .collect();
        // sort and deduplication
        iterations.sort();
        iterations.dedup();
        println!("iters {:?}", iterations);
        let learning_rates: Vec<Option<f64>> = matches
            .values_of("learning_rates")
            .unwrap()
            .map(|s| parse_arg_f64_none(s))
            .collect();

        // FIXME: remove later by getting boost_type from wgt config
        let boost_param = BoostParam::new_type2();

        boosting_rust::run_boosting_score(
            &dout_score,
            &fin,
            &iterations,
            dout_wgt.as_deref(), // use enum?
            fout_wgt.as_deref(),
            fin_cov.as_deref(),
            fin_sample.as_deref(),
            boost_param,
            &learning_rates,
            use_iter,
        );
    }

    let end = start.elapsed();
    println!("It took {} seconds.", end.as_secs());
    println!("Done!!");
}
