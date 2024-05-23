// use AsRef<Path>
// https://www.reddit.com/r/rust/comments/7mu7q1/is_working_with_paths_always_this_painful/
// should use &Path not &Path

use super::WgtBoost;
use super::{Coef, Model, Wgt, WgtKind};
use crate::{BoostType, DoutParaFile};
use genetics::{textfile, vec, SnvId, SnvInteractionWgt, SnvWgt};
use std::collections::HashMap;
//use std::fs;
//use std::fs::File;
//use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::Path;

//pub fn get_dir_score(dout: &Path, learning_rate: f64) -> PathBuf {
//    get_dname_para(dout, learning_rate)
//}
//
//pub fn get_dir_cv(dout: &Path, cvi: usize) -> PathBuf {
//    let mut d = dout.to_owned();
//    let dpara = String::from("cv-") + &cvi.to_string();
//    d.push(dpara);
//    d
//}
//
//pub fn get_file_wgt(dout: &Path, learning_rate: f64) -> PathBuf {
//    let p = get_dname_para(dout, learning_rate);
//    get_fname_wgt(&p)
//}
///* pub fn get_file_wgt(dout: &Path, learning_rate: &Option<f64>) -> PathBuf {
//    let p = get_dname_para(dout, learning_rate);
//    get_fname_wgt(&p)
//} */
//
//// best across snvs and lr
//pub fn get_file_wgt_best(dout: &Path) -> PathBuf {
//    get_fname_wgt(dout)
//}
//
//// store best snvs and lr
//pub fn get_file_wgt_best_para(dout: &Path) -> PathBuf {
//    get_fname_wgt_best_para(dout)
//}

// when increase para, add them
//pub fn get_dname_para(dout: &Path, learning_rate: f64) -> PathBuf {
//    // TODO: use .join
//    let mut d = dout.to_owned();
//    let mut dpara = String::from("para");
//    //if let Some(lr) = learning_rate {
//    dpara += &(String::from("_lr-") + &learning_rate.to_string());
//    //dpara += &(String::from(".lr") + &learning_rate.to_string());
//    //}
//    d.push(dpara);
//    d
//}

//pub fn get_dname_para_integrate(dout: &Path, learning_rate: f64, boost_type: BoostType) -> PathBuf {
//    // TODO: use .join
//    let mut d = dout.to_owned();
//    let mut dpara = String::from("para");
//    dpara += &(String::from("_type-") + &boost_type.genetic_model());
//    //if let Some(lr) = learning_rate {
//    dpara += &(String::from("_lr-") + &learning_rate.to_string());
//    //dpara += &(String::from(".lr") + &learning_rate.to_string());
//    //}
//    d.push(dpara);
//    d
//}

//pub fn get_dname_wgt(dout: &Path) -> PathBuf {
//    let d = dout.to_owned();
//    //let dwgt = "wgt/";
//    //d.push(dwgt);
//    d
//}

// moved to DoutFile::get_file_wgt()
// TODO: concat to get_file_wgt()
//pub fn get_fname_wgt(dout_para: &Path) -> PathBuf {
//    let mut fwgt = get_dname_wgt(dout_para);
//
//    fwgt.push("boosting.wgt");
//    //fwgt.set_file_name("boosting.wgt");
//    fwgt
//    //fout.to_owned() + ".wgt"
//}

//pub fn get_fname_wgt_best_para(dout_para: &Path) -> PathBuf {
//    let mut fwgt = get_dname_wgt(dout_para);
//
//    fwgt.push("boosting.wgt.para");
//    //fwgt.set_file_name("boosting.wgt");
//    fwgt
//    //fout.to_owned() + ".wgt"
//}

/* pub fn get_fname_wgt(fout: &str) -> String {
    fout.to_owned() + ".wgt"
} */

//pub fn create_dir(dout: &Path) {
//    fs::create_dir_all(&dout).unwrap();
//}
//
//pub fn get_dname_loss(dout: &Path) -> PathBuf {
//    let mut f = dout.to_owned();
//    let dscore = "loss/";
//    f.push(dscore);
//    f
//}
//
//pub fn get_fname_loss(dout: &Path, ti: usize) -> PathBuf {
//    let mut f = get_dname_loss(dout);
//    let fname = "iter-".to_string() + &ti.to_string() + ".loss";
//    f.push(fname);
//    f
//    //fout.to_owned() + ".iter" + &ti.to_string() + ".loss"
//}
//
//pub fn get_fname_loss_adjmax(dout: &Path, ti: usize) -> PathBuf {
//    let mut f = get_dname_loss(dout);
//    let fname = "iter-".to_string() + &ti.to_string() + ".adjmax.loss";
//    f.push(fname);
//    f
//    //fout.to_owned() + ".iter" + &ti.to_string() + ".loss"
//}
//
//pub fn get_fname_acc(dout: &Path) -> PathBuf {
//    let mut f = dout.to_owned();
//    let fname = "monitor.acc".to_string();
//    f.push(fname);
//    f
//}

/* pub fn get_fname_ws(fout: &str) -> String {
    fout.to_owned() + ".ws"
} */

// moved to DoutFile
//pub fn check_file_wgt_not_exist(dout: &Path, learning_rate: f64) {
//    let fwgt = get_file_wgt(dout, learning_rate);
//    //log::debug!("fwgt {:?}", fwgt);
//    let exist_fwgt = textfile::exist_file(&fwgt);
//    if exist_fwgt {
//        panic!(
//            "Weight file already exists: {:?}. Delete it or use --resume option.",
//            &fwgt
//        );
//    }
//}

/* pub fn check_file_wgt_not_exist(dout: &Path, learning_rate: &Option<f64>) {
    let fwgt = get_file_wgt(dout, learning_rate);
    let exist_fwgt = textfile::exist_file(&fwgt);
    if exist_fwgt {
        panic!(
            "Weight file already exists: {:?}. Delete it or use --resume option.",
            &fwgt
        );
    }
} */

//pub fn check_file_wgt_not_exist_integrate(dout: &Path) {
//    let fwgt = get_fname_wgt(dout);
//    log::debug!("fwgt {:?}", fwgt);
//    let exist_fwgt = textfile::exist_file(&fwgt);
//    if exist_fwgt {
//        panic!(
//            "Weight file already exists: {:?}. Delete it or use --resume option.",
//            &fwgt
//        );
//    }
//}

//pub fn check_file_wgt_exist(fwgt: &Path) {
//    let exist_fwgt = textfile::exist_file(&fwgt);
//    if !exist_fwgt {
//        panic!("fwgt does not exist: {:?}.", &fwgt);
//    }
//}

//#[allow(dead_code)]
//pub fn check_file_wgt_exist_dir(dout: &Path) {
//    let fwgt = get_fname_wgt(dout);
//    check_file_wgt_exist(&fwgt)
//}

//pub fn bufwriter_floss(dout: &Path, ti: usize) -> BufWriter<File> {
//    let dloss = get_dname_loss(dout);
//    create_dir(&dloss);
//    let floss = get_fname_loss(dout, ti);
//    //if text::exist_file(&floss) {
//    //    panic!("File already exists: {}.", fout);
//    //}
//
//    let file = match File::create(&floss) {
//        Ok(file) => file,
//        Err(_) => panic!(
//            "Cannot create file, possibly directory does not exist: {:?}",
//            &floss
//        ),
//    };
//
//    BufWriter::new(file)
//    /*
//    let mut writer = BufWriter::new(file);
//    let str =
//        "iteration\tkind\tvar\tmodel\tthreshold\talpha\tconst\tchrom\tpos\ta1\ta2\n".to_owned();
//    writer.write(str.as_bytes()).unwrap();
//     */
//}

//pub fn bufwriter_floss_adjmax(dout: &Path, ti: usize) -> BufWriter<File> {
//    let dloss = get_dname_loss(dout);
//    create_dir(&dloss);
//    let floss = get_fname_loss_adjmax(dout, ti);
//    //if text::exist_file(&floss) {
//    //    panic!("File already exists: {}.", fout);
//    //}
//
//    let file = match File::create(&floss) {
//        Ok(file) => file,
//        Err(_) => panic!(
//            "Cannot create file, possibly directory does not exist: {:?}",
//            &floss
//        ),
//    };
//
//    BufWriter::new(file)
//    /*
//    let mut writer = BufWriter::new(file);
//    let str =
//        "iteration\tkind\tvar\tmodel\tthreshold\talpha\tconst\tchrom\tpos\ta1\ta2\n".to_owned();
//    writer.write(str.as_bytes()).unwrap();
//     */
//}

//pub fn bufwriter_acc(dout: &Path) -> BufWriter<File> {
//    let facc = get_fname_acc(dout);
//
//    let file = match File::create(&facc) {
//        Ok(file) => file,
//        Err(_) => panic!(
//            "Cannot create file, possibly directory does not exist: {:?}",
//            &facc
//        ),
//    };
//
//    BufWriter::new(file)
//}

//pub fn is_exist_wgt(dout: &Path) -> bool {
//    let fwgt = get_fname_wgt(dout);
//    textfile::exist_file(&fwgt)
//}

//pub fn is_nonzero(dout: &Path) -> bool {
//    let fwgt = get_fname_wgt(dout);
//
//    match File::open(&fwgt) {
//        Ok(f) => f.metadata().unwrap().len() != 0,
//        Err(_) => false,
//    }
//}

//pub fn is_exist_wgt_nonzero(dout: &Path) -> bool {
//    is_exist_wgt(dout) && is_nonzero(dout)
//}

/*
pub fn bufwriter_file_wgt(fout: &str) -> (File, bool) {
    let fwgt = get_fname_wgt(fout);
    //if text::exist_file(&fwgt) {
    //    panic!("File already exists: {}.", fout);
    //}
    let is_exist = text::exist_file(&fwgt);

    let file = match OpenOptions::new()
        .read(true)
        .append(true)
        .create(true)
        .open(&fout)
    {
        Ok(file) => file,
        Err(_) => panic!("Cannot open file: {}", &fout),
    };

    //let file = match File::create(&fwgt) {
    //    Ok(file) => file,
    //    Err(_) => panic!(
    //        "Cannot create file, possibly directory does not exist: {}",
    //        &fwgt
    //    ),
    //};
    // TO append
    // https://masahiko-ofgp-ja-notebook.blogspot.com/2019/05/rust-lang-file-append-mode.html
    //let file = OpenOptions::new().append(true).open(&fwgt).unwrap();

    (file, is_exist)
}
*/

/*
pub fn bufwriter_file(file: File) -> BufWriter<File> {
    BufWriter::new(file)
}
*/

/* pub fn bufwriter_fwgt(dout: &Path) -> BufWriter<File> {
    let fwgt = get_fname_wgt(dout);
    if text::exist_file(&fwgt) {
        panic!("File already exists: {:?}.", fwgt);
    }

    let file = match File::create(&fwgt) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {:?}",
            &fwgt
        ),
    };
    // TO append
    // https://masahiko-ofgp-ja-notebook.blogspot.com/2019/05/rust-lang-file-append-mode.html
    //let file = OpenOptions::new().append(true).open(&fwgt).unwrap();

    BufWriter::new(file)
    /*
    let mut writer = BufWriter::new(file);
    let str =
        "iteration\tkind\tvar\tmodel\tthreshold\talpha\tconst\tchrom\tpos\ta1\ta2\n".to_owned();
    writer.write(str.as_bytes()).unwrap();
     */
} */

//pub fn bufwriter_fwgt_append(dout: &Path) -> BufWriter<File> {
//    // create dwgt
//    //let dwgt = get_dname_wgt(dout);
//    //create_dir(&dwgt);
//
//    let fwgt = get_fname_wgt(dout);
//    // allow exist for resume
//    //if text::exist_file(&fwgt) {
//    //    panic!("File already exists: {}.", fwgt);
//    //}
//
//    let file = match OpenOptions::new()
//        .read(true) // unnecessary
//        .append(true)
//        .create(true)
//        .open(&fwgt)
//    {
//        Ok(file) => file,
//        Err(_) => panic!("file does not exist: {:?}", &fwgt),
//    };
//
//    // TO append
//    // https://masahiko-ofgp-ja-notebook.blogspot.com/2019/05/rust-lang-file-append-mode.html
//    //let file = OpenOptions::new().append(true).open(&fwgt).unwrap();
//
//    BufWriter::new(file)
//}

/* pub fn bufreader_fwgt(dout: &Path) -> BufReader<File> {
    let fwgt = get_fname_wgt(dout);
    if !text::exist_file(&fwgt) {
        panic!("File does not exist: {:?}.", fwgt);
    }

    BufReader::new(File::open(fwgt).unwrap())
} */

pub fn write_cols<W: std::io::Write>(writer: &mut BufWriter<W>, boost_type: BoostType) {
    //let str =
    //    "iteration\tkind\tvar\tmodel\tthreshold\talpha\tconst\tchrom\tpos\ta1\ta2\n".to_owned();
    let strings = wgt_columns(boost_type);
    let str = strings.join("\t") + "\n";
    writer.write(str.as_bytes()).unwrap();
}

/* pub fn init_fwgt(fout: &str) {
    let fwgt = get_fname_wgt(fout);
    if text::exist_file(&fwgt) {
        panic!("File already exists: {}.", fout);
    }

    let file = match File::create(&fwgt) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {}",
            &fwgt
        ),
    };

    let mut writer = BufWriter::new(file);
    //let mut writer = BufWriter::new(File::create(fwgt).unwrap());
    let str =
        "iteration\tkind\tvar\tmodel\tthreshold\talpha\tconst\tchrom\tpos\ta1\ta2\n".to_owned();
    writer.write(str.as_bytes()).unwrap();
} */

pub fn wgt_columns(boost_type: BoostType) -> Vec<String> {
    let mut cols_boost = match boost_type {
        // score is the same as alpha, const
        BoostType::Ada | BoostType::ConstAda => vec!["threshold"],
        BoostType::FreeModelMissing | BoostType::Logit => {
            vec!["score0", "score1", "score2", "scorem"]
        }
        BoostType::LogitNoMissing => {
            vec!["score0", "score1", "score2"]
        }
        BoostType::LogitMhcNoMissing | BoostType::LogitCommon => {
            vec!["score0", "score1", "score2"]
        }
        BoostType::LogitAdd => {
            // alpha, const
            // same cols for cov
            vec![]
        }
        BoostType::LogitAddInteraction => {
            // alpha, const
            vec![]
        }
    };

    let mut cols_snv_2 = match boost_type {
        BoostType::LogitAddInteraction => {
            vec!["var_2", "chrom_2", "pos_2", "a1_2", "a2_2", "a2_2_frq"]
        }
        _ => vec![],
    };

    let mut columns = vec![
        "iteration",
        "kind",
        "var",
        "model", //"threshold",
        "eps",
        "eff_eps",
        "alpha",
        "const",
    ];
    columns.append(&mut cols_boost);

    let mut cols_snv = vec!["chrom", "pos", "a1", "a2", "a2_frq"];
    //let mut cols_snv = vec!["chrom", "pos", "a1", "a2"];
    columns.append(&mut cols_snv);

    columns.append(&mut cols_snv_2);
    columns.iter().map(|x| (*x).to_owned()).collect()
}

/* pub fn content_to_col(wgt_kind: &WgtKind, boost_type: BoostType) -> HashMap<String, String> {
    let column_convert = HashMap::from([
        ("alpha".to_owned(), "alpha".to_owned()),
        ("const".to_owned(), "const".to_owned()),
    ]);
    column_convert
} */

/* /// Even if fws exists, overwrite.
pub fn init_fws(fout: &str) {
    let fws = get_fname_ws(fout);
    let file = match File::create(&fws) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {}",
            &fws
        ),
    };
    let mut writer = BufWriter::new(file);
    //let mut writer = BufWriter::new(File::create(fws).unwrap());
    let str = "".to_owned();
    writer.write(str.as_bytes()).unwrap();
} */

/* pub fn write_ws(fout: &str, ws: &[f64], iter: usize) {
    let fws = get_fname_ws(fout);

    let file = OpenOptions::new().append(true).open(fws).unwrap();
    let mut writer = BufWriter::new(file);

    let str = ws
        .iter()
        .map(|s| format!("{:.2e}", s))
        .collect::<Vec<String>>()
        .join(" ");
    //let str = "a\tb\n".to_owned();
    let str: String = iter.to_string() + " " + &str + "\n";

    writer.write(str.as_bytes()).unwrap();
} */

/* fn create_var_name_to_coef(fin_wgt_cov: &Path) -> HashMap<String, f64> {
    let mut var_name_to_coef = HashMap::new();

    let cols = [0usize, 1];
    let vss: Vec<Vec<String>> = text::load_table_cols(fin_wgt_cov, &cols, true).unwrap();

    for var_i in 0..vss[0].len() {
        let coef = vss[1][var_i].parse::<f64>().unwrap();
        var_name_to_coef.insert(vss[0][var_i].to_owned(), coef);
    }

    var_name_to_coef
} */

/*
/// consume covs
//pub fn load_wgt_cov(fin_wgt_cov: Option<&str>, covs: Option<&[Var]>, n: usize) -> Vec<WgtBoost> {
pub fn load_wgt_cov(fin_wgt_cov: Option<&str>, covs: Option<Vec<Var>>, n: usize) -> Vec<WgtBoost> {
    if fin_wgt_cov.is_none() || covs.is_none() {
        let wgts_cov: Vec<WgtBoost> = Vec::new();
        //let wgts_cov: Vec<Wgt> = Vec::new();
        return wgts_cov;
    }
    /*
    if fin_wgt_cov.is_none() || covs.len() == 0 {
        let wgts_cov: Vec<Wgt> = Vec::new();
        return wgts_cov;
    }
    */

    let covs = covs.unwrap();

    let var_name_to_coef = create_var_name_to_coef(fin_wgt_cov.unwrap());

    log::debug!("var_name_to_coef {:?}", var_name_to_coef);

    let mut wgts_cov: Vec<WgtBoost> = Vec::new();

    let iteration = 0;
    let coef = *var_name_to_coef
        .get("const")
        .expect("No 'const' row in fin_wgt_cov");

    let wgt_const = Wgt::construct_const_linear(Coef::Linear(coef));
    /*
    // set as Cov since this const is in fin_cov
    let mut cov = Var::construct_var(CovKind::Cov, "const".to_owned());
    //let mut cov = Var::construct_var(VarKind::Const, "const".to_owned());
    cov.set_vals_const(n);
    let cov_wgt = CovWgt::construct(cov);
    let wgt_const = Wgt::construct_wgt(
        WgtKind::Cov(cov_wgt),
        Model::new_linear(coef)
        //Model::Linear(LinearModel::new(Coef::Linear(coef))),
        //iter,
        //f64::NAN,
        //(f64::NAN, f64::NAN, f64::NAN, f64::NAN),
    );
     */

    let wgt_boost_const = WgtBoost::construct_wgt_iteration(wgt_const, iteration);

    /*
    // first add const
    let const_wgt = ConstWgt::construct_const_wgt(LinearModel::construct(coef));
    let wgt_const = Wgt::construct_wgt(
        WgtKind::Const(const_wgt),
        iter,
        f64::NAN,
        (f64::NAN, f64::NAN, f64::NAN, f64::NAN),
    );
    */
    wgts_cov.push(wgt_boost_const);

    for (wgt_i, cov) in covs.into_iter().enumerate() {
        let iteration = wgt_i + 1;

        let coef = *var_name_to_coef
            .get(cov.name())
            .expect(&format!("No '{}' row in fin_wgt_cov", cov.name()));

        let cov_wgt = CovWgt::construct(cov);
        let wgt_cov = Wgt::construct_wgt(
            WgtKind::Cov(cov_wgt),
            Model::new_linear(coef)
            //Model::Linear(LinearModel::new(Coef::Linear(coef))),
            //iter,
            //f64::NAN,
            //(f64::NAN, f64::NAN, f64::NAN, f64::NAN),
        );
        let wgt_boost_cov = WgtBoost::construct_wgt_iteration(wgt_cov, iteration);
        wgts_cov.push(wgt_boost_cov);
    }

    log::debug!("wgts_cov: {:?}", wgts_cov);

    wgts_cov
}
     */

//#[allow(dead_code)]
//fn load_num_iteration_dir(dout: &Path) -> usize {
//    let fwgt = get_fname_wgt(dout);
//    let iter_col = 0;
//    let iters = textfile::load_table_col(&fwgt, iter_col, true);
//    let iters: Vec<usize> = iters.iter().map(|s| s.parse::<usize>().unwrap()).collect();
//    let num_iteration = iters.iter().max().unwrap() + 1;
//    num_iteration
//}

fn load_num_iteration(fwgt: &Path) -> usize {
    //let fwgt = get_fname_wgt(dout);
    let iter_col = 0;
    let iters = textfile::load_table_col(&fwgt, iter_col, true);
    let iters: Vec<usize> = iters.iter().map(|s| s.parse::<usize>().unwrap()).collect();
    let num_iteration = iters.iter().max().unwrap() + 1;
    num_iteration
}

//#[allow(dead_code)]
//pub fn valid_iterations_dir(iterations_in: &[usize], dout: &Path) -> Vec<usize> {
//    //let fwgt = get_fname_wgt(dout);
//
//    let mut iterations = iterations_in.to_vec();
//    iterations.sort();
//
//    let num_iteration = load_num_iteration_dir(dout);
//
//    let iterations = iterations
//        .into_iter()
//        .filter(|&v| v <= num_iteration)
//        .collect();
//
//    iterations
//}

pub fn valid_iterations(iterations_in: &[usize], fwgt: &Path) -> Vec<usize> {
    //let fwgt = get_fname_wgt(dout);

    let mut iterations = iterations_in.to_vec();
    iterations.sort();

    let num_iteration = load_num_iteration(fwgt);

    let iterations = iterations
        .into_iter()
        .filter(|&v| v <= num_iteration)
        .collect();

    iterations
}

//pub fn load_wgts(dout: &Path, boost_type: BoostType) -> Vec<WgtBoost> {
//pub fn load_wgts(dout: &Path) -> Vec<WgtBoost> {
pub fn load_wgts(dout: &DoutParaFile) -> Vec<WgtBoost> {
    //let fwgt = get_fname_wgt(dout);
    let fwgt = dout.get_file_wgt();

    load_wgts_file(&fwgt)
    //load_wgts_file(&fwgt, boost_type)
}

/// if v1 is in v
//pub fn vector_issubset(v1: &[&str], v: &[&str]) -> bool {
//    v1.iter().map(|x| v.contains(x)).all(|x| x)
//}

// coef_type could be different among rows
// judge from model column in each row
//pub fn load_wgts_boost_coef_type(fwgt: &Path) -> Coef {
//    let header = textfile::load_table_header(fwgt);
//    let header_str: Vec<&str> = header.iter().map(AsRef::as_ref).collect();
//
//    let cols_score4 = ["score0", "score1", "score2", "scorem"];
//
//    if vector_issubset(&cols_score4, &header_str) {
//        return Coef::Score4((f64::NAN, f64::NAN, f64::NAN, f64::NAN));
//    }
//
//    let cols_score3 = ["score0", "score1", "score2"];
//
//    if vector_issubset(&cols_score3, &header_str) {
//        return Coef::Score3((f64::NAN, f64::NAN, f64::NAN));
//    }
//
//    let cols_score2 = ["score0", "score1"];
//
//    if vector_issubset(&cols_score2, &header_str) {
//        return Coef::Score2((f64::NAN, f64::NAN));
//    }
//
//    // TODO: others
//
//    // alpha, const is also for cov
//    // check if all NAN for snvs
//
//    // const, alpha
//    // LogitAdd
//    return Coef::LinearConst((f64::NAN, f64::NAN));
//}

pub fn non_nan_all(
    cols: &[&str],
    vss: &Vec<Vec<String>>,
    col_to_i: &HashMap<String, usize>,
    wgt_i: usize,
) -> bool {
    let is_exist = cols.iter().map(|x| col_to_i.contains_key(*x)).all(|x| x);

    if !is_exist {
        return false;
    }

    let is_non_nan = cols
        .iter()
        .map(|x| vss[col_to_i[*x]][wgt_i].parse::<f64>().unwrap().is_nan())
        .all(|x| !x);

    is_non_nan
}

// TODO: binary of cov and linear of snv
//pub fn load_wgts_file(fwgt: &Path, boost_type: BoostType) -> Vec<WgtBoost> {
pub fn load_wgts_file(fwgt: &Path) -> Vec<WgtBoost> {
    //let fwgt = get_fname_wgt(dout);

    // also load header
    let vss: Vec<Vec<String>> = textfile::load_table(fwgt, false);

    load_wgts_file_string(&vss)

    //let col_n = vss.len();
    //let mut columns: Vec<String> = Vec::new();
    //for col_i in 0..col_n {
    //    columns.push(vss[col_i][0].clone());
    //}

    //let mut col_to_i: HashMap<String, usize> = HashMap::new();
    //for col_i in 0..col_n {
    //    col_to_i.insert(columns[col_i].clone(), col_i);
    //}
    //log::debug!("col_to_i {:?}", col_to_i);

    //let mut wgts: Vec<WgtBoost> = Vec::new();
    //for wgt_i in 1..vss[0].len() {
    //    // FIXME: when running without deleting wgt file beforehand,
    //    // vss[col_to_i["iteration"]][wgt_i] could be "iteration" and .parse() would fail
    //    assert_eq!(
    //        wgt_i - 1,
    //        (vss[col_to_i["iteration"]][wgt_i].parse::<usize>().unwrap()),
    //        "Number of iterator in wgt is wrong at iteration {} but found {}.",
    //        wgt_i - 1,
    //        (vss[col_to_i["iteration"]][wgt_i].parse::<usize>().unwrap()),
    //    );

    //    // kind
    //    let wgt_kind = &vss[col_to_i["kind"]][wgt_i];
    //    let model_name = &vss[col_to_i["model"]][wgt_i];
    //    let wgt = if wgt_kind == "CONST" {
    //        //const
    //        //let threshold = vss[4][wgt_i].parse::<f64>().unwrap();
    //        // should not use?
    //        let threshold = 0.0;

    //        //let name = vss[col_to_i["var"]][wgt_i].clone();
    //        //let const_var = Var::construct_var(CovKind::Const, name);
    //        // later
    //        //const_var.set_vals_const(n);
    //        // FIXME: how to use CovWgt?
    //        //let const_wgt = CovWgt::construct(const_var);

    //        let model = if model_name == "LINEAR" {
    //            let coef = vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap();
    //            Model::new_linear(coef)
    //        } else if model_name == "BINARY" {
    //            let coef = Coef::Binary((
    //                vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
    //                vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
    //            ));
    //            Model::new(coef, Some(threshold))
    //            //Model::Binary(BinaryModel::new(coef, threshold))
    //        } else {
    //            panic!("Unknown model.");
    //        };
    //        /*
    //        let model = if threshold.is_nan() {
    //            let coef = vss[5][wgt_i].parse::<f64>().unwrap();
    //            Model::Linear(LinearModel::construct(coef))
    //        } else {
    //            Model::Binary(BinaryModel::construct_threshold(threshold))
    //        };
    //        */

    //        let wgt = Wgt::construct_const(model);
    //        //let wgt = Wgt::construct_wgt(WgtKind::Cov(const_wgt), model);
    //        //let wgt = Wgt::construct_wgt(WgtKind::Cov(const_wgt), model, wgt_i);
    //        let wgt_boost = WgtBoost::construct_wgt_iteration(wgt, wgt_i);
    //        wgt_boost
    //    } else if wgt_kind == "COV" {
    //        // cov
    //        // only implement for linear
    //        let name = vss[col_to_i["var"]][wgt_i].clone();
    //        let coef = vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap();
    //        let wgt = Wgt::construct_cov_name(name, Coef::Linear(coef));
    //        /*
    //        let cov_var = Var::construct_var(CovKind::Cov, name.clone());
    //        let cov_wgt = CovWgt::construct(cov_var);
    //        let coef = vss[5][wgt_i].parse::<f64>().unwrap();
    //        let model = Model::new_linear(coef);
    //        let wgt = Wgt::construct_wgt(WgtKind::Cov(cov_wgt), model);
    //         */
    //        let wgt_boost = WgtBoost::construct_wgt_iteration(wgt, wgt_i);
    //        wgt_boost
    //    } else if wgt_kind == "SNV" {
    //        let snv = SnvId::new_snv_index(
    //            vss[col_to_i["var"]][wgt_i].clone(),
    //            &vss[col_to_i["chrom"]][wgt_i],
    //            &vss[col_to_i["pos"]][wgt_i],
    //            vss[col_to_i["a1"]][wgt_i].clone(),
    //            vss[col_to_i["a2"]][wgt_i].clone(),
    //        );
    //        //let snv_wgt = SnvWgt::construct_score(snv);

    //        let model = match boost_coef_type {
    //            Coef::Binary(_) => {
    //                let coef = Coef::Binary((
    //                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
    //                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
    //                ));
    //                let threshold = vss[col_to_i["threshold"]][wgt_i].parse::<f64>().unwrap();
    //                let model = Model::new(coef, Some(threshold));
    //                model
    //            }
    //            Coef::LinearConst(_) => {
    //                let coef = Coef::LinearConst((
    //                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
    //                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
    //                ));
    //                let model = Model::new_coef(coef);
    //                model
    //            }
    //            //BoostType::FreeModelMissing | BoostType::Logit => {
    //            Coef::Score4(_) => {
    //                let coef = Coef::Score4((
    //                    vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
    //                    vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
    //                    vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
    //                    vss[col_to_i["scorem"]][wgt_i].parse::<f64>().unwrap(),
    //                ));
    //                let model = Model::new_coef(coef);
    //                model
    //            }
    //            Coef::Score3(_) => {
    //                let coef = Coef::Score3((
    //                    vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
    //                    vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
    //                    vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
    //                    //vss[col_to_i["scorem"]][wgt_i].parse::<f64>().unwrap(),
    //                ));
    //                let model = Model::new_coef(coef);
    //                model
    //            }
    //            _ => {
    //                unimplemented!();
    //            }
    //        };

    //        //let model = match boost_type {
    //        //    BoostType::ConstAda => {
    //        //        let coef = Coef::Binary((
    //        //            vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
    //        //            vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
    //        //        ));
    //        //        let threshold = vss[col_to_i["threshold"]][wgt_i].parse::<f64>().unwrap();
    //        //        let model = Model::new(coef, Some(threshold));
    //        //        model
    //        //    }
    //        //    BoostType::FreeModelMissing | BoostType::Logit => {
    //        //        let coef = Coef::Score4((
    //        //            vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
    //        //            vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
    //        //            vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
    //        //            vss[col_to_i["scorem"]][wgt_i].parse::<f64>().unwrap(),
    //        //        ));
    //        //        let model = Model::new_coef(coef);
    //        //        model
    //        //    }
    //        //    BoostType::Ada => {
    //        //        unimplemented!();
    //        //    }
    //        //};

    //        //let model = Model::Binary(BinaryModel::new(coef, threshold));
    //        let wgt = Wgt::construct_wgt(WgtKind::Snv(snv, None, None), model);
    //        let wgt_boost = WgtBoost::construct_wgt_iteration(wgt, wgt_i);
    //        wgt_boost
    //    } else {
    //        panic!("Unknown wgt kind: {}", vss[1][wgt_i]);
    //    };

    //    wgts.push(wgt);
    //}

    //// check iters are sorted.
    //assert!(vec::is_sorted(
    //    &wgts.iter().map(|v| v.iteration()).collect::<Vec<usize>>()
    //));
    // /*
    //// is_sorted() is nightly
    //assert!(wgts
    //    .iter()
    //    .map(|v| v.get_iter())
    //    .collect::<Vec<usize>>()
    //    .is_sorted());
    //     */

    //wgts
}

pub fn load_wgts_file_string(vss: &Vec<Vec<String>>) -> Vec<WgtBoost> {
    let col_n = vss.len();
    let mut columns: Vec<String> = Vec::new();
    for col_i in 0..col_n {
        columns.push(vss[col_i][0].clone());
    }

    let mut col_to_i: HashMap<String, usize> = HashMap::new();
    for col_i in 0..col_n {
        col_to_i.insert(columns[col_i].clone(), col_i);
    }
    log::debug!("col_to_i {:?}", col_to_i);

    let wgt_ver = judge_wgt_version(&vss, &col_to_i);
    log::debug!("wgt_ver {:?}", wgt_ver);
    let load_wgt_func = match wgt_ver.as_str() {
        "ver1" => load_wgt_row_ver1,
        "ver2" => load_wgt_row_ver2,
        _ => unreachable!(),
    };

    let mut wgts: Vec<WgtBoost> = Vec::new();
    for row_i in 1..vss[0].len() {
        // FIXME: when running without deleting wgt file beforehand,
        // vss[col_to_i["iteration"]][wgt_i] could be string "iteration" and .parse() would fail
        if vss[col_to_i["iteration"]][row_i].parse::<usize>().is_err() {
            panic!(
                "Iteration number at row {} is not usize: {}",
                row_i, vss[col_to_i["iteration"]][row_i]
            );
        }
        if row_i - 1 != vss[col_to_i["iteration"]][row_i].parse::<usize>().unwrap() {
            panic!(
                "Number of iterator at row {} should be {} but {}.",
                row_i,
                row_i - 1,
                (vss[col_to_i["iteration"]][row_i].parse::<usize>().unwrap())
            );
        }

        let wgt = load_wgt_func(vss, &col_to_i, row_i);
        //let wgt = load_wgt_row(vss, &col_to_i, wgt_i);

        wgts.push(wgt);
    }

    // check iters are sorted.
    assert!(vec::is_sorted(
        &wgts.iter().map(|v| v.iteration()).collect::<Vec<usize>>()
    ));

    wgts
}

fn judge_wgt_version(vss: &Vec<Vec<String>>, col_to_i: &HashMap<String, usize>) -> String {
    for row_i in 1..vss[0].len() {
        let model_name = &vss[col_to_i["model"]][row_i];
        // These are new in ver. 2
        if ["LINEARC", "LINEARCI", "FREE2", "FREE3", "FREE4"].contains(&model_name.as_str()) {
            return "ver2".to_owned();
        }
    }

    // Even if wgt is ver2 but does not include those rows, treat as "ver1"
    return "ver1".to_owned();
}

//pub fn load_wgt_row(
//    vss: &Vec<Vec<String>>,
//    col_to_i: &HashMap<String, usize>,
//    wgt_i: usize,
//) -> WgtBoost {
//    //let model_name = &vss[col_to_i["model"]][wgt_i];
//    load_wgt_row_not_safe(vss, col_to_i, wgt_i)
//    //if ["LINEARI", "FREE2", "FREE3", "FREE4"].contains(model_name) {
//    //    load_wgt_row_safe(vss, col_to_i, wgt_i)
//    //} else {
//    //    load_wgt_row_old_240116(vss, col_to_i, wgt_i)
//    //}
//}

// model name is unique for Coef
//pub fn load_wgt_row_safe()

/// ver2: safe ver.
/// model_name is unique to Coef
pub fn load_wgt_row_ver2(
    vss: &Vec<Vec<String>>,
    col_to_i: &HashMap<String, usize>,
    wgt_i: usize,
) -> WgtBoost {
    // kind
    let wgt_kind = &vss[col_to_i["kind"]][wgt_i];
    let model_name = &vss[col_to_i["model"]][wgt_i];

    if wgt_kind == "CONST" {
        // const from cov should be in COV?
        // but cannot distinguish when no cov
        // -> unnecessary

        let model = if model_name == "LINEAR" {
            let coef = vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap();
            Model::new_linear(coef)
        } else if model_name == "BINARY" {
            let coef = Coef::Binary((
                vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
            ));
            Model::new(coef, Some(f64::NAN))
        } else {
            panic!("Unknown model.");
        };

        let wgt = Wgt::new_const(model);
        let wgt_boost = WgtBoost::new_wgt_iteration(wgt, wgt_i);
        wgt_boost
    } else if wgt_kind == "COV" {
        // cov
        // only implement for linear
        let name = vss[col_to_i["var"]][wgt_i].clone();
        let coef = vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap();
        let wgt = Wgt::new_cov_name(name, Coef::Linear(coef));
        let wgt_boost = WgtBoost::new_wgt_iteration(wgt, wgt_i);
        wgt_boost
    } else if wgt_kind == "SNV" {
        // judge coef type for each row
        let boost_coef_type = load_wgt_boost_coef_type_ver2(vss, &col_to_i, wgt_i);

        let model = match boost_coef_type {
            Coef::Binary(_) => {
                let coef = Coef::new_binary_check((
                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let threshold = vss[col_to_i["threshold"]][wgt_i].parse::<f64>().unwrap();
                let model = Model::new(coef, Some(threshold));
                model
            }
            Coef::LinearConst(_) => {
                let coef = Coef::new_linearconst_check((
                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let model = Model::new_coef(coef);
                model
            }
            Coef::LinearConstInteraction(_) => {
                let coef = Coef::new_linearconstinteraction_check((
                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let model = Model::new_coef(coef);
                model
            }

            Coef::Score3(_) => {
                let coef = Coef::new_score3_check((
                    vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let model = Model::new_coef(coef);
                model
            }
            //BoostType::FreeModelMissing | BoostType::Logit => {
            Coef::Score4(_) => {
                let coef = Coef::new_score4_check((
                    vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["scorem"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let model = Model::new_coef(coef);
                model
            }
            Coef::Linear(_) | Coef::Single(_) | Coef::Score2(_) | Coef::NaN => {
                unimplemented!();
            }
        };

        let snv = SnvId::new(
            vss[col_to_i["var"]][wgt_i].clone(),
            &vss[col_to_i["chrom"]][wgt_i],
            &vss[col_to_i["pos"]][wgt_i],
            vss[col_to_i["a1"]][wgt_i].clone(),
            vss[col_to_i["a2"]][wgt_i].clone(),
        );
        let a2_frq = vss[col_to_i["a2_frq"]][wgt_i].parse::<f64>().ok();
        let maf = a2_frq.map(|x| 1.0 - x);

        // SnvWgt or SnvInteractionWgt for WgtKind
        let wgt = match boost_coef_type {
            Coef::LinearConstInteraction(_) => {
                let snv_2 = SnvId::new(
                    vss[col_to_i["var_2"]][wgt_i].clone(),
                    &vss[col_to_i["chrom_2"]][wgt_i],
                    &vss[col_to_i["pos_2"]][wgt_i],
                    vss[col_to_i["a1_2"]][wgt_i].clone(),
                    vss[col_to_i["a2_2"]][wgt_i].clone(),
                );

                //let a2_frq = vss[col_to_i["a2_frq"]][wgt_i].parse::<f64>().unwrap();
                let a2_2_frq = vss[col_to_i["a2_2_frq"]][wgt_i].parse::<f64>().ok();
                let maf_2 = a2_2_frq.map(|x| 1.0 - x);
                let snv_interaction_wgt = SnvInteractionWgt::new_score(snv, snv_2, maf, maf_2);
                Wgt::new(WgtKind::new_snv_interaction(snv_interaction_wgt), model)
            }
            _ => {
                let snv_wgt = SnvWgt::new_score(snv, maf);
                //let snv_wgt = SnvWgt::new_snv_id(snv);
                Wgt::new(WgtKind::new_snv(snv_wgt), model)
                //let wgt = Wgt::new(WgtKind::new_snv(snv_wgt), model);
            }
        };
        //let wgt = Wgt::new(WgtKind::Snv(snv, None, None), model);
        let wgt_boost = WgtBoost::new_wgt_iteration(wgt, wgt_i);
        wgt_boost
    } else {
        panic!("Unknown wgt kind: {}", vss[1][wgt_i]);
    }
}

pub fn load_wgt_boost_coef_type_ver2(
    vss: &Vec<Vec<String>>,
    col_to_i: &HashMap<String, usize>,
    wgt_i: usize,
) -> Coef {
    let model = vss[col_to_i["model"]][wgt_i].as_str();

    let coef = Coef::new_from_model_name(model);

    let coef = coef.unwrap_or_else(|| panic!("Unknown coef type. row: {}", wgt_i));

    coef

    //let coef = match model {
    //    "LINEAR" => Coef::Linear((f64::NAN)),
    //    "LINEARC" => Coef::LinearConst((f64::NAN, f64::NAN)),
    //    "LINEARCI" => Coef::LinearConstInteraction((f64::NAN, f64::NAN)),
    //    "SINGLE" => Coef::Single((f64::NAN)),
    //    "BINARY" => Coef::Binary((f64::NAN, f64::NAN)),
    //    "FREE2" => Coef::Score2((f64::NAN, f64::NAN)),
    //    "FREE3" => Coef::Score3((f64::NAN, f64::NAN, f64::NAN)),
    //    "FREE4" => Coef::Score4((f64::NAN, f64::NAN, f64::NAN, f64::NAN)),
    //    _ => panic!("Unknown coef type.: {}", model),
    //};

    // TODO: others

    // alpha, const is also for cov
    // check if all NAN for snvs

    // const, alpha
    // LogitAdd
    //return Coef::LinearConst((f64::NAN, f64::NAN));
}

/// Ver1: not safe
/// Several Coefs have same model name, so might be confusing and unsafe.
/// However, the change will affect older version of wgt and python etc.
pub fn load_wgt_row_ver1(
    vss: &Vec<Vec<String>>,
    col_to_i: &HashMap<String, usize>,
    wgt_i: usize,
) -> WgtBoost {
    // kind
    let wgt_kind = &vss[col_to_i["kind"]][wgt_i];
    let model_name = &vss[col_to_i["model"]][wgt_i];

    if wgt_kind == "CONST" {
        // const from cov should be in COV?
        // but cannot distinguish when no cov
        // -> unnecessary

        let model = if model_name == "LINEAR" {
            let coef = vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap();
            Model::new_linear(coef)
        } else if model_name == "BINARY" {
            let coef = Coef::Binary((
                vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
            ));
            Model::new(coef, Some(f64::NAN))
        } else {
            panic!("Unknown model.");
        };

        let wgt = Wgt::new_const(model);
        let wgt_boost = WgtBoost::new_wgt_iteration(wgt, wgt_i);
        wgt_boost
    } else if wgt_kind == "COV" {
        // cov
        // only implement for linear
        let name = vss[col_to_i["var"]][wgt_i].clone();
        let coef = vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap();
        let wgt = Wgt::new_cov_name(name, Coef::Linear(coef));
        let wgt_boost = WgtBoost::new_wgt_iteration(wgt, wgt_i);
        wgt_boost
    } else if wgt_kind == "SNV" {
        let snv = SnvId::new(
            vss[col_to_i["var"]][wgt_i].clone(),
            &vss[col_to_i["chrom"]][wgt_i],
            &vss[col_to_i["pos"]][wgt_i],
            vss[col_to_i["a1"]][wgt_i].clone(),
            vss[col_to_i["a2"]][wgt_i].clone(),
        );

        // judge coef type for each row
        let boost_coef_type = load_wgt_boost_coef_type_ver1(vss, &col_to_i, wgt_i);

        let model = match boost_coef_type {
            Coef::Binary(_) => {
                let coef = Coef::new_binary_check((
                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let threshold = vss[col_to_i["threshold"]][wgt_i].parse::<f64>().unwrap();
                let model = Model::new(coef, Some(threshold));
                model
            }
            Coef::LinearConst(_) => {
                let coef = Coef::new_linearconst_check((
                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let model = Model::new_coef(coef);
                model
            }
            Coef::LinearConstInteraction(_) => {
                let coef = Coef::new_linearconst_check((
                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let model = Model::new_coef(coef);
                model
            }

            Coef::Score3(_) => {
                let coef = Coef::new_score3_check((
                    vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let model = Model::new_coef(coef);
                model
            }
            //BoostType::FreeModelMissing | BoostType::Logit => {
            Coef::Score4(_) => {
                let coef = Coef::new_score4_check((
                    vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["scorem"]][wgt_i].parse::<f64>().unwrap(),
                ));
                let model = Model::new_coef(coef);
                model
            }
            Coef::Linear(_) | Coef::Single(_) | Coef::Score2(_) | Coef::NaN => {
                unimplemented!();
            }
        };

        let snv_wgt = SnvWgt::new_snv_id(snv);
        let wgt = Wgt::new(WgtKind::new_snv(snv_wgt), model);
        //let wgt = Wgt::new(WgtKind::Snv(snv, None, None), model);
        let wgt_boost = WgtBoost::new_wgt_iteration(wgt, wgt_i);
        wgt_boost
    } else {
        panic!("Unknown wgt kind: {}", vss[1][wgt_i]);
    }
}

pub fn load_wgt_boost_coef_type_ver1(
    vss: &Vec<Vec<String>>,
    col_to_i: &HashMap<String, usize>,
    wgt_i: usize,
) -> Coef {
    let model = vss[col_to_i["model"]][wgt_i].as_str();

    // check order is important
    let cols = ["score0", "score1", "score2", "scorem"];
    if model == "FREE" && non_nan_all(&cols, vss, col_to_i, wgt_i) {
        return Coef::Score4((f64::NAN, f64::NAN, f64::NAN, f64::NAN));
    }

    let cols = ["score0", "score1", "score2"];
    if model == "FREE" && non_nan_all(&cols, vss, col_to_i, wgt_i) {
        return Coef::Score3((f64::NAN, f64::NAN, f64::NAN));
    }

    let cols = ["score0", "score1", "threshold"];
    if model == "BINARY" && non_nan_all(&cols, vss, col_to_i, wgt_i) {
        return Coef::Score2((f64::NAN, f64::NAN));
    }

    let cols = ["const", "alpha", "threshold"];
    if model == "BINARY" && non_nan_all(&cols, vss, col_to_i, wgt_i) {
        return Coef::Binary((f64::NAN, f64::NAN));
    }

    // TODO: check threshold is NaN
    let cols = ["const", "alpha"];
    if model == "LINEAR" && non_nan_all(&cols, vss, col_to_i, wgt_i) {
        return Coef::LinearConst((f64::NAN, f64::NAN));
    }

    // const, threshold is NaN
    let cols = ["alpha"];
    if model == "LINEAR" && non_nan_all(&cols, vss, col_to_i, wgt_i) {
        return Coef::Linear(f64::NAN);
    }

    panic!("Unknown coef type. row: {}", wgt_i);

    // TODO: others

    // alpha, const is also for cov
    // check if all NAN for snvs

    // const, alpha
    // LogitAdd
    //return Coef::LinearConst((f64::NAN, f64::NAN));
}

/* // TODO: binary of cov and linear of snv
pub fn load_wgts(fin_wgt: &str) -> Vec<WgtBoost> {
    let fwgt = get_fname_wgt(fin_wgt);

    let mut wgts: Vec<WgtBoost> = Vec::new();

    let vss: Vec<Vec<String>> =
        text::load_table(&fwgt, true).expect(&format!("Cannot load wgt from {}", fwgt));

    for wgt_i in 0..vss[0].len() {
        assert_eq!(
            wgt_i,
            (vss[0][wgt_i].parse::<usize>().unwrap()),
            "Number of iterator in wgt is wrong at iteration {} but found {}.",
            wgt_i,
            (vss[0][wgt_i].parse::<usize>().unwrap()),
        );

        // kind
        let wgt_kind = &vss[1][wgt_i];
        let model_name = &vss[3][wgt_i];
        let wgt = if wgt_kind == "CONST" {
            //const
            let threshold = vss[4][wgt_i].parse::<f64>().unwrap();

            let name = vss[2][wgt_i].clone();
            let const_var = Var::construct_var(CovKind::Const, name);
            // later
            //const_var.set_vals_const(n);
            let const_wgt = CovWgt::construct(const_var);

            let model = if model_name == "LINEAR" {
                let coef = vss[5][wgt_i].parse::<f64>().unwrap();
                Model::new_linear(coef)
            } else if model_name == "BINARY" {
                let coef = Coef::Binary((
                    vss[5][wgt_i].parse::<f64>().unwrap(),
                    vss[6][wgt_i].parse::<f64>().unwrap(),
                ));
                Model::new(coef, Some(threshold))
                //Model::Binary(BinaryModel::new(coef, threshold))
            } else {
                panic!("Unknown model.");
            };
            /*
            let model = if threshold.is_nan() {
                let coef = vss[5][wgt_i].parse::<f64>().unwrap();
                Model::Linear(LinearModel::construct(coef))
            } else {
                Model::Binary(BinaryModel::construct_threshold(threshold))
            };
            */

            let wgt = Wgt::construct_const(model);
            //let wgt = Wgt::construct_wgt(WgtKind::Cov(const_wgt), model);
            //let wgt = Wgt::construct_wgt(WgtKind::Cov(const_wgt), model, wgt_i);
            let wgt_boost = WgtBoost::construct_wgt_iteration(wgt, wgt_i);
            wgt_boost
        } else if wgt_kind == "COV" {
            // cov
            // only implement for linear
            let name = vss[2][wgt_i].clone();
            let coef = vss[5][wgt_i].parse::<f64>().unwrap();
            let wgt = Wgt::construct_cov_name(name, Coef::Linear(coef));
            /*
            let cov_var = Var::construct_var(CovKind::Cov, name.clone());
            let cov_wgt = CovWgt::construct(cov_var);
            let coef = vss[5][wgt_i].parse::<f64>().unwrap();
            let model = Model::new_linear(coef);
            let wgt = Wgt::construct_wgt(WgtKind::Cov(cov_wgt), model);
             */
            let wgt_boost = WgtBoost::construct_wgt_iteration(wgt, wgt_i);
            wgt_boost
        } else if wgt_kind == "SNV" {
            let snv = SnvId::construct_snv_index_string(
                vss[2][wgt_i].clone(),
                &vss[7][wgt_i],
                &vss[8][wgt_i],
                vss[9][wgt_i].clone(),
                vss[10][wgt_i].clone(),
            );
            //let snv_wgt = SnvWgt::construct_score(snv);
            let coef = Coef::Binary((
                vss[5][wgt_i].parse::<f64>().unwrap(),
                vss[6][wgt_i].parse::<f64>().unwrap(),
            ));
            let threshold = vss[4][wgt_i].parse::<f64>().unwrap();
            let model = Model::new(coef, Some(threshold));
            //let model = Model::Binary(BinaryModel::new(coef, threshold));
            let wgt = Wgt::construct_wgt(WgtKind::Snv(snv, None, None), model);
            let wgt_boost = WgtBoost::construct_wgt_iteration(wgt, wgt_i);
            wgt_boost
        } else {
            panic!("Unknown wgt kind: {}", vss[1][wgt_i]);
        };

        wgts.push(wgt);
    }

    // check iters are sorted.
    assert!(vec::is_sorted(
        &wgts.iter().map(|v| v.iteration()).collect::<Vec<usize>>()
    ));
    /*
    // is_sorted() is nightly
    assert!(wgts
        .iter()
        .map(|v| v.get_iter())
        .collect::<Vec<usize>>()
        .is_sorted());
         */

    wgts
} */

/*
pub fn set_covs(wgts: &mut [WgtBoost], covs_in: Option<&[Var]>, n: usize) {
    if let Some(covs) = covs_in {
        for wgt in wgts.iter_mut() {
            if let WgtKind::Cov(ref mut cov_wgt) = wgt.wgt_mut().kind_mut() {
                let cov_in_wgt = cov_wgt.get_var_mut();
                if cov_in_wgt.name() == "const" {
                    cov_in_wgt.set_vals_const(n);
                } else {
                    for cov in covs.iter() {
                        if cov.name() == cov_in_wgt.name() {
                            cov_in_wgt.set_vals(cov.vals().to_vec());
                            break;
                        }
                    }
                }
            }
        }
    }
}
 */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_non_nan_all() {
        let cols = vec!["score0", "score1"];
        let col_to_i = vec![("score0".to_string(), 0), ("score1".to_string(), 1)];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let valss = vec![
            "0.0,1.0"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
            "2.0,3.0"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
        ];

        let wgt_i = 1;

        assert_eq!(non_nan_all(&cols, &valss, &col_to_i, wgt_i), true);
    }

    /// false if NAN exists
    #[test]
    fn test_non_nan_all_false() {
        let cols = vec!["score0", "score1"];
        let col_to_i = vec![("score0".to_string(), 0), ("score1".to_string(), 1)];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let valss = vec![
            "0.0,NaN"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
            "2.0,NaN"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
        ];

        let wgt_i = 1;

        assert_eq!(non_nan_all(&cols, &valss, &col_to_i, wgt_i), false);
    }

    /// false when one of cols is not in col_to_i
    #[test]
    fn test_non_nan_all_false2() {
        let cols = vec!["score0", "score1", "score2", "scorem"];
        let col_to_i = vec![
            ("score0".to_string(), 0),
            ("score1".to_string(), 1),
            ("score2".to_string(), 2),
        ];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let valss = vec![
            "0.0,1.0"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
            "2.0,3.0"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
            "4.0,5.0"
                .split(",")
                .map(|x| x.to_string())
                .collect::<Vec<String>>(),
        ];

        let wgt_i = 1;

        assert_eq!(non_nan_all(&cols, &valss, &col_to_i, wgt_i), false);
    }

    #[test]
    fn test_load_wgt_boost_coef_type_score4() {
        let col_to_i = vec![
            ("model".to_string(), 0),
            ("score0".to_string(), 1),
            ("score1".to_string(), 2),
            ("score2".to_string(), 3),
            ("scorem".to_string(), 4),
        ];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let vss = vec![
            vec!["FREE".to_string()],
            vec!["0.0".to_string()],
            vec!["1.0".to_string()],
            vec!["2.0".to_string()],
            vec!["3.0".to_string()],
        ];
        let wgt_i = 0;
        let coef = load_wgt_boost_coef_type_ver1(&vss, &col_to_i, wgt_i);
        let exp = Coef::Score4((f64::NAN, f64::NAN, f64::NAN, f64::NAN));

        assert!(coef.match_type(exp));
    }

    #[test]
    fn test_load_wgt_boost_coef_type_score3() {
        let col_to_i = vec![
            ("model".to_string(), 0),
            ("score0".to_string(), 1),
            ("score1".to_string(), 2),
            ("score2".to_string(), 3),
        ];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let vss = vec![
            vec!["FREE".to_string()],
            vec!["0.0".to_string()],
            vec!["1.0".to_string()],
            vec!["2.0".to_string()],
        ];
        let wgt_i = 0;
        let coef = load_wgt_boost_coef_type_ver1(&vss, &col_to_i, wgt_i);
        let exp = Coef::Score3((f64::NAN, f64::NAN, f64::NAN));

        assert!(coef.match_type(exp));
    }

    #[test]
    fn test_load_wgt_boost_coef_type_score2() {
        let col_to_i = vec![
            ("model".to_string(), 0),
            ("score0".to_string(), 1),
            ("score1".to_string(), 2),
            ("threshold".to_string(), 3),
        ];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let vss = vec![
            vec!["BINARY".to_string()],
            vec!["0.0".to_string()],
            vec!["1.0".to_string()],
            vec!["0.5".to_string()],
        ];
        let wgt_i = 0;
        let coef = load_wgt_boost_coef_type_ver1(&vss, &col_to_i, wgt_i);
        let exp = Coef::Score2((f64::NAN, f64::NAN));

        assert!(coef.match_type(exp));
    }

    #[test]
    fn test_load_wgt_boost_coef_type_binary() {
        let col_to_i = vec![
            ("model".to_string(), 0),
            ("const".to_string(), 1),
            ("alpha".to_string(), 2),
            ("threshold".to_string(), 3),
        ];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let vss = vec![
            vec!["BINARY".to_string()],
            vec!["0.0".to_string()],
            vec!["1.0".to_string()],
            vec!["0.5".to_string()],
        ];
        let wgt_i = 0;
        let coef = load_wgt_boost_coef_type_ver1(&vss, &col_to_i, wgt_i);
        let exp = Coef::Binary((f64::NAN, f64::NAN));

        assert!(coef.match_type(exp));
    }

    #[test]
    fn test_load_wgt_boost_coef_type_linearconst() {
        let col_to_i = vec![
            ("model".to_string(), 0),
            ("const".to_string(), 1),
            ("alpha".to_string(), 2),
        ];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let vss = vec![
            vec!["LINEAR".to_string()],
            vec!["0.0".to_string()],
            vec!["1.0".to_string()],
        ];
        let wgt_i = 0;
        let coef = load_wgt_boost_coef_type_ver1(&vss, &col_to_i, wgt_i);
        let exp = Coef::LinearConst((f64::NAN, f64::NAN));

        assert!(coef.match_type(exp));
    }

    #[test]
    fn test_load_wgt_boost_coef_type_linear() {
        let col_to_i = vec![("model".to_string(), 0), ("alpha".to_string(), 1)];
        let col_to_i: HashMap<String, usize> = col_to_i.into_iter().collect();
        let vss = vec![vec!["LINEAR".to_string()], vec!["1.0".to_string()]];
        let wgt_i = 0;
        let coef = load_wgt_boost_coef_type_ver1(&vss, &col_to_i, wgt_i);
        let exp = Coef::Linear(f64::NAN);

        assert!(coef.match_type(exp));
    }
}
