// use AsRef<Path>
// https://www.reddit.com/r/rust/comments/7mu7q1/is_working_with_paths_always_this_painful/
// should use &Path not &Path

use super::WgtBoost;
use super::{Coef, Model, Wgt, WgtKind};
use crate::BoostType;
use genetics::{text, vec, SnvId};
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn get_file_wgt(dout: &Path, learning_rate: &Option<f64>) -> PathBuf {
    let p = get_dname_para(dout, learning_rate);
    get_fname_wgt(&p)
}

// when increase para, add them
pub fn get_dname_para(dout: &Path, learning_rate: &Option<f64>) -> PathBuf {
    let mut d = dout.to_owned();
    let mut dpara = String::from("para");
    if let Some(lr) = learning_rate {
        dpara += &(String::from(".lr") + &lr.to_string());
    }
    d.push(dpara);
    d
}

pub fn get_dname_wgt(dout: &Path) -> PathBuf {
    let d = dout.to_owned();
    //let dwgt = "wgt/";
    //d.push(dwgt);
    d
}

pub fn get_fname_wgt(dout_para: &Path) -> PathBuf {
    let mut fwgt = get_dname_wgt(dout_para);

    fwgt.push("boosting.wgt");
    //fwgt.set_file_name("boosting.wgt");
    fwgt
    //fout.to_owned() + ".wgt"
}

/* pub fn get_fname_wgt(fout: &str) -> String {
    fout.to_owned() + ".wgt"
} */

pub fn create_dir(dout: &Path) {
    fs::create_dir_all(&dout).unwrap();
}

pub fn get_dname_loss(dout: &Path) -> PathBuf {
    let mut f = dout.to_owned();
    let dscore = "loss/";
    f.push(dscore);
    f
}

pub fn get_fname_loss(dout: &Path, ti: usize) -> PathBuf {
    let mut f = get_dname_loss(dout);
    let fname = "iter".to_string() + &ti.to_string() + ".loss";
    f.push(fname);
    f
    //fout.to_owned() + ".iter" + &ti.to_string() + ".loss"
}

/* pub fn get_fname_ws(fout: &str) -> String {
    fout.to_owned() + ".ws"
} */

pub fn check_file_wgt_not_exist(dout: &Path, learning_rate: &Option<f64>) {
    let fwgt = get_file_wgt(dout, learning_rate);
    let exist_fwgt = text::exist_file(&fwgt);
    if exist_fwgt {
        panic!(
            "Weight file already exists: {:?}. Delete it or use --resume option.",
            &fwgt
        );
    }
}

pub fn check_file_wgt_exist(fwgt: &Path) {
    let exist_fwgt = text::exist_file(&fwgt);
    if !exist_fwgt {
        panic!("fwgt does not exist: {:?}.", &fwgt);
    }
}

#[allow(dead_code)]
pub fn check_file_wgt_exist_dir(dout: &Path) {
    let fwgt = get_fname_wgt(dout);
    check_file_wgt_exist(&fwgt)
}

pub fn bufwriter_floss(dout: &Path, ti: usize) -> BufWriter<File> {
    let dloss = get_dname_loss(dout);
    create_dir(&dloss);
    let floss = get_fname_loss(dout, ti);
    //if text::exist_file(&floss) {
    //    panic!("File already exists: {}.", fout);
    //}

    let file = match File::create(&floss) {
        Ok(file) => file,
        Err(_) => panic!(
            "Cannot create file, possibly directory does not exist: {:?}",
            &floss
        ),
    };

    BufWriter::new(file)
    /*
    let mut writer = BufWriter::new(file);
    let str =
        "iteration\tkind\tvar\tmodel\tthreshold\talpha\tconst\tchrom\tpos\ta1\ta2\n".to_owned();
    writer.write(str.as_bytes()).unwrap();
     */
}

pub fn is_exist_wgt(dout: &Path) -> bool {
    let fwgt = get_fname_wgt(dout);
    text::exist_file(&fwgt)
}

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

pub fn bufwriter_fwgt_append(dout: &Path) -> BufWriter<File> {
    // create dwgt
    //let dwgt = get_dname_wgt(dout);
    //create_dir(&dwgt);

    let fwgt = get_fname_wgt(dout);
    // allow exist for resume
    //if text::exist_file(&fwgt) {
    //    panic!("File already exists: {}.", fwgt);
    //}

    let file = match OpenOptions::new()
        .read(true) // unnecessary
        .append(true)
        .create(true)
        .open(&fwgt)
    {
        Ok(file) => file,
        Err(_) => panic!("file does not exist: {:?}", &fwgt),
    };

    // TO append
    // https://masahiko-ofgp-ja-notebook.blogspot.com/2019/05/rust-lang-file-append-mode.html
    //let file = OpenOptions::new().append(true).open(&fwgt).unwrap();

    BufWriter::new(file)
}

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
        BoostType::FreeModelMissing => vec!["score0", "score1", "score2", "scorem"],
    };
    let mut columns = vec![
        "iteration",
        "kind",
        "var",
        "model", //"threshold",
        "eps",
        "alpha",
        "const",
    ];
    let mut cols_snv = vec!["chrom", "pos", "a1", "a2"];
    columns.append(&mut cols_boost);
    columns.append(&mut cols_snv);
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

    println!("var_name_to_coef {:?}", var_name_to_coef);

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

    println!("wgts_cov: {:?}", wgts_cov);

    wgts_cov
}
     */

#[allow(dead_code)]
fn load_num_iteration_dir(dout: &Path) -> usize {
    let fwgt = get_fname_wgt(dout);
    let iter_col = 0;
    let iters = text::load_table_col(&fwgt, iter_col, true).unwrap();
    let iters: Vec<usize> = iters.iter().map(|s| s.parse::<usize>().unwrap()).collect();
    let num_iteration = iters.iter().max().unwrap() + 1;
    num_iteration
}

fn load_num_iteration(fwgt: &Path) -> usize {
    //let fwgt = get_fname_wgt(dout);
    let iter_col = 0;
    let iters = text::load_table_col(&fwgt, iter_col, true).unwrap();
    let iters: Vec<usize> = iters.iter().map(|s| s.parse::<usize>().unwrap()).collect();
    let num_iteration = iters.iter().max().unwrap() + 1;
    num_iteration
}

#[allow(dead_code)]
pub fn valid_iterations_dir(iterations_in: &[usize], dout: &Path) -> Vec<usize> {
    //let fwgt = get_fname_wgt(dout);

    let mut iterations = iterations_in.to_vec();
    iterations.sort();

    let num_iteration = load_num_iteration_dir(dout);

    let iterations = iterations
        .into_iter()
        .filter(|&v| v <= num_iteration)
        .collect();

    iterations
}

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

pub fn load_wgts(dout: &Path, boost_type: BoostType) -> Vec<WgtBoost> {
    let fwgt = get_fname_wgt(dout);

    load_wgts_file(&fwgt, boost_type)
}

// TODO: binary of cov and linear of snv
pub fn load_wgts_file(fwgt: &Path, boost_type: BoostType) -> Vec<WgtBoost> {
    //let fwgt = get_fname_wgt(dout);

    let mut wgts: Vec<WgtBoost> = Vec::new();

    // also load header
    let vss: Vec<Vec<String>> =
        text::load_table(fwgt, false).expect(&format!("Cannot load wgt from {:?}", fwgt));

    let col_n = vss.len();
    let mut columns: Vec<String> = Vec::new();
    for col_i in 0..col_n {
        columns.push(vss[col_i][0].clone());
    }

    let mut col_to_i: HashMap<String, usize> = HashMap::new();
    for col_i in 0..col_n {
        col_to_i.insert(columns[col_i].clone(), col_i);
    }
    println!("col_to_i {:?}", col_to_i);

    for wgt_i in 1..vss[0].len() {
        assert_eq!(
            wgt_i - 1,
            (vss[col_to_i["iteration"]][wgt_i].parse::<usize>().unwrap()),
            "Number of iterator in wgt is wrong at iteration {} but found {}.",
            wgt_i - 1,
            (vss[col_to_i["iteration"]][wgt_i].parse::<usize>().unwrap()),
        );

        // kind
        let wgt_kind = &vss[col_to_i["kind"]][wgt_i];
        let model_name = &vss[col_to_i["model"]][wgt_i];
        let wgt = if wgt_kind == "CONST" {
            //const
            //let threshold = vss[4][wgt_i].parse::<f64>().unwrap();
            // should not use?
            let threshold = 0.0;

            //let name = vss[col_to_i["var"]][wgt_i].clone();
            //let const_var = Var::construct_var(CovKind::Const, name);
            // later
            //const_var.set_vals_const(n);
            // FIXME: how to use CovWgt?
            //let const_wgt = CovWgt::construct(const_var);

            let model = if model_name == "LINEAR" {
                let coef = vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap();
                Model::new_linear(coef)
            } else if model_name == "BINARY" {
                let coef = Coef::Binary((
                    vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                    vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
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
            let name = vss[col_to_i["var"]][wgt_i].clone();
            let coef = vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap();
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
                vss[col_to_i["var"]][wgt_i].clone(),
                &vss[col_to_i["chrom"]][wgt_i],
                &vss[col_to_i["pos"]][wgt_i],
                vss[col_to_i["a1"]][wgt_i].clone(),
                vss[col_to_i["a2"]][wgt_i].clone(),
            );
            //let snv_wgt = SnvWgt::construct_score(snv);

            let model = match boost_type {
                BoostType::ConstAda => {
                    let coef = Coef::Binary((
                        vss[col_to_i["const"]][wgt_i].parse::<f64>().unwrap(),
                        vss[col_to_i["alpha"]][wgt_i].parse::<f64>().unwrap(),
                    ));
                    let threshold = vss[col_to_i["threshold"]][wgt_i].parse::<f64>().unwrap();
                    let model = Model::new(coef, Some(threshold));
                    model
                }
                BoostType::FreeModelMissing => {
                    let coef = Coef::Score4((
                        vss[col_to_i["score0"]][wgt_i].parse::<f64>().unwrap(),
                        vss[col_to_i["score1"]][wgt_i].parse::<f64>().unwrap(),
                        vss[col_to_i["score2"]][wgt_i].parse::<f64>().unwrap(),
                        vss[col_to_i["scorem"]][wgt_i].parse::<f64>().unwrap(),
                    ));
                    let model = Model::new_coef(coef);
                    model
                }
                BoostType::Ada => {
                    unimplemented!();
                }
            };

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
    fn test_get_dname_para() {
        let dout = PathBuf::from("./abc");
        let lr = Some(0.1);
        let dout_para = get_dname_para(&dout, &lr);
        assert_eq!(dout_para, PathBuf::from("./abc/para.lr0.1/"));
    }

    #[test]
    fn test_get_dname_para_lrnone() {
        let dout = PathBuf::from("./abc");
        let lr = None;
        let dout_para = get_dname_para(&dout, &lr);
        assert_eq!(dout_para, PathBuf::from("./abc/para/"));
    }
}
