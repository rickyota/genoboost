// TODO?: wrap up all struct to one submodule?
// -> when create struct Sample,
// -> create mod structures/snv.rs and structures/sample.rs ?
// -> unnecessary now
// We should not use var since only Cov and const are enough
// rename to cov?

use super::{samples, text};
use std::collections::HashMap;
use std::path::Path;
use std::str::FromStr;

// rename to CovKind
#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug)]
pub enum CovKind {
    Const, //const; necessary since value of the var is not provided from input
    Cov,   // assumed to be excluded, same as variates in fin_cov
}

//#[derive(Debug, Clone, Copy)]
//#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Copy, Hash, Debug)]
//pub enum VarKind {
//    Const, //const
//    Var,   // same layer as SNV
//    Cov,   // assumed to be excluded, same as variates in fin_cov
//}

impl Default for CovKind {
    fn default() -> Self {
        CovKind::Const
    }
}

impl FromStr for CovKind {
    // FromStr shoudl be `Err`
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            "const" => Ok(CovKind::Const),
            "cov" => Ok(CovKind::Cov),
            _ => Err(format!(
                "VarKind string should be one of const, cov not {}",
                str
            )),
        }

        /*
        if str == "const" {
            return Ok(CovKind::Const);
        //} else if str == "var" {
        //    return Ok(VarKind::Var);
        } else if str == "cov" {
            return Ok(CovKind::Cov);
        } else {
            return Err(format!(
                "VarKind string should be one of const, cov not {}",
                // "VarKind string should be one of const, var, cov not {}",
                str
            ));
        }
         */
    }
}

impl std::fmt::Display for CovKind {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let var_kind_str = match self {
            CovKind::Const => "const",
            //VarKind::Var => "var",
            CovKind::Cov => "cov",
        };
        write!(f, "{}", var_kind_str)
    }
}

// vals should be here?
// -> use VarIndex if vals are not necessary like in wgt
#[derive(Debug, Clone)]
pub struct Var {
    // TODO: substitute to VarIndex
    kind: CovKind,
    name: String,
    //n: usize, // unnecessary if vals.len()==n
    vals: Vec<f64>, // values of samples
                    /*
                    vals_sorted: Vec<f64>, // for values Binary
                    index_sorted: Vec<usize>,
                     */
}

impl Var {
    /// use this func to call Snv
    /// Snv::construct_snv()

    // TODO: change name: String -> &str
    pub fn construct_var(kind: CovKind, name: String) -> Var {
        Var {
            kind,
            name,
            vals: Vec::new(),
            /*
            vals_sorted: Vec::new(),
            index_sorted: Vec::new(),
             */
        }
    }

    /*
    // use Default instead
    pub fn construct_empty() -> Var {
        Var {
            kind: VarKind::Const,
            name: "".to_owned(),
            vals: Vec::new(),
        }
    }
    */

    //fn sort_vals_index(vals: &[f64]) -> (Vec<f64>, Vec<usize>) {}

    pub fn set_vals(&mut self, vals: Vec<f64>) {
        self.vals = vals;

        /*
        let (v_sorted, v_index) = vec::sort_and_argsort(&self.vals);
        self.vals_sorted = v_sorted;
        self.index_sorted = v_index;
         */
    }

    // how to combine to above?
    pub fn set_vals_const(&mut self, n: usize) {
        let v = vec![1.0; n];
        //let mut v = Vec::with_capacity(n);
        //for _ in 0..n {
        //    v.push(1.0);
        //}

        self.vals = v.clone();
        /*
        self.vals_sorted = v;

        let mut v_index = Vec::with_capacity(n);
        for v_i in 0..n {
            v_index.push(v_i)
        }
        self.index_sorted = v_index;
        */

        //let (v_sorted, v_index) = vec::sort_and_argsort(&self.vals);
        //self.vals_sorted = v_sorted;
        //self.index_sorted = v_index;
    }

    /*
    fn make_sida(&self) -> String {
        self.chrom.to_string() + ":" + &self.pos.to_string() + ":" + &self.a1 + ":" + &self.a2
    }
    */

    /*
    fn set_sida(&mut self) {
        self.sida =
            self.chrom.to_string() + ":" + &self.pos.to_string() + ":" + &self.a1 + ":" + &self.a2;
    }
    */

    /*
    /// if implement as String, just write `rs` is fine
    pub fn construct_snv_string(rs: &str, chrom: &str, pos: &str, a1: &str, a2: &str) -> Snv {
        Snv {
            rs: rs.to_owned(),
            chrom: chrom.to_owned().parse::<usize>().unwrap(),
            pos: pos.to_owned().parse::<usize>().unwrap(),
            a1: a1.to_owned(),
            a2: a2.to_owned(),
        }
    }
    */

    pub fn kind(&self) -> CovKind {
        self.kind
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn vals(&self) -> &[f64] {
        &self.vals
    }
    pub fn vals_consume(self) -> Vec<f64> {
        self.vals
    }
}

pub fn vals_align_id(
    vals_str: &[String],
    ids: (&[String], &[String]),
    sample_id_to_n: &HashMap<String, usize>,
    n: usize,
) -> Vec<f64> {
    let vals_in: Vec<f64> = vals_str.iter().map(|x| x.parse::<f64>().unwrap()).collect();

    let mut vals = vec![f64::NAN; n];
    //let mut vals: Vec<f64> = Vec::with_capacity(n);
    //vec::push(&mut vals, f64::NAN, n);

    for (n_in_i, val) in vals_in.iter().enumerate() {
        let sample_id = samples::sample_id(ids.0[n_in_i].clone(), &ids.1[n_in_i]);

        if let Some(ni) = sample_id_to_n.get(&sample_id) {
            vals[*ni] = *val;
        }
        // else, sample id in .cov is not used.

        /*
        match sample_id_to_n.get(&sample_id) {
            Some(ni) => {
                vals[*ni] = *val;
            }
            // FIXME: sample_id conly contains use_samples
            None => {
                panic!("Sample id in fin_var did not find in plink: {}.", sample_id);
            }
        }
        */
    }

    // assert all samples in using samples have non-nan value
    if vals.iter().any(|v| (*v).is_nan()) {
        // find one of nan index and panic!
        for n_in_i in 0..vals_in.len() {
            let sample_id = samples::sample_id(ids.0[n_in_i].clone(), &ids.1[n_in_i]);

            if let Some(ni) = sample_id_to_n.get(&sample_id) {
                if vals[*ni].is_nan() {
                    panic!("Var of some sample is not found in fin_cov: {}.", sample_id);
                }
            }
        }
    }

    vals
}

/// Add const var here
/// Assume that fin_var contains #samples the same as fin_bed not fin_samples.
/// -> both are fine
/// Assume that fin_var uses "fid" and "iid" to identify individual
pub fn load_vars(
    fin_var: Option<&Path>,
    n_in: usize,
    sample_id_to_n: &HashMap<String, usize>,
    //use_samples: &[bool],
    kind: CovKind,
) -> Option<Vec<Var>> {
    if fin_var.is_none() {
        return None;
        /*
        let vars: Vec<Var> = Vec::new();
        return vars;
         */
    }

    let cols_name = text::load_table_header(fin_var.unwrap()).unwrap();

    let valss = text::load_table(fin_var.unwrap(), true).unwrap();

    assert_eq!(valss[0].len(), n_in);

    let n = sample_id_to_n.keys().len();

    let mut vars: Vec<Var> = Vec::with_capacity(n);

    //let p = cols_name.len() - 2;

    let col_fid = 0;
    let col_iid = 1;

    for pi in 0..cols_name.len() {
        if (pi == col_fid) || (pi == col_iid) {
            continue;
        }

        let mut var = Var::construct_var(kind, cols_name[pi].clone());

        let vals = vals_align_id(
            &valss[pi],
            (&valss[col_fid], &valss[col_iid]),
            sample_id_to_n,
            n,
        );

        var.set_vals(vals);

        vars.push(var);
    }

    /*
    for n_in_i in 0..n_in {
        let sample_id = sample::sample_id(valss[0][n_in_i].clone(), &valss[1][n_in_i]);

        match sample_id_to_n.get(&sample_id) {
            Some(ni) => {}
            None => {
                panic!("Sample id in fin_var did not find in plink: {}.", sample_id);
            }
        }
    }
    */

    Some(vars)
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    use crate::plink;

    fn setup_test() -> (Option<PathBuf>, usize, HashMap<String, usize>) {
        let fin = PathBuf::from("../../test/data/toy1/genot");
        let fin_var = Some(PathBuf::from("../../test/data/toy1/genot.cov"));
        let fin_sample = None;

        let n_in: usize = plink::compute_num_sample(&fin).unwrap();
        println!("n_in: {}", n_in);

        let (n, use_samples) = plink::make_use_samples(fin_sample, &fin, n_in);
        println!("n: {}", n);

        let sample_id_to_n = samples::create_sample_id_to_n(&fin, &use_samples);

        (fin_var, n_in, sample_id_to_n)
    }

    #[test]
    fn test_load_vars() {
        let (fin_var, n_in, sample_id_to_n) = setup_test();

        let vars: Vec<Var> =
            load_vars(fin_var.as_deref(), n_in, &sample_id_to_n, CovKind::Cov).unwrap();
        //load_vars(fin_var.as_deref(), n_in, &sample_id_to_n, CovKind::Cov).unwrap();
        println!("vars {:?}", vars);
    }
}
