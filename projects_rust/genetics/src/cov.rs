// TODO?: wrap up all struct to one submodule?
// -> when create struct Sample,
// -> create mod structures/snv.rs and structures/sample.rs ?
// -> unnecessary now
//
// We should not use CovKind::Var since only Covkind::Cov and Covkind::Const are enough

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

//// DEPRECATED: Use Covs instead
//// vals should be here?
//// -> use VarIndex if vals are not necessary like in wgt
//#[derive(Debug, Clone)]
//pub struct Cov {
//    // TODO: substitute to VarIndex
//    kind: CovKind,
//    name: String,
//    //n: usize, // unnecessary if vals.len()==n
//    vals: Vec<f64>, // values of samples
//                    /*
//                    vals_sorted: Vec<f64>, // for values Binary
//                    index_sorted: Vec<usize>,
//                     */
//}
//
//impl Cov {
//    pub fn construct_cov(kind: CovKind, name: String) -> Self {
//        Self {
//            kind,
//            name,
//            vals: Vec::new(),
//        }
//    }
//
//    // use Default instead
//    //pub fn construct_empty() -> Var {
//    //    Var {
//    //        kind: VarKind::Const,
//    //        name: "".to_owned(),
//    //        vals: Vec::new(),
//    //    }
//    //}
//
//    pub fn set_vals(&mut self, vals: Vec<f64>) {
//        self.vals = vals;
//    }
//
//    pub fn set_vals_const(&mut self, n: usize) {
//        let v = vec![1.0; n];
//
//        self.set_vals(v);
//        //self.vals = v.clone();
//    }
//
//    pub fn kind(&self) -> CovKind {
//        self.kind
//    }
//
//    pub fn name(&self) -> &str {
//        &self.name
//    }
//
//    pub fn vals(&self) -> &[f64] {
//        &self.vals
//    }
//    pub fn vals_consume(self) -> Vec<f64> {
//        self.vals
//    }
//}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_covkind_default() {
        assert_eq!(CovKind::default(), CovKind::Const);
    }

    #[test]
    fn test_covkind_from_str() {
        assert_eq!(CovKind::from_str("const").unwrap(), CovKind::Const);
        assert_eq!(CovKind::from_str("cov").unwrap(), CovKind::Cov);
    }
}
