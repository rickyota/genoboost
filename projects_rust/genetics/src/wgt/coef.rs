use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Coef {
    // alpha or const
    Linear(f64),
    //Linear((f64,)),
    //  (const,alpha)
    LinearConst((f64, f64)),
    //  (const,alpha)
    LinearConstInteraction((f64, f64)),
    // alpha
    // for AdaBoost
    // use Binary((alpha, 0.0))
    Single(f64),
    // (const, alpha); under genetic model
    // x (alpha, const)
    // if recessive model, -alpha+const for MA=0,1; +alpha+const for MA=2
    // if dominant model, -alpha+const for MA=0; +alpha+const for MA=1,2
    Binary((f64, f64)),
    // (score0, score1); under genetic model
    // if recessive model, score0 for MA=0,1; score1 for MA=2
    // if dominant model, score0 for MA=0; score1 for MA=1,2
    Score2((f64, f64)),
    // (MA0, MA1, MA2)
    Score3((f64, f64, f64)),
    // (MA0, MA1, MA2, missing)
    Score4((f64, f64, f64, f64)),
    // unset
    NaN,
}

impl Coef {
    pub fn linearconst_f64(self) -> (f64, f64) {
        match self {
            Coef::LinearConst(ss) => ss,
            _ => panic!("not LinearConst"),
        }
    }

    pub fn linearconstinteraction_f64(self) -> (f64, f64) {
        match self {
            Coef::LinearConstInteraction(ss) => ss,
            _ => panic!("not LinearConstInteraction"),
        }
    }

    pub fn score4_f64(self) -> (f64, f64, f64, f64) {
        match self {
            Coef::Score4(ss) => ss,
            _ => panic!("not score4"),
        }
    }
    pub fn score3_f64(self) -> (f64, f64, f64) {
        match self {
            Coef::Score3(ss) => ss,
            _ => panic!("not score3"),
        }
    }

    pub fn new_score3(scores: (f64, f64, f64)) -> Self {
        Coef::Score3(scores)
    }

    pub fn new_score4(scores: (f64, f64, f64, f64)) -> Self {
        Coef::Score4(scores)
    }

    pub fn new_binary_check(scores: (f64, f64)) -> Self {
        let c = Coef::Binary(scores);
        c.check();
        c
    }

    pub fn new_linearconst_check(scores: (f64, f64)) -> Self {
        let c = Coef::LinearConst(scores);
        c.check();
        c
    }

    pub fn new_linearconstinteraction_check(scores: (f64, f64)) -> Self {
        let c = Coef::LinearConstInteraction(scores);
        c.check();
        c
    }

    pub fn new_score3_check(scores: (f64, f64, f64)) -> Self {
        let c = Coef::Score3(scores);
        c.check();
        c
    }

    pub fn new_score4_check(scores: (f64, f64, f64, f64)) -> Self {
        let c = Coef::Score4(scores);
        c.check();
        c
    }

    pub fn check(&self) {
        // check all f64 is not NaN
        let is_nan = match self {
            Coef::Linear(x) => x.is_nan(),
            Coef::LinearConst((x1, x2)) => x1.is_nan() || x2.is_nan(),
            Coef::LinearConstInteraction((x1, x2)) => x1.is_nan() || x2.is_nan(),
            Coef::Single(x) => x.is_nan(),
            Coef::Binary((x1, x2)) => x1.is_nan() || x2.is_nan(),
            Coef::Score2((x1, x2)) => x1.is_nan() || x2.is_nan(),
            Coef::Score3((x1, x2, x3)) => x1.is_nan() || x2.is_nan() || x3.is_nan(),
            Coef::Score4((x1, x2, x3, x4)) => {
                x1.is_nan() || x2.is_nan() || x3.is_nan() || x4.is_nan()
            }
            Coef::NaN => false,
        };

        if is_nan {
            panic!("coef is NaN: {:?}", self);
        }
    }

    // to compare type only
    pub fn match_type(self, other: Self) -> bool {
        match (self, other) {
            (Coef::Linear(_), Coef::Linear(_)) => true,
            (Coef::LinearConst(_), Coef::LinearConst(_)) => true,
            (Coef::LinearConstInteraction(_), Coef::LinearConstInteraction(_)) => true,
            (Coef::Single(_), Coef::Single(_)) => true,
            (Coef::Binary(_), Coef::Binary(_)) => true,
            (Coef::Score2(_), Coef::Score2(_)) => true,
            (Coef::Score3(_), Coef::Score3(_)) => true,
            (Coef::Score4(_), Coef::Score4(_)) => true,
            _ => false,
        }
    }

    // better imple
    // use coef_lr()
    pub fn apply_lr(self, lr: f64) -> Self {
        match self {
            Coef::Linear(x) => Coef::Linear(lr * x),
            Coef::LinearConst((x1, x2)) => Coef::LinearConst((lr * x1, lr * x2)),
            Coef::LinearConstInteraction((x1, x2)) => {
                Coef::LinearConstInteraction((lr * x1, lr * x2))
            }
            Coef::Single(x) => Coef::Single(lr * x),
            Coef::Binary((x1, x2)) => Coef::Binary((lr * x1, lr * x2)),
            Coef::Score2((x1, x2)) => Coef::Score2((lr * x1, lr * x2)),
            Coef::Score3((x1, x2, x3)) => Coef::Score3((lr * x1, lr * x2, lr * x3)),
            Coef::Score4((x1, x2, x3, x4)) => Coef::Score4((lr * x1, lr * x2, lr * x3, lr * x4)),
            Coef::NaN => Coef::NaN,
        }
    }
}

/// For model in wgt file
//#[derive(Debug, Clone, Copy)]
//pub enum ModelType {
//    // FIXME: make modeltype corresponding to all Coef
//    //
//    // so that able to distinguish Coef from .wgt
//    // Coef::Linear, LinearConst
//    Linear,
//    Single,
//    // Coef::Binary, Score2
//    Binary,
//    // Coef::Score3, Score4
//    Free,
//    NaN,
//}

/*
// Which sep should be used, '\t', ' ' ?
impl std::fmt::Display for Coef {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.sida())
    }
}
*/

pub fn coef_format_string(v: f64) -> String {
    format!("{:.7}", v)
}

impl Coef {
    //pub fn model_type(self) -> ModelType {
    //    match self {
    //        Coef::Linear(_) | Coef::LinearConst(_) => ModelType::Linear,
    //        Coef::Single(_) => ModelType::Single,
    //        Coef::Binary(_) | Coef::Score2(_) => ModelType::Binary,
    //        Coef::Score3(_) | Coef::Score4(_) => ModelType::Free,
    //        Coef::NaN => ModelType::NaN,
    //    }
    //}
    //pub fn is_linear(self) -> bool {
    //    match self.model_type() {
    //        ModelType::Linear => true,
    //        //Coef::Linear(_) | Coef::LinearConst(_) => true,
    //        _ => false,
    //    }
    //}
    //pub fn is_binary(self) -> bool {
    //    match self.model_type() {
    //        ModelType::Binary => true,
    //        //Coef::Binary(_) | Coef::Score2(_) => true,
    //        _ => false,
    //    }
    //}
    //pub fn is_free(self) -> bool {
    //    match self.model_type() {
    //        ModelType::Free => true,
    //        _ => false,
    //    }
    //}

    pub fn model_string_for_wgt(self) -> String {
        // v2: safe ver.
        // hopefully <8 chars for tsv
        // should be unique
        match self {
            Coef::Linear(_) => "LINEAR".to_owned(),
            Coef::LinearConst(_) => "LINEARC".to_owned(), // or "LINEARI"="LINEARINTERCEPT"
            Coef::LinearConstInteraction(_) => "LINEARCI".to_owned(), //
            Coef::Single(_) => "SINGLE".to_owned(),
            Coef::Binary(_) => "BINARY".to_owned(),
            Coef::Score2(_) => "FREE2".to_owned(),
            Coef::Score3(_) => "FREE3".to_owned(),
            Coef::Score4(_) => "FREE4".to_owned(),
            Coef::NaN => panic!("Coef is NaN"),
        }

        // ver1, unsafe: model string is not unique and might be confusing
        //match self {
        //    Coef::Linear(_) | Coef::LinearConst(_)|Coef::LinearConstInteraction(_) => "LINEAR".to_owned(),
        //    Coef::Single(_) => "SINGLE".to_owned(),
        //    Coef::Binary(_) | Coef::Score2(_) => "BINARY".to_owned(),
        //    Coef::Score3(_) | Coef::Score4(_) => "FREE".to_owned(),
        //    Coef::NaN => panic!("Coef is NaN"),
        //}
    }

    // check order is not important in ver2
    pub fn new_from_model_name(model_name: &str) -> Option<Self> {
        match model_name {
            "LINEAR" => Some(Coef::Linear(f64::NAN)),
            "LINEARC" => Some(Coef::LinearConst((f64::NAN, f64::NAN))),
            "LINEARCI" => Some(Coef::LinearConstInteraction((f64::NAN, f64::NAN))),
            "SINGLE" => Some(Coef::Single(f64::NAN)),
            "BINARY" => Some(Coef::Binary((f64::NAN, f64::NAN))),
            "FREE2" => Some(Coef::Score2((f64::NAN, f64::NAN))),
            "FREE3" => Some(Coef::Score3((f64::NAN, f64::NAN, f64::NAN))),
            "FREE4" => Some(Coef::Score4((f64::NAN, f64::NAN, f64::NAN, f64::NAN))),
            _ => None,
        }

        //match model_name {
        //    "LINEAR" => Coef::Linear(f64::NAN),
        //    "LINEARC" => Coef::LinearConst((f64::NAN, f64::NAN)),
        //    "LINEARCI" => Coef::LinearConstInteraction((f64::NAN, f64::NAN)),
        //    "SINGLE" => Coef::Single(f64::NAN),
        //    "BINARY" => Coef::Binary((f64::NAN, f64::NAN)),
        //    "FREE2" => Coef::Score2((f64::NAN, f64::NAN)),
        //    "FREE3" => Coef::Score3((f64::NAN, f64::NAN, f64::NAN)),
        //    "FREE4" => Coef::Score4((f64::NAN, f64::NAN, f64::NAN, f64::NAN)),
        //    _ => None,
        //}
    }

    // for safe ver. of model_string()
    //pub fn model_from_string(model_string: &str) -> Self {
    //    match model_string {
    //        "LINEAR" => Coef::Linear(f64::NAN),
    //        "LINEARI" => Coef::LinearConst((f64::NAN, f64::NAN)),
    //        "SINGLE" => Coef::Single(f64::NAN),
    //        "BINARY" => Coef::Binary((f64::NAN, f64::NAN)),
    //        "FREE2" => Coef::Score2((f64::NAN, f64::NAN)),
    //        "FREE3" => Coef::Score3((f64::NAN, f64::NAN, f64::NAN)),
    //        "FREE4" => Coef::Score4((f64::NAN, f64::NAN, f64::NAN, f64::NAN)),
    //        _ => panic!("unknown model: {}", model_string),
    //    }
    //}

    pub fn to_string_hash(self) -> HashMap<String, String> {
        let mut hash = HashMap::new();
        match self {
            Coef::Single(alpha_t) => {
                hash.insert("alpha".to_owned(), coef_format_string(alpha_t));
            }
            Coef::Linear(alpha_t) => {
                hash.insert("alpha".to_owned(), coef_format_string(alpha_t));
            }
            Coef::Binary((const_t, alpha_t)) => {
                hash.insert("const".to_owned(), coef_format_string(const_t));
                hash.insert("alpha".to_owned(), coef_format_string(alpha_t));
            }
            Coef::Score3((s0, s1, s2)) => {
                hash.insert("score0".to_owned(), coef_format_string(s0));
                hash.insert("score1".to_owned(), coef_format_string(s1));
                hash.insert("score2".to_owned(), coef_format_string(s2));
                //hash.insert("scorem".to_owned(), coef_format_string(m));
            }
            Coef::Score4((s0, s1, s2, m)) => {
                hash.insert("score0".to_owned(), coef_format_string(s0));
                hash.insert("score1".to_owned(), coef_format_string(s1));
                hash.insert("score2".to_owned(), coef_format_string(s2));
                hash.insert("scorem".to_owned(), coef_format_string(m));
            }
            Coef::LinearConst((const_t, alpha_t)) => {
                hash.insert("const".to_owned(), coef_format_string(const_t));
                hash.insert("alpha".to_owned(), coef_format_string(alpha_t));
            }
            Coef::LinearConstInteraction((const_t, alpha_t)) => {
                hash.insert("const".to_owned(), coef_format_string(const_t));
                hash.insert("alpha".to_owned(), coef_format_string(alpha_t));
            }
            _ => unimplemented!(),
        }
        hash
    }
}
