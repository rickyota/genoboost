use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Coef {
    // alpha
    // should use Linear(f64)?
    Linear(f64),
    //Linear((f64,)),
    //  (const,alpha)
    LinearConst((f64, f64)),
    // alpha
    // for AdaBoost
    // use Binary((alpha, 0.0))
    Single(f64),
    // FIXME: DO NOT USE Bianry() now; use Score2 instead
    // since cannot distinguish from LinearConst from .wgt file
    // -> if Model column exist, ok to use
    //
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

    // better imple
    // use coef_lr()
    pub fn apply_lr(self, lr: f64) -> Self {
        match self {
            Coef::Linear(x) => Coef::Linear(lr * x),
            Coef::LinearConst((x1, x2)) => Coef::LinearConst((lr * x1, lr * x2)),
            Coef::Single(x) => Coef::Single(lr * x),
            Coef::Binary((x1, x2)) => Coef::Binary((lr * x1, lr * x2)),
            Coef::Score2((x1, x2)) => Coef::Score2((lr * x1, lr * x2)),
            Coef::Score3((x1, x2, x3)) => Coef::Score3((lr * x1, lr * x2, lr * x3)),
            Coef::Score4((x1, x2, x3, x4)) => Coef::Score4((lr * x1, lr * x2, lr * x3, lr * x4)),
            Coef::NaN => Coef::NaN,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum ModelType {
    // FIXME: make modeltype corresponding to all Coef
    //
    // so that able to distinguish Coef from .wgt
    // Coef::Linear, LinearConst
    Linear,
    Single,
    // Coef::Binary, Score2
    Binary,
    // Coef::Score3, Score4
    Free,
    NaN,
}

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
    pub fn model_type(self) -> ModelType {
        match self {
            Coef::Linear(_) | Coef::LinearConst(_) => ModelType::Linear,
            Coef::Single(_) => ModelType::Single,
            Coef::Binary(_) | Coef::Score2(_) => ModelType::Binary,
            Coef::Score3(_) | Coef::Score4(_) => ModelType::Free,
            Coef::NaN => ModelType::NaN,
        }
    }
    pub fn is_linear(self) -> bool {
        match self.model_type() {
            ModelType::Linear => true,
            //Coef::Linear(_) | Coef::LinearConst(_) => true,
            _ => false,
        }
    }
    pub fn is_binary(self) -> bool {
        match self.model_type() {
            ModelType::Binary => true,
            //Coef::Binary(_) | Coef::Score2(_) => true,
            _ => false,
        }
    }
    pub fn is_free(self) -> bool {
        match self.model_type() {
            ModelType::Free => true,
            _ => false,
        }
    }

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
            _ => unimplemented!(),
        }
        hash
    }

    /*
    pub fn to_string_vec(self) -> Vec<String> {
        match self {
            Coef::Linear(alpha_t) => vec![coef_format_string(alpha_t)],
            Coef::Binary((const_t, alpha_t)) => {
                vec![coef_format_string(const_t), coef_format_string(alpha_t)]
            }
            Coef::Score4((s0, s1, s2, m)) => {
                vec![
                    coef_format_string(s0),
                    coef_format_string(s1),
                    coef_format_string(s2),
                    coef_format_string(m),
                ]
            }
            _ => unimplemented!(),
        }
    }
     */
}
