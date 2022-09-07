use std::str::FromStr;

/*
pub trait BoostingParamTrait: Copy {
    fn boosting_param(&self) -> BoostingParam;
    fn iter(&self) -> usize {
        self.boosting_param().iter()
    }
    fn boost_type(&self) -> BoostType {
        self.boosting_param().boost_type()
    }
    fn loss_func(&self) -> LossFunc {
        self.boosting_param().loss_func()
    }
    fn learning_rate(&self) -> Option<f64> {
        self.boosting_param().learning_rate()
    }
    fn eps(&self) -> Eps {
        self.boosting_param().eps()
    }
}
*/

#[derive(PartialEq, Copy, Clone, Debug)]
pub struct BoostParam {
    iter: usize,
    boost_type: BoostType,
    loss_func: LossFunc,
    learning_rate: Option<f64>,
    eps: Eps,
    sample_weight_clip: Option<SampleWeightClip>,
    is_dom_rec: bool,
}

impl BoostParam {
    fn check(self) {
        if (self.boost_type() == BoostType::FreeModelMissing) & self.is_dom_rec {
            panic!("Cannot assign is_dom_and_rec in freemodelmissing")
        }
    }

    pub fn new_str(
        iter: usize,
        boost_type: &str,
        loss_func: &str,
        learning_rate: Option<f64>,
        eps: &str,
        sample_weight_clip: Option<&str>,
        is_dom_rec: bool,
    ) -> Self {
        let bp = BoostParam {
            iter,
            boost_type: BoostType::from_str(boost_type).unwrap(),
            loss_func: LossFunc::from_str(loss_func).unwrap(),
            learning_rate,
            eps: Eps::from_str(eps).unwrap(),
            sample_weight_clip: sample_weight_clip.map(|x| SampleWeightClip::from_str(x).unwrap()),
            is_dom_rec,
        };
        bp.check();
        bp
    }
    // for test
    pub fn new_type1() -> Self {
        BoostParam {
            iter: 100,
            boost_type: BoostType::ConstAda,
            loss_func: LossFunc::Logistic,
            learning_rate: None,
            eps: Eps::MedLarge2AllCell,
            sample_weight_clip: None,
            is_dom_rec: false,
        }
    }
    // for test
    pub fn new_type2() -> Self {
        BoostParam {
            iter: 100,
            boost_type: BoostType::FreeModelMissing,
            loss_func: LossFunc::Logistic,
            learning_rate: None,
            eps: Eps::MedLarge2AllCell,
            sample_weight_clip: None,
            is_dom_rec: false,
        }
    }
    /*     fn boost_param(&self) -> BoostParam {
        *self
    } */
    pub fn iter(&self) -> usize {
        self.iter
    }
    pub fn boost_type(&self) -> BoostType {
        self.boost_type
    }
    pub fn loss_func(&self) -> LossFunc {
        self.loss_func
    }
    pub fn learning_rate(&self) -> Option<f64> {
        self.learning_rate
    }
    pub fn eps(&self) -> Eps {
        self.eps
    }
    pub fn sample_weight_clip(&self) -> Option<SampleWeightClip> {
        self.sample_weight_clip
    }
    pub fn is_dom_rec(&self) -> bool {
        self.is_dom_rec
    }
}

/*
impl BoostingParamTrait for BoostingParam {
    fn boosting_param(&self) -> BoostingParam {
        self.boosting_param()
    }
}
 */

#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
pub enum BoostType {
    Ada, // original adaboost
    ConstAda,
    //Logit, // indicate LossFunc
    // both Modelfree and ModelfreeMissing use ModelfreeClassifier
    //FreeModel,
    FreeModelMissing,
    // AdaLogit
    //Eta(f64), // arg: eta:f64
    //etalogit(f64)
}

impl BoostType {
    pub fn use_missing(self) -> bool {
        match self {
            BoostType::Ada | BoostType::ConstAda => false,
            BoostType::FreeModelMissing => true,
        }
    }
}

impl FromStr for BoostType {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            "ada" => Ok(BoostType::Ada),
            "constada" => Ok(BoostType::ConstAda),
            "genoboost" | "freemodelmissing" => Ok(BoostType::FreeModelMissing),
            _ => Err(format!("Unknown BoostType: {}", str)),
        }
    }
}

#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
pub enum LossFunc {
    Exp,
    Logistic,
}

impl FromStr for LossFunc {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            "exp" => Ok(LossFunc::Exp),
            "logistic" => Ok(LossFunc::Logistic),
            _ => Err(format!("Unknown LossFunc: {}", str)),
        }
    }
}

// add MedLarge2AllCellM : for missing
// add MedLarge2AllCellAlways : for all snvs
/// In Fromstr, assume #case < #cont, but the results are the same even if label is flipped
#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
pub enum Eps {
    MedLarge2AllCell, // ~ median of case; choose larger median
    MedLarge2,        // ~ median of case; choose larger median
    MedSmall,         // ~ meidan of cont; choose smaller median
    Med,
}

impl Eps {
    pub fn allcell(self) -> bool {
        match self {
            Eps::MedLarge2AllCell => true,
            Eps::MedLarge2 | Eps::MedSmall | Eps::Med => false,
        }
    }
}

impl FromStr for Eps {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            "medcase2allcell" | "medlarge2allcell" => Ok(Eps::MedLarge2AllCell),
            "medcont" | "medsmall" => Ok(Eps::MedSmall),
            _ => Err(format!("Unknown Eps: {}", str)),
        }
    }
}

// do not create None, as a member since Option<SampleWeigCLip> is the best
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum SampleWeightClip {
    // clip top x proportion
    Top(f64),
    Both(f64),
}

impl FromStr for SampleWeightClip {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        if str.starts_with("top") {
            let x = str[3..].parse::<f64>().unwrap();
            // make fn?
            assert!(0.0 < x && x < 1.0);
            return Ok(SampleWeightClip::Top(x));
        } else if str.starts_with("both") {
            let x = str[4..].parse::<f64>().unwrap();
            assert!(0.0 < x && x < 1.0);
            return Ok(SampleWeightClip::Both(x));
        }
        Err(format!("Unknown SampleWeightClip: {}", str))
    }
}

#[derive(Debug, Copy, Clone)]
pub enum BoostMethod {
    Classic,      // usual boosting
    Pruning(f64), // arg: clump_r2
    Ss(f64),
}

#[cfg(test)]
mod tests {
    //use super::*;
}
