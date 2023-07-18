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

// for boost_types and learning rates
#[derive(PartialEq, Clone, Debug, Default)]
pub struct BoostParamsTypes {
    // TODO: None for integrate
    //iteration: Option<IterationNumber>,
    iteration: IterationNumber,
    boost_types: Vec<BoostType>,
    //boost_type: BoostType,
    loss_func: LossFunc,
    // TODO: make learning_rate out of BoostParam
    // since this value varies among parameters
    learning_rates: Vec<f64>,
    //learning_rate: Option<f64>,
    eps: Option<Eps>,
    sample_weight_clip: Option<SampleWeightClip>,
    sample_weight_wls_clip: Option<SampleWeightWlsClip>,
    is_dom_rec: bool,
    cov_way: Option<CovWay>,
    batch_way: Option<BatchWay>,
    eff_eps: Option<EffEps>,
}

impl BoostParamsTypes {
    // or itegrate with set_iteration_snv?
    pub fn set_iteration(self, iteration: usize) -> Self {
        Self {
            iteration: IterationNumber::Iteration(iteration),
            ..self
        }
    }

    pub fn set_iteration_snv(self, iteration: usize) -> Self {
        Self {
            iteration: IterationNumber::Snv(iteration),
            ..self
        }
    }

    pub fn set_loss_func(self, loss_func: &str) -> Self {
        Self {
            loss_func: LossFunc::from_str(loss_func).unwrap(),
            ..self
        }
    }

    pub fn set_boost_types(self, boost_types: Vec<String>) -> Self {
        let boost_types = boost_types
            .iter()
            .map(|x| BoostType::from_str(x).unwrap())
            .collect::<Vec<BoostType>>();
        Self {
            boost_types: boost_types,
            ..self
        }
    }

    //pub fn set_boost_type(self, boost_type: &str) -> Self {
    //    Self {
    //        boost_types: vec![BoostType::from_str(boost_type).unwrap()],
    //        ..self
    //    }
    //}

    pub fn set_learning_rates(self, learning_rates: Vec<f64>) -> Self {
        Self {
            learning_rates,
            ..self
        }
    }

    pub fn set_eps(self, eps: Option<&str>) -> Self {
        Self {
            eps: eps.map(|x| Eps::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_eff_eps(self, eff_eps: Option<&str>) -> Self {
        Self {
            eff_eps: eff_eps.map(|x| EffEps::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_sample_weight_clip(self, sample_weight_clip: Option<&str>) -> Self {
        Self {
            sample_weight_clip: sample_weight_clip.map(|x| SampleWeightClip::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_sample_weight_wls_clip(self, sample_weight_wls_clip: Option<&str>) -> Self {
        Self {
            sample_weight_wls_clip: sample_weight_wls_clip
                .map(|x| SampleWeightWlsClip::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_is_dom_rec(self, is_dom_rec: bool) -> Self {
        Self { is_dom_rec, ..self }
    }

    pub fn set_cov_way(self, cov_way: Option<&str>) -> Self {
        Self {
            cov_way: cov_way.map(|x| CovWay::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_batch_way(self, batch_way: Option<&str>) -> Self {
        Self {
            batch_way: batch_way.map(|x| BatchWay::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn param_num(&self) -> usize {
        self.boost_types().len()
        //self.learning_rates().len()
    }

    pub fn param_i(&self, i: usize) -> BoostParams {
        if i >= self.param_num() {
            panic!("Exceed the limit of argument {}", i);
        }
        BoostParams {
            iteration: self.iteration,
            boost_type: self.boost_types()[i],
            loss_func: self.loss_func,
            learning_rates: self.learning_rates().to_vec(),
            eps: self.eps,
            sample_weight_clip: self.sample_weight_clip,
            sample_weight_wls_clip: self.sample_weight_wls_clip,
            is_dom_rec: self.is_dom_rec,
            cov_way: self.cov_way,
            batch_way: self.batch_way,
            eff_eps: self.eff_eps,
        }
    }

    pub fn iteration(&self) -> IterationNumber {
        self.iteration
    }

    pub fn boost_types(&self) -> &[BoostType] {
        &self.boost_types
    }
    /*     pub fn boost_type(&self) -> BoostType {
        self.boost_type
    } */
    pub fn loss_func(&self) -> LossFunc {
        self.loss_func
    }
    //pub fn learning_rate(&self) -> Option<f64> {
    pub fn learning_rates(&self) -> &[f64] {
        &self.learning_rates
    }
    pub fn eps(&self) -> Option<Eps> {
        self.eps
    }
    pub fn sample_weight_clip(&self) -> Option<SampleWeightClip> {
        self.sample_weight_clip
    }
    pub fn sample_weight_wls_clip(&self) -> Option<SampleWeightWlsClip> {
        self.sample_weight_wls_clip
    }
    pub fn is_dom_rec(&self) -> bool {
        self.is_dom_rec
    }
    pub fn cov_way(&self) -> Option<CovWay> {
        self.cov_way
    }
    pub fn batch_way(&self) -> Option<BatchWay> {
        self.batch_way
    }
    pub fn eff_eps(&self) -> Option<EffEps> {
        self.eff_eps
    }

    pub fn into_iter(self) -> BoostParamsTypesIntoIter {
        BoostParamsTypesIntoIter::new(self)
    }
}

pub struct BoostParamsTypesIntoIter {
    boost_params: BoostParamsTypes,
    idx: usize,
}

impl BoostParamsTypesIntoIter {
    pub fn new(boost_params: BoostParamsTypes) -> Self {
        Self {
            boost_params,
            idx: 0,
        }
    }
}

impl Iterator for BoostParamsTypesIntoIter {
    type Item = BoostParams;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx < self.boost_params.param_num() {
            let i = self.idx;
            self.idx += 1;
            return Some(self.boost_params.param_i(i));
        }
        None
    }
}

// for learning rates
#[derive(PartialEq, Clone, Debug, Default)]
pub struct BoostParams {
    iteration: IterationNumber,
    //boost_types: Vec<BoostType>,
    boost_type: BoostType,
    loss_func: LossFunc,
    // TODO: make learning_rate out of BoostParam
    // since this value varies among parameters
    learning_rates: Vec<f64>,
    //learning_rate: Option<f64>,
    eps: Option<Eps>,
    sample_weight_clip: Option<SampleWeightClip>,
    sample_weight_wls_clip: Option<SampleWeightWlsClip>,
    is_dom_rec: bool,
    cov_way: Option<CovWay>,
    batch_way: Option<BatchWay>,
    eff_eps: Option<EffEps>,
}

impl BoostParams {
    pub fn check(&self) {
        //if (self.boost_type() == BoostType::FreeModelMissing) & self.is_dom_rec {
        if self.boost_type().is_freemodelmissing() & self.is_dom_rec {
            panic!("Cannot assign is_dom_and_rec in freemodelmissing")
        }

        if self.boost_type().is_type_ada() & self.eps().is_none() {
            panic!("Cannot assign Eps::None in freemodelmissing")
        }

        if self.boost_type().is_logit() & self.sample_weight_clip().is_none() {
            log::info!("WARNING: Cannot assign SampleWeightClip::None in Logit")
            //panic!("Cannot assign SampleWeightClip::None in Logit")
        }

        if self.boost_type().is_type_logit() & !self.cov_way().unwrap().is_first() {
            // since now ps is not renewed for iteration
            panic!("Cannot assign except CovWay::First in Logit")
        }

        //if self.eff_eps().is_some() {
        //    if (self.eff_eps().unwrap().is_lim_s2())
        //        & (self.boost_type() != BoostType::Logit)
        //        & (self.boost_type() != BoostType::LogitNoMissing)
        //    {
        //        panic!("Cannot assign EffEps::LimS2 other than Logit or LogitNoMissing")
        //    }
        //}
    }


    // or itegrate with set_iteration_snv?
    pub fn set_iteration(self, iteration: usize) -> Self {
        Self {
            iteration: IterationNumber::Iteration(iteration),
            ..self
        }
    }

    pub fn set_iteration_snv(self, iteration: usize) -> Self {
        Self {
            iteration: IterationNumber::Snv(iteration),
            ..self
        }
    }

    pub fn set_loss_func(self, loss_func: &str) -> Self {
        Self {
            loss_func: LossFunc::from_str(loss_func).unwrap(),
            ..self
        }
    }

    //pub fn set_boost_types(self, boost_types: Vec<String>) -> Self {
    //    let boost_types = boost_types
    //        .iter()
    //        .map(|x| BoostType::from_str(x).unwrap())
    //        .collect::<Vec<BoostType>>();
    //    Self {
    //        boost_types: boost_types,
    //        ..self
    //    }
    //}

    // TODO: set_boost_types()
    pub fn set_boost_type(self, boost_type: &str) -> Self {
        Self {
            boost_type: BoostType::from_str(boost_type).unwrap(),
            ..self
        }
    }

    pub fn set_learning_rates(self, learning_rates: Vec<f64>) -> Self {
        Self {
            learning_rates,
            ..self
        }
    }

    pub fn set_eps(self, eps: Option<&str>) -> Self {
        Self {
            eps: eps.map(|x| Eps::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_eff_eps(self, eff_eps: Option<&str>) -> Self {
        Self {
            eff_eps: eff_eps.map(|x| EffEps::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_sample_weight_clip(self, sample_weight_clip: Option<&str>) -> Self {
        Self {
            sample_weight_clip: sample_weight_clip.map(|x| SampleWeightClip::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_sample_weight_wls_clip(self, sample_weight_wls_clip: Option<&str>) -> Self {
        Self {
            sample_weight_wls_clip: sample_weight_wls_clip
                .map(|x| SampleWeightWlsClip::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_is_dom_rec(self, is_dom_rec: bool) -> Self {
        Self { is_dom_rec, ..self }
    }

    pub fn set_cov_way(self, cov_way: Option<&str>) -> Self {
        Self {
            cov_way: cov_way.map(|x| CovWay::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn set_batch_way(self, batch_way: Option<&str>) -> Self {
        Self {
            batch_way: batch_way.map(|x| BatchWay::from_str(x).unwrap()),
            ..self
        }
    }

    pub fn param_num(&self) -> usize {
        self.learning_rates().len()
    }

    pub fn param_i(&self, i: usize) -> BoostParam {
        if i >= self.param_num() {
            panic!("Exceed the limit of argument {}", i);
        }
        BoostParam {
            iteration: self.iteration,
            boost_type: self.boost_type,
            loss_func: self.loss_func,
            learning_rate: self.learning_rates()[i],
            eps: self.eps,
            sample_weight_clip: self.sample_weight_clip,
            sample_weight_wls_clip: self.sample_weight_wls_clip,
            is_dom_rec: self.is_dom_rec,
            cov_way: self.cov_way,
            batch_way: self.batch_way,
            eff_eps: self.eff_eps,
        }
    }

    pub fn param_lr_none(&self) -> BoostParam {
        BoostParam {
            iteration: self.iteration,
            boost_type: self.boost_type,
            loss_func: self.loss_func,
            learning_rate: 1.0,
            eps: self.eps,
            sample_weight_clip: self.sample_weight_clip,
            sample_weight_wls_clip: self.sample_weight_wls_clip,
            is_dom_rec: self.is_dom_rec,
            cov_way: self.cov_way,
            batch_way: self.batch_way,
            eff_eps: self.eff_eps,
        }
    }

    pub fn iteration(&self) -> IterationNumber {
        self.iteration
    }

    /*     pub fn boost_types(&self) -> &[BoostType] {
        &self.boost_types
    } */
    pub fn boost_type(&self) -> BoostType {
        self.boost_type
    }
    pub fn loss_func(&self) -> LossFunc {
        self.loss_func
    }
    //pub fn learning_rate(&self) -> Option<f64> {
    pub fn learning_rates(&self) -> &[f64] {
        &self.learning_rates
    }
    pub fn eps(&self) -> Option<Eps> {
        self.eps
    }
    pub fn sample_weight_clip(&self) -> Option<SampleWeightClip> {
        self.sample_weight_clip
    }
    pub fn sample_weight_wls_clip(&self) -> Option<SampleWeightWlsClip> {
        self.sample_weight_wls_clip
    }
    pub fn is_dom_rec(&self) -> bool {
        self.is_dom_rec
    }
    pub fn cov_way(&self) -> Option<CovWay> {
        self.cov_way
    }
    pub fn batch_way(&self) -> Option<BatchWay> {
        self.batch_way
    }
    pub fn eff_eps(&self) -> Option<EffEps> {
        self.eff_eps
    }

    pub fn into_iter(self) -> BoostParamsIntoIter {
        BoostParamsIntoIter::new(self)
    }
}

pub struct BoostParamsIntoIter {
    boost_params: BoostParams,
    idx: usize,
}

impl BoostParamsIntoIter {
    pub fn new(boost_params: BoostParams) -> Self {
        Self {
            boost_params,
            idx: 0,
        }
    }
}

impl Iterator for BoostParamsIntoIter {
    type Item = BoostParam;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx < self.boost_params.param_num() {
            let i = self.idx;
            self.idx += 1;
            return Some(self.boost_params.param_i(i));
        }
        None
    }
}

#[derive(PartialEq, Copy, Clone, Debug, Default)]
pub struct BoostParam {
    iteration: IterationNumber,
    boost_type: BoostType,
    loss_func: LossFunc,
    // TODO: make learning_rate out of BoostParam
    // since this value varies among parameters
    learning_rate: f64,
    //learning_rate: Option<f64>,
    eps: Option<Eps>,
    sample_weight_clip: Option<SampleWeightClip>,
    sample_weight_wls_clip: Option<SampleWeightWlsClip>,
    is_dom_rec: bool,
    cov_way: Option<CovWay>,
    batch_way: Option<BatchWay>,
    eff_eps: Option<EffEps>,
}

impl BoostParam {
    // TODO: use
    pub fn check(&self) {
        //if (self.boost_type() == BoostType::FreeModelMissing) & self.is_dom_rec {
        if self.boost_type().is_freemodelmissing() & self.is_dom_rec {
            panic!("Cannot assign is_dom_and_rec in freemodelmissing")
        }

        if self.boost_type().is_type_ada() & self.eps().is_none() {
            panic!("Cannot assign Eps::None in freemodelmissing")
        }

        if self.boost_type().is_logit() & self.sample_weight_clip().is_none() {
            log::info!("WARNING: Cannot assign SampleWeightClip::None in Logit")
            //panic!("Cannot assign SampleWeightClip::None in Logit")
        }

        if self.boost_type().is_type_logit() & !self.cov_way().unwrap().is_first() {
            // since now ps is not renewed for iteration
            panic!("Cannot assign except CovWay::First in Logit")
        }

        //if self.eff_eps().is_some() {
        //    if (self.eff_eps().unwrap().is_lim_s2())
        //        & (self.boost_type() != BoostType::Logit)
        //        & (self.boost_type() != BoostType::LogitNoMissing)
        //    {
        //        panic!("Cannot assign EffEps::LimS2 other than Logit or LogitNoMissing")
        //    }
        //}
    }

    // for test
    pub fn new_type1() -> Self {
        BoostParam {
            iteration: IterationNumber::Iteration(100),
            boost_type: BoostType::ConstAda,
            loss_func: LossFunc::Logistic,
            learning_rate: 1.0,
            //learning_rate: None,
            //eps: None,
            eps: Some(Eps::MedLarge2AllCell),
            sample_weight_clip: None,
            sample_weight_wls_clip: None,
            is_dom_rec: false,
            cov_way: None,
            batch_way: None,
            eff_eps: None,
        }
    }
    // for test
    pub fn new_type2() -> Self {
        BoostParam {
            iteration: IterationNumber::Iteration(100),
            boost_type: BoostType::FreeModelMissing,
            loss_func: LossFunc::Logistic,
            learning_rate: 1.0,
            //learning_rate: None,
            //eps: None,
            eps: Some(Eps::MedLarge2AllCell),
            sample_weight_clip: None,
            sample_weight_wls_clip: None,
            is_dom_rec: false,
            cov_way: None,
            batch_way: None,
            eff_eps: None,
        }
    }
    /*     fn boost_param(&self) -> BoostParam {
        *self
    } */

    pub fn iteration(&self) -> IterationNumber {
        self.iteration
    }

    //pub fn iter(&self) -> usize {
    //    self.iter
    //}
    pub fn boost_type(&self) -> BoostType {
        self.boost_type
    }
    pub fn loss_func(&self) -> LossFunc {
        self.loss_func
    }
    //pub fn learning_rate(&self) -> Option<f64> {
    pub fn learning_rate(&self) -> f64 {
        self.learning_rate
    }
    pub fn eps(&self) -> Option<Eps> {
        self.eps
    }
    pub fn sample_weight_clip(&self) -> Option<SampleWeightClip> {
        self.sample_weight_clip
    }
    pub fn sample_weight_wls_clip(&self) -> Option<SampleWeightWlsClip> {
        self.sample_weight_wls_clip
    }
    pub fn is_dom_rec(&self) -> bool {
        self.is_dom_rec
    }
    pub fn cov_way(&self) -> Option<CovWay> {
        self.cov_way
    }
    pub fn batch_way(&self) -> Option<BatchWay> {
        self.batch_way
    }
    pub fn eff_eps(&self) -> Option<EffEps> {
        self.eff_eps
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
    Logit, // LogitBoost
    LogitNoMissing,
    LogitAdd, // Additive Logit
    // both Modelfree and ModelfreeMissing use ModelfreeClassifier
    //FreeModel,
    FreeModelMissing,
    // AdaLogit
    //Eta(f64), // arg: eta:f64
    //etalogit(f64)
}

impl BoostType {
    pub fn is_freemodelmissing(self) -> bool {
        match self {
            BoostType::FreeModelMissing => true,
            BoostType::Ada
            | BoostType::ConstAda
            | BoostType::Logit
            | BoostType::LogitNoMissing
            | BoostType::LogitAdd => false,
        }
    }

    pub fn is_logit(self) -> bool {
        match self {
            BoostType::Logit | BoostType::LogitNoMissing => true,
            BoostType::Ada
            | BoostType::ConstAda
            | BoostType::LogitAdd
            | BoostType::FreeModelMissing => false,
        }
    }

    // allow missing value in dataset
    pub fn use_missing(self) -> bool {
        match self {
            BoostType::Ada
            | BoostType::ConstAda
            | BoostType::LogitAdd
            | BoostType::LogitNoMissing => false,
            BoostType::FreeModelMissing | BoostType::Logit => true,
        }
    }

    // TODO: better
    pub fn is_type_ada(self) -> bool {
        match self {
            BoostType::Ada | BoostType::ConstAda | BoostType::FreeModelMissing => true,
            BoostType::Logit | BoostType::LogitNoMissing | BoostType::LogitAdd => false,
        }
    }

    pub fn is_type_logit(self) -> bool {
        match self {
            BoostType::Ada | BoostType::ConstAda | BoostType::FreeModelMissing => false,
            BoostType::Logit | BoostType::LogitNoMissing | BoostType::LogitAdd => true,
        }
    }

    // for directory name of integrate
    pub fn genetic_model(self) -> String {
        match self {
            BoostType::Logit | BoostType::LogitNoMissing => "nonadd".to_string(),
            BoostType::LogitAdd => "add".to_string(),
            _ => unimplemented!(),
        }
    }
}

impl Default for BoostType {
    fn default() -> Self {
        Self::LogitNoMissing
    }
}

impl FromStr for BoostType {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            "ada" => Ok(BoostType::Ada),
            "constada" => Ok(BoostType::ConstAda),
            "freemodelmissing" => Ok(BoostType::FreeModelMissing),
            //"genoboost" | "freemodelmissing" => Ok(BoostType::FreeModelMissing),
            "logit" => Ok(BoostType::Logit),
            "nonadd" | "logitnomissing" => Ok(BoostType::LogitNoMissing),
            "add" | "logitadd" => Ok(BoostType::LogitAdd),
            _ => Err(format!("Unknown BoostType: {}", str)),
        }
    }
}

#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
pub enum LossFunc {
    Exp,
    Logistic,
    //Logit,
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

impl Default for LossFunc {
    fn default() -> Self {
        Self::Logistic
    }
}

// add MedLarge2AllCellM : for missing
// add MedLarge2AllCellAlways : for all snvs
/// In Fromstr, assume #case < #cont, but the results are the same even if label is flipped
#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
pub enum Eps {
    MedLargeAllCellAllSnvs, // ~ median of case; for all snvs even if no cells are 0
    MedLarge2AllCellAllSnvs, // ~ median of case; for all snvs even if no cells are 0
    MedLarge2AllCell,       // ~ median of case; choose larger median
    MedLarge2,              // ~ median of case; choose larger median
    MedSmall,               // ~ meidan of cont; choose smaller median
    Med,
    Dom, // make w2 and w1 together
         //None, // for Logit, TODO: Option<Eps>
}

impl Eps {
    // TODO -> is_none?
    /*     pub fn none(self) -> bool {
        match self {
            Eps::None => true,
            _ => false,
        }
    } */

    pub fn dom(self) -> bool {
        match self {
            Eps::Dom => true,
            _ => false,
        }
    }

    pub fn allcell(self) -> bool {
        match self {
            Eps::MedLarge2AllCell | Eps::MedLargeAllCellAllSnvs | Eps::MedLarge2AllCellAllSnvs => {
                true
            }
            //Eps::MedLarge2 | Eps::MedSmall | Eps::Med => false,
            _ => false,
        }
    }

    pub fn allsnv(self) -> bool {
        match self {
            Eps::MedLarge2AllCellAllSnvs | Eps::MedLargeAllCellAllSnvs => true,
            //Eps::MedLarge2AllCell | Eps::MedLarge2 | Eps::MedSmall | Eps::Med => false,
            _ => false,
        }
    }
}

impl FromStr for Eps {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        match str {
            "medcase2allcellsnvs" | "medlarge2allcellsnvs" => Ok(Eps::MedLarge2AllCellAllSnvs),
            "medcase2allcell" | "medlarge2allcell" => Ok(Eps::MedLarge2AllCell),
            "medcase2" | "medlarge2" => Ok(Eps::MedLarge2),
            "medcont" | "medsmall" => Ok(Eps::MedSmall),
            "medcaseallcellsnvs" | "medlargeallcellsnvs" => Ok(Eps::MedLargeAllCellAllSnvs),
            "dom" => Ok(Eps::Dom),
            //"none" => Ok(Eps::None),
            _ => Err(format!("Unknown Eps: {}", str)),
        }
    }
}

// TODO: rename S2->Score since limit not only s2
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum EffEps {
    LimScore(usize, f64), // lim abs(s2)< LimS2.1 * max(abs(s1), abs(s0)) when d2 or n2 <=LimS2.0
    LimScoreProp(f64, f64),
    LimS2AddProp(f64),
    LimS2GmodelProp(f64, f64, f64), // prop, rec_pos_ratio (rs1 and rs2 are same sign), rec_neg_ratio(rs1 and rs2 are different sign)
    LimS12GmodelProp(f64, f64, f64, f64, f64), // limit s1 too, rec_pos_ratio, rec_neg_ratio, het_pos_ratio, het_neg_ratio,
    LimS12GmodelPropOnUpdate(f64, f64, f64, f64, f64), // limit s1 too, rec_pos_ratio, rec_neg_ratio, het_pos_ratio, het_neg_ratio,
    LimS2GmodelOverProp(f64, f64, f64, f64), //  rec_pos_ratio, rec_neg_ratio, het_pos_ratio
    LimS2GmodelOverPropOnUpdate(f64, f64, f64, f64), // ng
    LimS2GmodelOverKeepSignProp(f64, f64, f64, f64), // diferent het_pos_ratio def from LimS2GmodelOverProp
    LimS2GmodelBorderProp(f64, f64, f64, f64, f64), // prop, thre: rec-add-dom-hetonly-rec // 0.25, 0.75, 1.25, -0.25
}

impl FromStr for EffEps {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        if str.starts_with("lims2_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2");

            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            return Ok(EffEps::LimScore(x, y));
        } else if str.starts_with("lims2prop_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2prop");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            return Ok(EffEps::LimScoreProp(x, y));
        } else if str.starts_with("lims2gmodelprop_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2gmodelprop");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            let z = str_iter.next().unwrap().parse::<f64>().unwrap();

            if y < 0.0 {
                panic!("y should be 0.0 or positive.")
            }
            if z < 0.0 {
                panic!("z should be 0.0 or positive.")
            }

            // so that original genetic model reserves
            let y = 1.001 * y;
            let z = 1.001 * z;

            assert_eq!(str_iter.next(), None);

            return Ok(EffEps::LimS2GmodelProp(x, y, z));
        } else if str.starts_with("lims12gmodelprop_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims12gmodelprop");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            let z = str_iter.next().unwrap().parse::<f64>().unwrap();
            let u = str_iter.next().unwrap().parse::<f64>().unwrap();
            let v = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            if y < 0.0 {
                panic!("y should be 0.0 or positive.")
            }
            if z < 0.0 {
                panic!("z should be 0.0 or positive.")
            }
            if u < 0.0 {
                panic!("u should be 0.0 or positive.")
            }
            if v < 0.0 {
                panic!("v should be 0.0 or positive.")
            }

            // so that original genetic model reserves
            let y = 1.001 * y;
            let z = 1.001 * z;
            let u = 1.001 * u;
            let v = 1.001 * v;

            return Ok(EffEps::LimS12GmodelProp(x, y, z, u, v));
        } else if str.starts_with("lims12gmodelproponupdate_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims12gmodelproponupdate");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            let z = str_iter.next().unwrap().parse::<f64>().unwrap();
            let u = str_iter.next().unwrap().parse::<f64>().unwrap();
            let v = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            if y < 0.0 {
                panic!("y should be 0.0 or positive.")
            }
            if z < 0.0 {
                panic!("z should be 0.0 or positive.")
            }
            if u < 0.0 {
                panic!("u should be 0.0 or positive.")
            }
            if v < 0.0 {
                panic!("v should be 0.0 or positive.")
            }

            // so that original genetic model reserves
            let y = 1.001 * y;
            let z = 1.001 * z;
            let u = 1.001 * u;
            let v = 1.001 * v;

            return Ok(EffEps::LimS12GmodelPropOnUpdate(x, y, z, u, v));
        } else if str.starts_with("lims2gmodeloverprop_") {
            // similar to LimS12GmodelProp but adjust s2 only for rec and overdom
            // ex. (4,4,1.25)
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2gmodeloverprop");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            let z = str_iter.next().unwrap().parse::<f64>().unwrap();
            let u = str_iter.next().unwrap().parse::<f64>().unwrap();
            //let v = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            //if y < 0.0 | y<1.0 {
            if y < 1.0 {
                panic!("y should be 0.0 or positive.")
            }
            //if z < 0.0 | z<1.0{
            if z < 1.0 {
                panic!("z should be 0.0 or positive.")
            }
            //if u < 0.0 | u<1.0{
            if u < 1.0 {
                panic!("u should be 0.0 or positive.")
            }
            //if v < 0.0 {
            //    panic!("v should be 0.0 or positive.")
            //}

            // so that original genetic model reserves
            let y = 1.001 * y;
            let z = 1.001 * z;
            let u = 1.001 * u;
            //let v = 1.001 * v;

            return Ok(EffEps::LimS2GmodelOverProp(x, y, z, u));
        } else if str.starts_with("lims2gmodeloverproponupdate_") {
            // only lim2 on updating score not calculating loss
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2gmodeloverproponupdate");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            let z = str_iter.next().unwrap().parse::<f64>().unwrap();
            let u = str_iter.next().unwrap().parse::<f64>().unwrap();
            //let v = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            if y < 0.0 {
                panic!("y should be 0.0 or positive.")
            }
            if z < 0.0 {
                panic!("z should be 0.0 or positive.")
            }
            if u < 0.0 {
                panic!("u should be 0.0 or positive.")
            }
            //if v < 0.0 {
            //    panic!("v should be 0.0 or positive.")
            //}

            // so that original genetic model reserves
            let y = 1.001 * y;
            let z = 1.001 * z;
            let u = 1.001 * u;
            //let v = 1.001 * v;

            return Ok(EffEps::LimS2GmodelOverPropOnUpdate(x, y, z, u));
        } else if str.starts_with("lims2gmodeloverkeepsignprop_") {
            // similar to LimS12GmodelProp but adjust s2 only for rec and overdom
            // CAUTION: different def of u from LimS2GmodelOverProp (=1/u)
            // ex. (4,4,0.8)
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2gmodeloverkeepsignprop");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            let z = str_iter.next().unwrap().parse::<f64>().unwrap();
            let u = str_iter.next().unwrap().parse::<f64>().unwrap();
            //let v = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            if y < 1.0 {
                panic!("y should be 0.0 or positive.")
            }
            if z < 1.0 {
                panic!("z should be 0.0 or positive.")
            }
            if (u < 0.0) | (u > 1.0) {
                panic!("u should be 0.0 or positive.")
            }
            //if v < 0.0 {
            //    panic!("v should be 0.0 or positive.")
            //}

            // so that original genetic model reserves
            let y = 1.001 * y;
            let z = 1.001 * z;
            let u = 0.999 * u;
            //let u = 1.001 * u;

            return Ok(EffEps::LimS2GmodelOverKeepSignProp(x, y, z, u));
        } else if str.starts_with("lims2gmodeloverkeepsignprop-") {
            // similar to LimS12GmodelProp but adjust s2 only for rec and overdom
            // CAUTION: different def of u from LimS2GmodelOverProp (=1/u)
            // ex. (4,4,0.8)
            let mut str_iter = str.split("-");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2gmodeloverkeepsignprop");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            let z = str_iter.next().unwrap().parse::<f64>().unwrap();
            let u = str_iter.next().unwrap().parse::<f64>().unwrap();
            //let v = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            if y < 1.0 {
                panic!("y should be 0.0 or positive.")
            }
            if z < 1.0 {
                panic!("z should be 0.0 or positive.")
            }
            if (u < 0.0) | (u > 1.0) {
                panic!("u should be 0.0 or positive.")
            }
            //if v < 0.0 {
            //    panic!("v should be 0.0 or positive.")
            //}

            // so that original genetic model reserves
            let y = 1.001 * y;
            let z = 1.001 * z;
            let u = 0.999 * u;
            //let u = 1.001 * u;

            return Ok(EffEps::LimS2GmodelOverKeepSignProp(x, y, z, u));
        } else if str.starts_with("lims2gmodelborderprop_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2gmodelborderprop");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            let y = str_iter.next().unwrap().parse::<f64>().unwrap();
            let z = str_iter.next().unwrap().parse::<f64>().unwrap();
            let u = str_iter.next().unwrap().parse::<f64>().unwrap();
            let w = str_iter.next().unwrap().parse::<f64>().unwrap();
            assert_eq!(str_iter.next(), None);

            return Ok(EffEps::LimS2GmodelBorderProp(x, y, z, u, w));
        } else if str.starts_with("lims2onlyaddprop_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "lims2onlyaddprop");

            let x = str_iter.next().unwrap().parse::<f64>().unwrap();
            // % -> prop
            let x = x / 100.0;
            assert_eq!(str_iter.next(), None);

            return Ok(EffEps::LimS2AddProp(x));
        } else {
            Err(format!("Unknown EffEps: {}", str))
        }
    }
}

impl EffEps {
    pub fn is_on_update(self) -> bool {
        match self {
            Self::LimS12GmodelPropOnUpdate(..) | Self::LimS2GmodelOverPropOnUpdate(..) => true,
            Self::LimScore(..)
            | Self::LimScoreProp(..)
            | Self::LimS2GmodelProp(..)
            | Self::LimS12GmodelProp(..)
            | Self::LimS2GmodelOverProp(..)
            | Self::LimS2GmodelOverKeepSignProp(..)
            | Self::LimS2GmodelBorderProp(..)
            | Self::LimS2AddProp(..) => false,
        }
    }
}

// do not create None, as a member since Option<SampleWeightClip> is the best
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
            if x <= 0.0 || x > 1.0 {
                panic!("Sample weight clip should be 0.0<x<=1.0.")
            }
            return Ok(SampleWeightClip::Top(x));
        } else if str.starts_with("both") {
            let x = str[4..].parse::<f64>().unwrap();
            if x <= 0.0 || x > 1.0 {
                panic!("Sample weight clip should be 0.0<x<=1.0.")
            }
            return Ok(SampleWeightClip::Both(x));
        }
        Err(format!("Unknown SampleWeightClip: {}", str))
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum SampleWeightWlsClip {
    // clip top x proportion
    Top(f64),
    Both(f64),
}

impl FromStr for SampleWeightWlsClip {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        if str.starts_with("top") {
            let x = str[3..].parse::<f64>().unwrap();
            if x <= 0.0 || x > 1.0 {
                panic!("Sample weight clip should be 0.0<x<=1.0.")
            }
            return Ok(SampleWeightWlsClip::Top(x));
        } else if str.starts_with("both") {
            let x = str[4..].parse::<f64>().unwrap();
            if x <= 0.0 || x > 1.0 {
                panic!("Sample weight clip should be 0.0<x<=1.0.")
            }
            return Ok(SampleWeightWlsClip::Both(x));
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

// TODO: mv in boostparam?
#[derive(PartialEq, Copy, Clone, Debug)]
pub enum IterationNumber {
    Iteration(usize),
    Snv(usize),
}

impl Default for IterationNumber {
    fn default() -> Self {
        IterationNumber::Snv(100)
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum CovWay {
    First,
    Every(usize),
}

impl CovWay {
    pub fn is_first(self) -> bool {
        match self {
            CovWay::First => true,
            CovWay::Every(_) => false,
        }
    }
}

impl FromStr for CovWay {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        if str == "first" {
            return Ok(CovWay::First);
        }

        // FIXME: make str="every10"
        let x = str.parse::<usize>().unwrap();
        if x == 0 {
            panic!("cov should be >0.")
        }
        return Ok(CovWay::Every(x));
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum BatchWay {
    Every(usize),
    EveryMaxIter(usize, usize, usize), // (batch_size, batch_max_iter_lr1, batch_max_lr_others)
    Fix(usize, usize, usize),          // (batch_size, batch_iter_lr1, batch_lr_others)
}

impl BatchWay {
    pub fn batch_size(self) -> usize {
        match self {
            Self::Every(x) => x,
            Self::EveryMaxIter(x, ..) | Self::Fix(x, ..) => x,
            //Self::EveryMaxIter(x, ..) => x,
            //Self::Fix(x,..) => x,
        }
    }
    pub fn use_comp_loss(self) -> bool {
        match self {
            Self::Every(..) | Self::EveryMaxIter(..) => true,
            Self::Fix(..) => false,
        }
    }
    pub fn batch_max_iter(self, lr: f64) -> usize {
        match self {
            Self::Every(x) => x,
            Self::EveryMaxIter(_, x1, x2) | Self::Fix(_, x1, x2) => {
                if lr > 0.6f64 {
                    x1
                } else {
                    x2
                }
            }
        }
    }
}

impl FromStr for BatchWay {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        // TODO: deprecate '_'
        if str.starts_with("every_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "every");

            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            return Ok(BatchWay::Every(x));
        } else if str.starts_with("everymaxiter_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "everymaxiter");

            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<usize>().unwrap();
            let z = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            return Ok(BatchWay::EveryMaxIter(x, y, z));
        } else if str.starts_with("fix_") {
            let mut str_iter = str.split("_");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "fix");

            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<usize>().unwrap();
            let z = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            return Ok(BatchWay::Fix(x, y, z));
        } else if str.starts_with("fix-") {
            let mut str_iter = str.split("-");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "fix");

            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<usize>().unwrap();
            let z = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            return Ok(BatchWay::Fix(x, y, z));
        } else {
            // TODO: deprecate
            // assume usize
            let x = str.parse::<usize>().unwrap();
            if x == 0 {
                panic!("Batchway should be >0.")
            }
            return Ok(BatchWay::Every(x));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boosting_params() {
        let lrs = vec![1.0, 0.1];
        let boost_params = BoostParams::default().set_learning_rates(lrs);

        let mut boost_iter = boost_params.into_iter();

        assert_eq!(boost_iter.next().unwrap().learning_rate(), 1.0);
        assert_eq!(boost_iter.next().unwrap().learning_rate(), 0.1);
    }
}
