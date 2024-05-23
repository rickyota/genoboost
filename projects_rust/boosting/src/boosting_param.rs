//  When adding one parameter, modify the following 4 places
// 1. BoostParamCommon.set_batch_way()
// 2. BatchWay::from_str()
// 3. BoostParamCommonTrait.batch_way()
// 4. BoostParamCommon.batch_way()

use core::panic;
use std::{collections::HashSet, str::FromStr};

use genetics::{LdCriteria, SnvId};

// for boost_types and learning rates
#[derive(PartialEq, Clone, Debug, Default)]
pub struct BoostParamsTypes {
    boost_types: Vec<BoostType>,
    learning_rates: Vec<f64>,
    boost_param_common: BoostParamCommon,
}

impl BoostParamCommonTrait for BoostParamsTypes {
    fn boost_param_common(&self) -> &BoostParamCommon {
        &self.boost_param_common
    }
}

impl BoostParamsTypes {
    pub fn set_boost_types(self, boost_types: Vec<String>) -> Self {
        if HashSet::<&String>::from_iter(boost_types.iter()).len() != boost_types.len() {
            panic!("BoostTypes must be unique");
        }

        let boost_types = boost_types
            .iter()
            .map(|x| BoostType::from_str(x).unwrap())
            .collect::<Vec<BoostType>>();

        Self {
            boost_types,
            ..self
        }
    }
    pub fn set_learning_rates(self, learning_rates: Vec<f64>) -> Self {
        Self {
            learning_rates,
            ..self
        }
    }
    pub fn set_boost_param_common(self, boost_param_common: BoostParamCommon) -> Self {
        Self {
            boost_param_common,
            ..self
        }
    }

    pub fn param_num(&self) -> usize {
        self.boost_types().len()
    }

    pub fn param_i(&self, i: usize) -> BoostParamLrs {
        if i >= self.param_num() {
            panic!("Exceed the limit of argument {}", i);
        }
        BoostParamLrs {
            boost_type: self.boost_types()[i],
            learning_rates: self.learning_rates().to_vec(),
            boost_param_common: self.boost_param_common.clone(),
        }
    }

    pub fn boost_types(&self) -> &[BoostType] {
        &self.boost_types
    }
    pub fn learning_rates(&self) -> &[f64] {
        &self.learning_rates
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
    type Item = BoostParamLrs;

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
pub struct BoostParamLrs {
    boost_type: BoostType,
    learning_rates: Vec<f64>,
    boost_param_common: BoostParamCommon,
}

impl BoostParamLrs {
    pub fn check(&self) {
        //if (self.boost_type() == BoostType::FreeModelMissing) & self.is_dom_rec {
        if self.boost_type().is_freemodelmissing() && self.is_dom_rec() {
            panic!("Cannot assign is_dom_and_rec in freemodelmissing")
        }

        // use effeps instead
        //if self.boost_type().is_type_ada() & self.eps().is_none() {
        //    panic!("Cannot assign Eps::None in freemodelmissing")
        //}

        //if self.boost_type().is_logit() & self.sample_weight_clip().is_none() {
        //    log::info!("WARNING: Cannot assign SampleWeightClip::None in Logit")
        //    //panic!("Cannot assign SampleWeightClip::None in Logit")
        //}

        if self.boost_type().is_type_logit() && !self.cov_way().unwrap().is_first() {
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

        if self.boost_type() == BoostType::LogitCommon
            && self.maf_threshold_logit_common().is_none()
        {
            panic!("maf_threshold_logit_common must be assigned in LogitCommon")
        }
    }

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
    pub fn set_boost_param_common(self, boost_param_common: BoostParamCommon) -> Self {
        Self {
            boost_param_common,
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
            boost_type: self.boost_type,
            learning_rate: self.learning_rates()[i],
            boost_param_common: self.boost_param_common.clone(),
        }
    }

    pub fn param_lr_none(&self) -> BoostParam {
        BoostParam {
            boost_type: self.boost_type,
            learning_rate: 1.0,
            boost_param_common: self.boost_param_common.clone(),
        }
    }

    pub fn boost_type(&self) -> BoostType {
        self.boost_type
    }
    pub fn learning_rates(&self) -> &[f64] {
        &self.learning_rates
    }

    pub fn into_iter(self) -> BoostParamsIntoIter {
        BoostParamsIntoIter::new(self)
    }
}

impl BoostParamCommonTrait for BoostParamLrs {
    fn boost_param_common(&self) -> &BoostParamCommon {
        &self.boost_param_common
    }
}

pub struct BoostParamsIntoIter {
    boost_params: BoostParamLrs,
    idx: usize,
}

impl BoostParamsIntoIter {
    pub fn new(boost_params: BoostParamLrs) -> Self {
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

#[derive(PartialEq, Clone, Debug, Default)]
pub struct BoostParam {
    boost_type: BoostType,
    learning_rate: f64,
    boost_param_common: BoostParamCommon,
}

impl BoostParamCommonTrait for BoostParam {
    fn boost_param_common(&self) -> &BoostParamCommon {
        &self.boost_param_common
    }
}

impl BoostParam {
    // TODO: use
    pub fn check(&self) {
        if self.boost_type().is_freemodelmissing() && self.is_dom_rec() {
            panic!("Cannot assign is_dom_and_rec in freemodelmissing")
        }

        if self.boost_type().is_type_ada() && self.eps().is_none() {
            panic!("Cannot assign Eps::None in freemodelmissing")
        }

        //if self.boost_type().is_logit() & self.sample_weight_clip().is_none() {
        //    log::info!("WARNING: Cannot assign SampleWeightClip::None in Logit")
        //    //panic!("Cannot assign SampleWeightClip::None in Logit")
        //}

        if self.boost_type().is_type_logit() && !self.cov_way().unwrap().is_first() {
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
    //pub fn new_type1() -> Self {
    //    BoostParam {
    //        boost_type: BoostType::ConstAda,
    //        learning_rate: 1.0,
    //        boost_param_common: BoostParamCommon {
    //            iteration: Some(IterationNumber::Iteration(100)),
    //            loss_func: LossFunc::Logistic,
    //            eps: Some(Eps::MedLarge2AllCell),
    //            sample_weight_clip: None,
    //            sample_weight_wls_clip: None,
    //            is_dom_rec: false,
    //            cov_way: None,
    //            batch_way: None,
    //            batch_interaction_way: None,
    //            eff_eps: None,
    //            is_monitor: false,
    //            monitor_nsnvs: None,
    //            mhc_region: None,
    //            maf_threshold_logit_common: None,
    //            interaction_way: None,
    //            ld_criteria: None,
    //        },
    //    }
    //}
    //// for test
    //pub fn new_type2() -> Self {
    //    BoostParam {
    //        boost_type: BoostType::FreeModelMissing,
    //        learning_rate: 1.0,
    //        boost_param_common: BoostParamCommon {
    //            iteration: Some(IterationNumber::Iteration(100)),
    //            loss_func: LossFunc::Logistic,
    //            eps: Some(Eps::MedLarge2AllCell),
    //            sample_weight_clip: None,
    //            sample_weight_wls_clip: None,
    //            is_dom_rec: false,
    //            cov_way: None,
    //            batch_way: None,
    //            batch_interaction_way: None,
    //            eff_eps: None,
    //            is_monitor: false,
    //            monitor_nsnvs: None,
    //            mhc_region: None,
    //            maf_threshold_logit_common: None,
    //            interaction_way: None,
    //            ld_criteria: None,
    //        },
    //    }
    //}

    // for test
    pub fn new_type3() -> Self {
        let boost_param_common = BoostParamCommon::default().set_iteration(Some(100), None);
        BoostParam {
            boost_type: BoostType::FreeModelMissing,
            learning_rate: 1.0,
            boost_param_common,
        }
    }

    pub fn boost_type(&self) -> BoostType {
        self.boost_type
    }
    pub fn learning_rate(&self) -> f64 {
        self.learning_rate
    }
}

#[derive(PartialEq, Clone, Debug, Default)]
pub struct BoostParamCommon {
    iteration: Option<IterationNumber>,
    loss_func: LossFunc,
    eps: Option<Eps>,
    sample_weight_clip: Option<SampleWeightClip>,
    sample_weight_wls_clip: Option<SampleWeightWlsClip>,
    is_dom_rec: bool,
    cov_way: Option<CovWay>,
    batch_way: Option<BatchWay>,
    batch_interaction_way: Option<BatchInteractionWay>,
    eff_eps: Option<EffEps>,
    is_monitor: bool,
    monitor_nsnvs: Option<Vec<usize>>,
    mhc_region: Option<(SnvId, SnvId)>,
    maf_threshold_logit_common: Option<f64>, // for BoostType::LogitCommon
    interaction_way: Option<InteractionWay>,
    ld_criteria: Option<LdCriteria>,
    acc_metic: Option<AccMetric>,
}

pub trait BoostParamCommonTrait {
    fn boost_param_common(&self) -> &BoostParamCommon;

    //fn boost_param_common_mut(&mut self) -> &mut BoostParamCommon;

    fn iteration(&self) -> Option<IterationNumber> {
        self.boost_param_common().iteration()
    }

    fn loss_func(&self) -> LossFunc {
        self.boost_param_common().loss_func()
    }
    fn eps(&self) -> Option<Eps> {
        self.boost_param_common().eps()
    }
    fn sample_weight_clip(&self) -> Option<SampleWeightClip> {
        self.boost_param_common().sample_weight_clip()
    }
    fn sample_weight_wls_clip(&self) -> Option<SampleWeightWlsClip> {
        self.boost_param_common().sample_weight_wls_clip()
    }
    fn is_dom_rec(&self) -> bool {
        self.boost_param_common().is_dom_rec()
    }
    fn cov_way(&self) -> Option<CovWay> {
        self.boost_param_common().cov_way()
    }
    fn batch_way(&self) -> Option<BatchWay> {
        self.boost_param_common().batch_way()
    }
    fn batch_interaction_way(&self) -> Option<BatchInteractionWay> {
        self.boost_param_common().batch_interaction_way()
    }
    fn eff_eps(&self) -> Option<EffEps> {
        self.boost_param_common().eff_eps()
    }
    fn is_monitor(&self) -> bool {
        self.boost_param_common().is_monitor()
    }
    fn monitor_nsnvs(&self) -> Option<&[usize]> {
        self.boost_param_common().monitor_nsnvs()
    }
    fn mhc_region(&self) -> Option<&(SnvId, SnvId)> {
        self.boost_param_common().mhc_region()
    }

    fn maf_threshold_logit_common(&self) -> Option<f64> {
        self.boost_param_common().maf_threshold_logit_common()
    }

    fn interaction_way(&self) -> Option<InteractionWay> {
        self.boost_param_common().interaction_way()
    }

    fn ld_criteria(&self) -> Option<LdCriteria> {
        self.boost_param_common().ld_criteria()
    }

    fn acc_metric(&self) -> Option<AccMetric> {
        self.boost_param_common().acc_metric()
    }

    // ng: it is hard to use set_() in trait; use in BoostParamCommon
    //fn set_eps(self, eps: Option<&str>) -> BoostParamCommon {
    //    BoostParamCommon {
    //        eps: eps.map(|x| Eps::from_str(x).unwrap()),
    //        ..self
    //    }
    //}
}

impl BoostParamCommonTrait for BoostParamCommon {
    fn boost_param_common(&self) -> &BoostParamCommon {
        self
    }
}

impl BoostParamCommon {
    pub fn iteration(&self) -> Option<IterationNumber> {
        self.iteration
    }

    pub fn loss_func(&self) -> LossFunc {
        self.loss_func
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
    pub fn batch_interaction_way(&self) -> Option<BatchInteractionWay> {
        self.batch_interaction_way
    }
    pub fn eff_eps(&self) -> Option<EffEps> {
        self.eff_eps
    }
    pub fn is_monitor(&self) -> bool {
        self.is_monitor
    }
    pub fn monitor_nsnvs(&self) -> Option<&[usize]> {
        self.monitor_nsnvs.as_deref()
    }
    pub fn mhc_region(&self) -> Option<&(SnvId, SnvId)> {
        self.mhc_region.as_ref()
    }
    pub fn maf_threshold_logit_common(&self) -> Option<f64> {
        self.maf_threshold_logit_common
    }

    pub fn interaction_way(&self) -> Option<InteractionWay> {
        self.interaction_way
    }

    pub fn ld_criteria(&self) -> Option<LdCriteria> {
        self.ld_criteria
    }

    pub fn acc_metric(&self) -> Option<AccMetric> {
        self.acc_metic
    }

    pub fn set_iteration(self, iteration: Option<usize>, iteration_snv: Option<usize>) -> Self {
        let iteration_opt = match (iteration, iteration_snv) {
            (None, None) => None,
            _ => Some(IterationNumber::new(iteration, iteration_snv)),
        };
        Self {
            iteration: iteration_opt,
            ..self
        }
    }

    // pub fn set_iteration(self, iteration: usize) -> Self {
    //     Self {
    //         iteration: Some(IterationNumber::Iteration(iteration)),
    //         ..self
    //     }
    // }
    // pub fn set_iteration_snv(self, iteration: usize) -> Self {
    //     Self {
    //         iteration: Some(IterationNumber::Snv(iteration)),
    //         ..self
    //     }
    // }
    pub fn set_loss_func(self, loss_func: &str) -> Self {
        Self {
            loss_func: LossFunc::from_str(loss_func).unwrap(),
            ..self
        }
    }
    /// Some("None") -> None
    pub fn set_eps(self, eps: Option<&str>) -> Self {
        Self {
            eps: match eps {
                Some("None") => None,
                _ => eps.map(|x| Eps::from_str(x).unwrap()),
            },
            //eps: eps.map(|x| Eps::from_str(x).unwrap()),
            ..self
        }
    }
    pub fn set_eff_eps(self, eff_eps: Option<&str>) -> Self {
        Self {
            eff_eps: match eff_eps {
                Some("None") => None,
                _ => eff_eps.map(|x| EffEps::from_str(x).unwrap()),
            },
            ..self
        }
    }
    pub fn set_sample_weight_clip(self, sample_weight_clip: Option<&str>) -> Self {
        Self {
            sample_weight_clip: match sample_weight_clip {
                Some("None") => None,
                _ => sample_weight_clip.map(|x| SampleWeightClip::from_str(x).unwrap()),
            },
            ..self
        }
    }
    pub fn set_sample_weight_wls_clip(self, sample_weight_wls_clip: Option<&str>) -> Self {
        Self {
            sample_weight_wls_clip: match sample_weight_wls_clip {
                Some("None") => None,
                _ => sample_weight_wls_clip.map(|x| SampleWeightWlsClip::from_str(x).unwrap()),
            },
            ..self
        }
    }
    pub fn set_is_dom_rec(self, is_dom_rec: bool) -> Self {
        Self { is_dom_rec, ..self }
    }
    pub fn set_cov_way(self, cov_way: Option<&str>) -> Self {
        Self {
            cov_way: match cov_way {
                Some("None") => None,
                _ => cov_way.map(|x| CovWay::from_str(x).unwrap()),
            },
            ..self
        }
    }
    pub fn set_interaction_way(self, interaction_way: Option<&str>) -> Self {
        Self {
            interaction_way: match interaction_way {
                Some("None") => None,
                _ => interaction_way.map(|x| InteractionWay::from_str(x).unwrap()),
            },
            ..self
        }
    }

    pub fn set_batch_way(self, batch_way: Option<&str>) -> Self {
        Self {
            batch_way: match batch_way {
                Some("None") => None,
                _ => batch_way.map(|x| BatchWay::from_str(x).unwrap()),
            },
            ..self
        }
    }

    pub fn set_batch_interaction_way(self, batch_interaction_way: Option<&str>) -> Self {
        Self {
            batch_interaction_way: match batch_interaction_way {
                Some("None") => None,
                _ => batch_interaction_way.map(|x| BatchInteractionWay::from_str(x).unwrap()),
            },
            ..self
        }
    }

    pub fn set_is_monitor(self, is_monitor: bool) -> Self {
        Self { is_monitor, ..self }
    }

    pub fn set_monitor_nsnvs(self, monitor_nsnvs: Option<Vec<usize>>) -> Self {
        Self {
            monitor_nsnvs,
            ..self
        }
    }

    pub fn set_mhc_region(self, mhc_region: Option<&str>) -> Self {
        let mhc_region_parse = match mhc_region {
            Some(mhc_region) => {
                //let chrom = mhc_region.split(':').next().unwrap().to_string();
                let chrom_poss = mhc_region
                    .split(':')
                    .map(|x| x.to_string())
                    .collect::<Vec<String>>();
                let chrom = chrom_poss[0].clone();
                let poss = chrom_poss[1].clone();

                let mhc_region_parse = poss
                    .split('-')
                    .map(|x| SnvId::new("".to_string(), &chrom, x, "".to_string(), "".to_string()))
                    .collect::<Vec<SnvId>>();
                if mhc_region_parse.len() != 2 {
                    panic!("MHC region must be specified by two snv ids");
                }
                Some((mhc_region_parse[0].clone(), mhc_region_parse[1].clone()))
            }
            None => None,
        };

        Self {
            mhc_region: mhc_region_parse,
            ..self
        }
    }

    pub fn set_maf_threshold_logit_common(self, maf_threshold_logit_common: Option<f64>) -> Self {
        Self {
            maf_threshold_logit_common,
            ..self
        }
    }

    pub fn set_ld_criteria(self, ld_r2: Option<f64>, ld_radius: Option<usize>) -> Self {
        let ld_criteria = match (ld_r2, ld_radius) {
            (None, None) => None,
            _ => Some(LdCriteria::new(ld_r2, ld_radius)),
        };
        Self {
            ld_criteria,
            ..self
        }
    }

    pub fn set_acc_metric(self, acc_metric: Option<&str>) -> Self {
        Self {
            acc_metic: match acc_metric {
                Some("None") => None,
                _ => acc_metric.map(|x| AccMetric::from_str(x).unwrap()),
            },
            ..self
        }
    }
}

#[derive(Eq, PartialEq, Copy, Clone, Hash, Debug)]
pub enum BoostType {
    Ada, // original adaboost
    ConstAda,
    Logit, // LogitBoost
    LogitNoMissing,
    LogitAdd,          // Additive Logit
    LogitMhcNoMissing, // Nonadd for MHC and Add for others, MHC region is in BoostParamCommon
    //LogitMhcNoMissing(SnvId, SnvId), // Nonadd for MHC and Add for others, (start_snv, end_snv) for MHC region
    LogitAddInteraction, // Additive + interaction Logit
    // both Modelfree and ModelfreeMissing use ModelfreeClassifier
    //FreeModel, -> FreeModelNoMissing
    FreeModelMissing,
    LogitCommon, // Nonadd for common variants only
                 // AdaLogit
                 //Eta(f64), // arg: eta:f64
                 //EtaLogit(f64)
}

impl BoostType {
    pub fn is_freemodelmissing(self) -> bool {
        match self {
            BoostType::FreeModelMissing => true,
            BoostType::Ada
            | BoostType::ConstAda
            | BoostType::Logit
            | BoostType::LogitNoMissing
            | BoostType::LogitAdd
            | BoostType::LogitMhcNoMissing
            | BoostType::LogitAddInteraction
            | BoostType::LogitCommon => false,
        }
    }

    // use below is_type_logti()
    //pub fn is_logit(self) -> bool {
    //    match self {
    //        BoostType::Logit | BoostType::LogitNoMissing => true,
    //        BoostType::Ada
    //        | BoostType::ConstAda
    //        | BoostType::LogitAdd  // all right?
    //        | BoostType::FreeModelMissing => false,
    //    }
    //}

    // FIXME: remove
    // allow missing value in dataset
    //pub fn use_missing(self) -> bool {
    //    match self {
    //        BoostType::Ada
    //        | BoostType::ConstAda
    //        | BoostType::LogitAdd
    //        | BoostType::LogitNoMissing => false,
    //        BoostType::FreeModelMissing | BoostType::Logit => true,
    //    }
    //}

    // fill missing value in dataset
    pub fn fill_missing(self) -> bool {
        match self {
            BoostType::Ada
            | BoostType::ConstAda
            | BoostType::LogitAdd
            | BoostType::LogitNoMissing
            | BoostType::LogitMhcNoMissing
            | BoostType::LogitAddInteraction
            | BoostType::LogitCommon => true,
            BoostType::FreeModelMissing | BoostType::Logit => false,
        }
    }

    // TODO: better
    pub fn is_type_ada(self) -> bool {
        match self {
            BoostType::Ada | BoostType::ConstAda | BoostType::FreeModelMissing => true,
            BoostType::Logit
            | BoostType::LogitNoMissing
            | BoostType::LogitAdd
            | BoostType::LogitMhcNoMissing
            | BoostType::LogitAddInteraction
            | BoostType::LogitCommon => false,
        }
    }

    pub fn is_type_logit(self) -> bool {
        match self {
            BoostType::Logit
            | BoostType::LogitNoMissing
            | BoostType::LogitAdd
            | BoostType::LogitMhcNoMissing
            | BoostType::LogitAddInteraction
            | BoostType::LogitCommon => true,
            BoostType::Ada | BoostType::ConstAda | BoostType::FreeModelMissing => false,
        }
    }

    pub fn is_interaction(self) -> bool {
        match self {
            BoostType::LogitAddInteraction => true,
            _ => false,
        }
    }

    // for directory name of integrate
    pub fn genetic_model(self) -> String {
        match self {
            BoostType::Logit | BoostType::LogitNoMissing => "nonadd".to_string(),
            BoostType::LogitAdd => "add".to_string(),
            BoostType::LogitMhcNoMissing => "mhc_nonadd".to_string(),
            BoostType::LogitCommon => "common_nonadd".to_string(),
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
            "logitmhcnomissing" => Ok(BoostType::LogitMhcNoMissing),
            "logitaddinteraction" => Ok(BoostType::LogitAddInteraction),
            "logitcommon" => Ok(BoostType::LogitCommon),
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

            //if y < 0.0 || y<1.0 {
            if y < 1.0 {
                panic!("y should be 0.0 or positive.")
            }
            //if z < 0.0 || z<1.0{
            if z < 1.0 {
                panic!("z should be 0.0 or positive.")
            }
            //if u < 0.0 || u<1.0{
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
            if (u < 0.0) || (u > 1.0) {
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
            if (u < 0.0) || (u > 1.0) {
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
        // true: only apply adjusting in PGS not selecting on loss
        // false: apply adjusting in PGS and loss
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

//// deprecated, use BoostType instead
//#[derive(Debug, Copy, Clone)]
//pub enum BoostMethod {
//    Classic,      // usual boosting
//    Pruning(f64), // arg: clump_r2
//    Ss(f64),
//}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum IterationNumber {
    Iteration(usize),
    Snv(usize),
}

impl IterationNumber {
    pub fn new(iteration: Option<usize>, iteration_snv: Option<usize>) -> Self {
        match (iteration, iteration_snv) {
            (Some(x), None) => Self::Iteration(x),
            (None, Some(x)) => Self::Snv(x),
            (None, None) => panic!("Cannot specify neither iteration nor iteration_snv."),
            (Some(_), Some(_)) => panic!("Cannot specify both iteration and iteration_snv."),
        }
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
    EveryMaxIter(usize, usize, usize), // (batch_size, batch_max_iter_lr1, batch_max_iter_lr_others)
    Fix(usize, usize, usize),          // (batch_size, batch_iter_lr1, batch_iter_lr_others)
    LdFix(usize, usize, usize), // (batch_size, batch_iter_lr1, batch_iter_lr_others), batch_size: top # snvs which decreased loss ~ # snvs in LD
}

impl BatchWay {
    pub fn batch_size(self) -> usize {
        match self {
            Self::Every(x) => x,
            Self::EveryMaxIter(x, ..) | Self::Fix(x, ..) => x,
            Self::LdFix(x, ..) => x,
        }
    }
    pub fn use_comp_loss(self) -> bool {
        match self {
            Self::Every(..) | Self::EveryMaxIter(..) => true,
            Self::Fix(..) | Self::LdFix(..) => false,
        }
    }
    pub fn batch_max_iter(self, lr: f64) -> usize {
        match self {
            Self::Every(x) => x,
            Self::EveryMaxIter(_, x1, x2) | Self::Fix(_, x1, x2) | Self::LdFix(_, x1, x2) => {
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
        } else if str.starts_with("ldfix-") {
            let mut str_iter = str.split("-");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "ldfix");

            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<usize>().unwrap();
            let z = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            return Ok(BatchWay::LdFix(x, y, z));
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

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum BatchInteractionWay {
    // batchsize for interaction set, # interaction to filter from batch,  # interaction to filter from batch_all
    // batchsize is not used to judge renewing batch or not
    FixFilter(usize, usize, usize),
    // max # for initial candidates
    // initial_size, (same as FixFilter)
    //InitialMaxFixFilter(usize, usize, usize, usize),
    //
    // ld_width, (same as FixFilter)
    // To create initial candidates, prepare snv pairs of SNVs which has largest loss in ld_width bp
    // For the first iteration, the smallest {#interaction all} will be extract from the candidates.
    // For the next iteration, the selected SNVs will be added.
    InitialLdWidthMaxFixFilter(usize, usize, usize, usize),
    InitialLdWidthRandomFixFilter(usize, usize, usize, usize),
    InitialAllFixFilter(usize, usize, usize),
    // TODO: rewite
    // InitialAllFixFilter{
    //             batch_size: usize,
    //              filter_size_batch: usize,
    //             filter_size_all: usize,
    //}
}

impl BatchInteractionWay {
    pub fn is_filter_all(self) -> bool {
        match self {
            Self::FixFilter(..) => true,
            Self::InitialLdWidthMaxFixFilter(..) => true,
            Self::InitialLdWidthRandomFixFilter(..) => true,
            Self::InitialAllFixFilter(..) => true,
        }
    }

    pub fn is_initial(self) -> bool {
        match self {
            Self::FixFilter(..) => false,
            Self::InitialLdWidthMaxFixFilter(..) => true,
            Self::InitialLdWidthRandomFixFilter(..) => true,
            Self::InitialAllFixFilter(..) => true,
        }
    }

    pub fn batch_size(self) -> usize {
        match self {
            Self::FixFilter(x, ..) => x,
            Self::InitialLdWidthMaxFixFilter(_, x, _, _) => x,
            Self::InitialLdWidthRandomFixFilter(_, x, _, _) => x,
            Self::InitialAllFixFilter(x, _, _) => x,
        }
    }

    pub fn filter_size_batch(self) -> usize {
        match self {
            Self::FixFilter(_, x, _) => x,
            Self::InitialLdWidthMaxFixFilter(_, _, x, _) => x,
            Self::InitialLdWidthRandomFixFilter(_, _, x, _) => x,
            Self::InitialAllFixFilter(_, x, _) => x,
        }
    }

    pub fn filter_size_all(self) -> usize {
        match self {
            Self::FixFilter(.., x) => x,
            Self::InitialLdWidthMaxFixFilter(_, _, _, x) => x,
            Self::InitialLdWidthRandomFixFilter(_, _, _, x) => x,
            Self::InitialAllFixFilter(_, _, x) => x,
        }
    }

    pub fn ld_width(self) -> usize {
        match self {
            Self::FixFilter(..) => panic!("ld_width is not used for FixFilter."),
            Self::InitialLdWidthMaxFixFilter(x, ..) => x,
            Self::InitialLdWidthRandomFixFilter(x, ..) => x,
            Self::InitialAllFixFilter(..) => {
                panic!("ld_width is not used for InitialAllFixFilter.")
            }
        }
    }
}

impl FromStr for BatchInteractionWay {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        if str.starts_with("fixfilter-") {
            let mut str_iter = str.split("-");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "fixfilter");

            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<usize>().unwrap();
            let z = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            if x > z {
                panic!("batchsize should be <= filter_size_all.")
            }

            return Ok(BatchInteractionWay::FixFilter(x, y, z));
        } else if str.starts_with("initialldwidthmaxfixfilter-") {
            let mut str_iter = str.split("-");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "initialldwidthmaxfixfilter");

            let u = str_iter.next().unwrap().parse::<usize>().unwrap();
            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<usize>().unwrap();
            let z = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            if x > z {
                panic!("batchsize should be <= filter_size_all.")
            }

            return Ok(BatchInteractionWay::InitialLdWidthMaxFixFilter(u, x, y, z));
        } else if str.starts_with("initialldwidthrandomfixfilter-") {
            let mut str_iter = str.split("-");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "initialldwidthrandomfixfilter");

            let u = str_iter.next().unwrap().parse::<usize>().unwrap();
            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<usize>().unwrap();
            let z = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            if x > z {
                panic!("batchsize should be <= filter_size_all.")
            }

            return Ok(BatchInteractionWay::InitialLdWidthRandomFixFilter(
                u, x, y, z,
            ));
        } else if str.starts_with("initialallfixfilter-") {
            let mut str_iter = str.split("-");

            let str0 = str_iter.next().unwrap();
            assert_eq!(str0, "initialallfixfilter");

            //let u = str_iter.next().unwrap().parse::<usize>().unwrap();
            let x = str_iter.next().unwrap().parse::<usize>().unwrap();
            let y = str_iter.next().unwrap().parse::<usize>().unwrap();
            let z = str_iter.next().unwrap().parse::<usize>().unwrap();
            assert_eq!(str_iter.next(), None);

            if x > z {
                panic!("batchsize should be <= filter_size_all.")
            }

            return Ok(BatchInteractionWay::InitialAllFixFilter(x, y, z));
        } else {
            panic!("Unknown BatchInteractionWay: {}", str);
        }
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum InteractionWay {
    //OneBySelected,
    //OneBySelectedLd, // use ld_radius
    SelectedBySelected,
    SelectedBySelectedLd,
    //OneBySelectedForce, //force selecte one interaction term after selecting one SNV
    //Every(usize),
}

impl InteractionWay {}

impl FromStr for InteractionWay {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        //if str == "onebyselected" {
        //    return Ok(InteractionWay::OneBySelected);
        //} else if str == "onebyselectedld" {
        //    return Ok(InteractionWay::OneBySelectedLd);
        if str == "selectedbyselected" {
            return Ok(InteractionWay::SelectedBySelected);
        } else if str == "selectedbyselectedld" {
            return Ok(InteractionWay::SelectedBySelectedLd);
        }
        Err(format!("Unknown InteractionWay: {}", str))
    }
}

#[derive(PartialEq, Copy, Clone, Debug)]
pub enum AccMetric {
    CovAdjustedPseudoR2,
    AUC,
}

impl FromStr for AccMetric {
    type Err = String;
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        if str == "cov-adjusted-pseudo-r2" {
            return Ok(AccMetric::CovAdjustedPseudoR2);
        } else if str == "auc" {
            return Ok(AccMetric::AUC);
        }
        Err(format!("Unknown AccMetric: {}", str))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boosting_params() {
        let lrs = vec![1.0, 0.1];
        let boost_params = BoostParamLrs::default().set_learning_rates(lrs);

        let mut boost_iter = boost_params.into_iter();

        assert_eq!(boost_iter.next().unwrap().learning_rate(), 1.0);
        assert_eq!(boost_iter.next().unwrap().learning_rate(), 0.1);
    }
}
