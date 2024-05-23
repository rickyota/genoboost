use super::coef::Coef;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub struct Model {
    coef: Coef,
    // put threshold in model_type?
    threshold: Option<f64>,
    // for now, use Coef
    //model_type: ModelType,
}

impl Model {
    pub fn new(coef: Coef, threshold: Option<f64>) -> Self {
        Self { coef, threshold }
    }
    pub fn new_linear(alpha: f64) -> Model {
        Self::new(Coef::Linear(alpha), None)
        //Model {
        //    coef: Coef::Linear(alpha),
        //    threshold: None,
        //}
    }
    pub fn new_coef(coef: Coef) -> Model {
        Self::new(coef, None)
        //Model {
        //    coef,
        //    threshold: None,
        //}
    }

    pub fn coef(self) -> Coef {
        self.coef
    }
    pub fn threshold(self) -> Option<f64> {
        self.threshold
    }
    pub fn model_string_for_wgt(self) -> String {
        self.coef().model_string_for_wgt()

        //match self.coef().model_type() {
        //    ModelType::Linear => "LINEAR".to_owned(),
        //    ModelType::Single => "SINGLE".to_owned(),
        //    ModelType::Binary => "BINARY".to_owned(),
        //    ModelType::Free => "FREE".to_owned(),
        //    ModelType::NaN => panic!("Coef is NaN"),
        //}
    }

    pub fn set_coef(&mut self, coef: Coef) {
        self.coef = coef;
    }

    pub fn threshold_string(self) -> String {
        match self.threshold() {
            Some(x) => format!("{:.7}", x),
            None => "NaN".to_owned(),
        }
    }

    pub fn to_string_hash(self) -> HashMap<String, String> {
        self.coef().to_string_hash()
    }

    /*
    pub fn to_string_vec(self) -> Vec<String> {
        self.coef().to_string_vec()
    }
     */
    /*
    pub fn const_string(&self) -> String {
        self.coef().const_to_string()
        // TMP
        //format!("{:.7}", 0.0)
        //format!("{:.7}", self.coef.0)
    }
    pub fn alpha_string(&self) -> String {
        self.coef().alpha_to_string()
        // TMP
        //format!("{:.7}", 0.0)
        //format!("{:.7}", self.coef.1)
    }
     */
}

//pub enum ModelType {
//    Const,
//    Linear,
//    LinearConst,
//    Free2(f64), //threshold
//    Free3,
//    Free4,
//    Interaction, // Coef::LienarCoef
//}
