pub mod coef;
pub mod io;
pub mod model;
pub mod wgt_kind;

pub use coef::{Coef, ModelType};
pub use model::Model;
pub use wgt_kind::{CovWgt, SnvWgt, WgtKind};

use crate::{CovId, SnvId};
//pub use wgt_kind::{ConstWgt, CovWgt, SnvWgt, WgtKind};

// now only for wgtboost
pub trait WgtTrait {
    fn wgt(&self) -> &Wgt;
    fn wgt_mut(&mut self) -> &mut Wgt;
    fn kind(&self) -> &WgtKind {
        self.wgt().kind()
    }
    fn model(&self) -> &Model {
        self.wgt().model()
    }
    //fn add_snv_index(&self, mi: usize) -> Wgt;

    fn set_snv_index(&mut self, mi: Option<usize>);
}

#[derive(Debug, Clone)]
pub struct Wgt {
    kind: WgtKind,
    model: Model,
}

impl WgtTrait for Wgt {
    fn wgt(&self) -> &Wgt {
        self.wgt()
    }
    fn wgt_mut(&mut self) -> &mut Wgt {
        self.wgt_mut()
    }

    fn set_snv_index(&mut self, mi: Option<usize>) {
        let kind = self.wgt_mut().kind_mut();
        // TODO: move to WgtKind::
        if let WgtKind::Snv(_, _, ref mut index) = kind {
            *index = mi;
        }
    }
}



impl Wgt {
    pub fn construct_wgt(kind: WgtKind, model: Model) -> Wgt {
        let wgt = Wgt { kind, model };
        wgt
    }

    pub fn construct_snv(snv_id: SnvId, coef: Coef) -> Wgt {
        Wgt {
            kind: WgtKind::Snv(snv_id, None, None),
            model: Model::new_coef(coef),
        }
    }

    pub fn construct_cov(cov_id: CovId, coef: Coef) -> Wgt {
        Wgt {
            kind: WgtKind::Cov(cov_id),
            model: Model::new_coef(coef),
        }
    }
    pub fn construct_cov_name(name: String, coef: Coef) -> Wgt {
        Wgt {
            kind: WgtKind::Cov(CovId::new_cov(name)),
            model: Model::new_coef(coef),
        }
    }
    pub fn construct_const(model: Model) -> Wgt {
        Wgt {
            kind: WgtKind::Cov(CovId::new_const()),
            model,
        }
    }
    pub fn construct_const_linear(coef: Coef) -> Wgt {
        Wgt {
            kind: WgtKind::Cov(CovId::new_const()),
            model: Model::new_coef(coef),
        }
    }
    pub fn construct_const_threshold(threshold: f64) -> Wgt {
        Wgt {
            kind: WgtKind::Cov(CovId::new_const()),
            model: Model::new(Coef::NaN, Some(threshold)),
        }
    }

    pub fn construct_snv_index(
        snv_id: SnvId,
        threshold: f64,
        indexs_predict: (usize, usize),
    ) -> Wgt {
        Wgt {
            kind: WgtKind::Snv(snv_id, Some(indexs_predict.0), Some(indexs_predict.1)),
            model: Model::new(Coef::NaN, Some(threshold)),
        }
    }

    pub fn construct_snv_index_freemodel(snv_id: SnvId, index_predict: usize) -> Wgt {
        Wgt {
            kind: WgtKind::Snv(snv_id, None, Some(index_predict)),
            model: Model::new(Coef::NaN, None),
        }
    }

    pub fn set_snv_index(&mut self, mi: Option<usize>) {
        let kind = self.kind_mut();
        // all right? -> I think so
        if let WgtKind::Snv(_, ref mut index, _) = kind {
            *index = mi;
        }
    }

    /*
    pub fn add_snv_index(&self, mi: usize) -> Wgt {
        //let mut wgt = self.clone();
        Wgt {
            kind: WgtKind::Snv(
                self.kind().snv_index().to_owned(),
                self.kind().index_inherit_snv(),
                Some(mi),
            ),
            model: self.model().to_owned(),
        }
    }
     */

    pub fn wgt(&self) -> &Wgt {
        self
    }

    pub fn wgt_mut(&mut self) -> &mut Wgt {
        self
    }


    pub fn kind(&self) -> &WgtKind {
        &self.kind
    }

    pub fn kind_mut(&mut self) -> &mut WgtKind {
        &mut self.kind
    }

    /*
    pub fn snv_wgt(&self) -> &SnvWgt {
        self.kind().snv_wgt()
    }

    pub fn snv_wgt_mut(&mut self) -> &mut SnvWgt {
        self.kind_mut().snv_wgt_mut()
    }
     */

    /*
    pub fn kind_snv(&self) -> &SnvWgt {
        if let WgtKind::Snv(snv_wgt) = self.kind() {
            return snv_wgt;
        }
        panic!("This wgt is not Snv kind.");
    }

    pub fn kind_snv_mut(&mut self) -> &mut SnvWgt {
        if let WgtKind::Snv(ref mut snv_wgt) = self.kind_mut() {
            return snv_wgt;
        }
        panic!("This wgt is not Snv kind.");
    }
    */

    pub fn is_snv(&self) -> bool {
        self.kind().is_snv()
    }

    pub fn is_cov(&self) -> bool {
        self.kind().is_cov()
    }

    /*
    pub fn is_snv(&self) -> bool {
        match &self.kind() {
            WgtKind::Cov(_) => {
                return false;
            }
            WgtKind::Snv(_) => {
                return true;
            }
        };
    }

    /// including const
    pub fn is_cov(&self) -> bool {
        match &self.kind() {
            WgtKind::Cov(_) => {
                return true;
            }
            WgtKind::Snv(_) => {
                return false;
            }
        };
    }
    */

    pub fn model(&self) -> &Model {
        &self.model
    }
    pub fn model_mut(&mut self) -> &mut Model {
        &mut self.model
    }

    /*
    pub fn model_binary(&self) -> &BinaryModel {
        self.model().model_binary()
    }

    pub fn model_binary_mut(&mut self) -> &mut BinaryModel {
        self.model_mut().model_binary_mut()
    }

    pub fn model_linear(&self) -> &LinearModel {
        self.model().model_linear()
    }

    pub fn model_linear_mut(&mut self) -> &mut LinearModel {
        self.model_mut().model_linear_mut()
    }
     */
}
