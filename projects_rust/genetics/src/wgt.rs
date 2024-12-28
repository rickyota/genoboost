pub mod coef;
pub mod io;
pub mod model;
pub mod wgt_kind;

pub use coef::Coef;
pub use model::Model;
pub use wgt_kind::{SnvInteractionWgt, SnvWgt, WgtKind};
//pub use wgt_kind::{CovWgt, SnvWgt, WgtKind};

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

    //fn set_snv_index(&mut self, mi: Option<usize>);
    fn set_snv_index_check(&mut self, mi: Option<usize>) {
        let kind = self.wgt_mut().kind_mut();
        kind.set_snv_index_check(mi);
        ////  moved to WgtKind::
        //if let WgtKind::Snv(_, _, ref mut index) = kind {
        //    *index = mi;
        //}
    }

    fn set_freq(&mut self, freq: Option<f64>) {
        let kind = self.wgt_mut().kind_mut();
        kind.set_a1_freq(freq);
    }

    //fn set_snv_index(&mut self, mi: Option<usize>) {
    //    let kind = self.wgt_mut().kind_mut();
    //    kind.set_snv_index(mi);
    //    ////  moved to WgtKind::
    //    //if let WgtKind::Snv(_, _, ref mut index) = kind {
    //    //    *index = mi;
    //    //}
    //}

    fn set_snv_index_interaction_check(&mut self, mi_1: Option<usize>, mi_2: Option<usize>) {
        let kind = self.wgt_mut().kind_mut();
        kind.set_snv_index_interaction_check(mi_1, mi_2);
    }
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

    //fn set_snv_index(&mut self, mi: Option<usize>) {
    //    let kind = self.wgt_mut().kind_mut();
    //    kind.set_snv_index(mi);
    //    //// TODO: move to WgtKind::
    //    //if let WgtKind::Snv(_, _, ref mut index) = kind {
    //    //    *index = mi;
    //    //}
    //}
}

impl Wgt {
    pub fn new(kind: WgtKind, model: Model) -> Self {
        let wgt = Self { kind, model };
        wgt
    }

    //pub fn new_snv(snv_id: SnvId, coef: Coef) -> Self {
    //    Self::new(WgtKind::Snv(snv_id, None, None), Model::new_coef(coef))
    //}

    pub fn new_cov(cov_id: CovId, coef: Coef) -> Self {
        Self::new(WgtKind::Cov(cov_id), Model::new_coef(coef))
        //Wgt {
        //    kind: WgtKind::Cov(cov_id),
        //    model: Model::new_coef(coef),
        //}
    }
    pub fn new_cov_name(name: String, coef: Coef) -> Self {
        Self::new(WgtKind::Cov(CovId::new_cov(name)), Model::new_coef(coef))

        //Wgt {
        //    kind: WgtKind::Cov(CovId::new_cov(name)),
        //    model: Model::new_coef(coef),
        //}
    }
    pub fn new_const(model: Model) -> Self {
        Self::new(WgtKind::Cov(CovId::new_const()), model)
        //Wgt {
        //    kind: WgtKind::Cov(CovId::new_const()),
        //    model,
        //}
    }
    pub fn new_const_linear(coef: Coef) -> Self {
        Self::new(WgtKind::Cov(CovId::new_const()), Model::new_coef(coef))
        //Wgt {
        //    kind: WgtKind::Cov(CovId::new_const()),
        //    model: Model::new_coef(coef),
        //}
    }
    pub fn new_const_threshold(threshold: f64) -> Self {
        Self::new(
            WgtKind::Cov(CovId::new_const()),
            Model::new(Coef::NaN, Some(threshold)),
        )
        //Wgt {
        //    kind: WgtKind::Cov(CovId::new_const()),
        //    model: Model::new(Coef::NaN, Some(threshold)),
        //}
    }

    pub fn new_snv_id(
        snv_id: SnvId,
        threshold: f64,
        indexs_predict: (usize, usize),
        maf: f64,
    ) -> Wgt {
        let snv_wgt = SnvWgt::new(
            snv_id,
            Some(indexs_predict.0),
            Some(indexs_predict.1),
            Some(maf),
        );
        Self::new(
            WgtKind::new_snv(snv_wgt),
            //WgtKind::Snv(snv_id, Some(indexs_predict.0), Some(indexs_predict.1)),
            Model::new(Coef::NaN, Some(threshold)),
        )
    }
    //pub fn new_snv_id(snv_id: SnvId, threshold: f64, indexs_predict: (usize, usize)) -> Wgt {
    //    Self::new(
    //        WgtKind::Snv(snv_id, Some(indexs_predict.0), Some(indexs_predict.1)),
    //        Model::new(Coef::NaN, Some(threshold)),
    //    )
    //}

    pub fn new_snv_id_freemodel(snv_id: SnvId, snv_index: usize, maf: f64) -> Wgt {
        let snv_wgt = SnvWgt::new(
            snv_id,
            None,
            //Some(indexs_predict.0),
            Some(snv_index),
            Some(maf),
        );

        Self::new(
            WgtKind::new_snv(snv_wgt),
            //WgtKind::Snv(snv_id, None, Some(index_predict)),
            Model::new(Coef::NaN, None),
        )
    }

    pub fn new_snv_id_interaction(
        snv_id_1: SnvId,
        snv_id_2: SnvId,
        index_predict_1: usize,
        index_predict_2: usize,
        maf_1: f64,
        maf_2: f64,
    ) -> Wgt {
        let snv_inter_wgt = SnvInteractionWgt::new(
            snv_id_1,
            snv_id_2,
            Some(index_predict_1),
            Some(index_predict_2),
            Some(maf_1),
            Some(maf_2),
        );
        Self::new(
            WgtKind::new_snv_interaction(snv_inter_wgt),
            //WgtKind::SnvInteraction(
            //    snv_id_1,
            //    Some(index_predict_1),
            //    snv_id_2,
            //    Some(index_predict_2),
            //),
            Model::new(Coef::NaN, None),
        )

        //Wgt {
        //    kind: WgtKind::Snv(snv_id, None, Some(index_predict)),
        //    model: Model::new(Coef::NaN, None),
        //}
    }

    // in WgtTrait
    //pub fn set_snv_index(&mut self, mi: Option<usize>) {
    //    let kind = self.kind_mut();
    //    // all right? -> I think so
    //    if let WgtKind::Snv(_, ref mut index, _) = kind {
    //        *index = mi;
    //    }
    //}

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

    // TODO: revert name to is_snv() ?
    pub fn is_snv_single(&self) -> bool {
        self.kind().is_snv_single()
    }

    pub fn is_snv_single_or_interaction(&self) -> bool {
        self.kind().is_snv_single_or_interaction()
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
