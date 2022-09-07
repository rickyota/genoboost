use crate::{CovId, SnvId, Var};

#[derive(Debug, Clone)]
pub enum WgtKind {
    // remove indexs later, since dataset should implement finding index
    // index_inherit(=si) and index(=mi)
    Snv(SnvId, Option<usize>, Option<usize>),
    Cov(CovId),
    //Snv(SnvWgt),
    //Cov(CovWgt), // rename to Var later, Const is in Cov
}

impl WgtKind {
    pub fn index_predict_snv(&self) -> (usize, usize) {
        if let WgtKind::Snv(_, si, mi) = self {
            return (si.unwrap(), mi.unwrap());
        }
        panic!("This wgtkind is not snv.");
    }

    pub fn index_inherit_snv(&self) -> Option<usize> {
        if let WgtKind::Snv(_, si, _) = self {
            return *si;
        }
        panic!("This wgtkind is not snv.");
    }
    pub fn index_snv(&self) -> Option<usize> {
        if let WgtKind::Snv(_, _, mi) = self {
            return *mi;
        }
        panic!("This wgtkind is not snv.");
    }

    pub fn snv_index(&self) -> &SnvId {
        if let WgtKind::Snv(snv_wgt, _, _) = self {
            return snv_wgt;
        }
        panic!("This wgtkind is not snv.");
    }

    pub fn is_snv(&self) -> bool {
        matches!(self, WgtKind::Snv(_, _, _))
    }

    /// including const
    pub fn is_cov(&self) -> bool {
        matches!(self, WgtKind::Cov(_))
    }

    // this is necessary when user knows only cov can happen.
    // so must not Option<CovId>
    pub fn cov(&self) -> &CovId {
        match self {
            WgtKind::Cov(x) => x,
            _ => panic!("Not cov."),
        }
    }

    /*
    pub fn snv_wgt(&self) -> &SnvWgt {
        if let WgtKind::Snv(snv_wgt) = self {
            return snv_wgt;
        }
        panic!("This wgt is not Snv Wgt.");
    }

    pub fn snv_wgt_mut(&mut self) -> &mut SnvWgt {
        if let WgtKind::Snv(ref mut snv_wgt) = self {
            return snv_wgt;
        }
        panic!("This wgt is not Snv Wgt.");
    }
     */
}

/*
#[derive(Debug, Clone)]
pub struct ConstWgt {
    model: LinearModel,
}

impl ConstWgt {
    pub fn construct_const_wgt(model: LinearModel) -> ConstWgt {
        ConstWgt { model }
    }
    pub fn get_model(&self) -> &LinearModel {
        &self.model
    }
    pub fn get_model_mut(&mut self) -> &mut LinearModel {
        &mut self.model
    }
}
*/

#[derive(Debug, Clone)]
pub struct SnvWgt {
    snv_index: SnvId,
    index_inherit: Option<usize>, // := si, 0(dom) or 1(rec)
    index: Option<usize>,         // := mi index
                                  /*
                                  index_predict: (usize, usize),
                                  index: Option<usize>, // for score
                                   */
}

// index_predict: (index_inherit, index)
// inherit=0: dominant or 1: recessive
impl SnvWgt {
    /*
    pub fn construct(snv: Snv, indexs_predict: (usize, usize)) -> SnvWgt {
        SnvWgt {
            snv_index: snv.snv_index().clone(),
            index_inherit: Some(indexs_predict.0),
            index: Some(indexs_predict.1),
        }
    }
     */
    pub fn construct_from_snv_index(snv_index: SnvId, indexs_predict: (usize, usize)) -> SnvWgt {
        SnvWgt {
            snv_index,
            index_inherit: Some(indexs_predict.0),
            index: Some(indexs_predict.1),
        }
    }
    // set index only
    pub fn construct_ss(snv: SnvId, index: usize) -> SnvWgt {
        SnvWgt {
            snv_index: snv.clone(),
            index_inherit: None,
            index: Some(index),
        }
    }
    /// set index later
    pub fn construct_score(snv_index: SnvId) -> SnvWgt {
        SnvWgt {
            snv_index: snv_index.clone(),
            index_inherit: None,
            index: None,
        }
    }
    /*
    pub fn construct_snv_wgt(snv: Snv, index: (usize, usize), model: BinaryModel) -> SnvWgt {
        SnvWgt { snv, index, model }
    }
    */
    pub fn snv_index(&self) -> &SnvId {
        &self.snv_index
    }
    pub fn index_predict(&self) -> (usize, usize) {
        (self.index_inherit.unwrap(), self.index.unwrap())
        //self.index_predict
    }

    pub fn index(&self) -> Option<usize> {
        self.index
    }
    pub fn set_index(&mut self, index: Option<usize>) {
        self.index = index;
    }
    /*
    pub fn get_model(&self) -> &BinaryModel {
        &self.model
    }
    pub fn get_model_mut(&mut self) -> &mut BinaryModel {
        &mut self.model
    }
    */
}

#[derive(Debug, Clone)]
pub struct CovWgt {
    var: Var,
    //index: usize, // index is used to match mistakes and SnvWgt. Cov does not need this since values are in Var.
    //model: Model,
}

impl CovWgt {
    pub fn construct(cov: Var) -> CovWgt {
        CovWgt { var: cov }
    }
    /*
    pub fn construct_empty() -> CovWgt {
        let cov = Var::construct_empty();
        //pub fn construct_var(kind: VarKind, name: String) -> Var {
        CovWgt { var: cov }
    }
     */
    /*
    pub fn construct_cov_wgt(cov: Var, index: usize, model: Model) -> CovWgt {
        CovWgt {
            var: cov,
            index,
            model,
        }
    }
    */
    pub fn var(&self) -> &Var {
        &self.var
    }

    pub fn get_var_mut(&mut self) -> &mut Var {
        &mut self.var
    }

    /*
    pub fn set_var(&mut self, var: &Var) {
        self.var = (*var).clone();
    }
    */
    /*
    pub fn get_index(&self) -> usize {
        self.index
    }
    */
    /*
    pub fn get_model(&self) -> &Model {
        &self.model
    }
    pub fn get_model_mut(&mut self) -> &mut Model {
        &mut self.model
    }
    */
}
