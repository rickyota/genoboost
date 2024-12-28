use crate::{CovId, SnvId};
//use crate::{CovId, SnvId, Var};

#[derive(Debug, Clone)]
pub enum WgtKind {
    // index_inherit(=si) and index(=mi)
    //Snv(SnvId, Option<usize>, Option<usize>),
    Snv(SnvWgt),
    //  snv1, index of snv1, snv2, index of snv2
    SnvInteraction(SnvInteractionWgt),
    //SnvInteraction(SnvId, Option<usize>, SnvId, Option<usize>),
    Cov(CovId),
    // TODO: add snv wgt type additive, dominant etc. here
    //Snv(SnvWgt),
    //Cov(CovWgt), // rename to Var later, Const is in Cov
}

impl WgtKind {
    pub fn new_snv(snv_wgt: SnvWgt) -> Self {
        Self::Snv(snv_wgt)
    }

    pub fn new_snv_interaction(snv_inter_wgt: SnvInteractionWgt) -> Self {
        Self::SnvInteraction(snv_inter_wgt)
    }

    pub fn snv_wgt(&self) -> &SnvWgt {
        if let WgtKind::Snv(snv_wgt) = self {
            return snv_wgt;
        }
        panic!("This wgtkind is not snv.");
    }

    pub fn index_predict_snv(&self) -> (usize, usize) {
        if let WgtKind::Snv(snv_wgt) = self {
            return snv_wgt.index_predict();
        }
        //if let WgtKind::Snv(_, si, mi) = self {
        //    return (si.unwrap(), mi.unwrap());
        //}
        panic!("This wgtkind is not snv.");
    }

    //pub fn index_inherit_snv(&self) -> Option<usize> {
    //    if let WgtKind::Snv(_, si, _) = self {
    //        return *si;
    //    }
    //    panic!("This wgtkind is not snv.");
    //}
    pub fn index_snv(&self) -> Option<usize> {
        if let WgtKind::Snv(snv_wgt) = self {
            return snv_wgt.index();
        }
        //if let WgtKind::Snv(_, _, mi) = self {
        //    return *mi;
        //}
        panic!("This wgtkind is not snv.");
    }

    pub fn snv_id(&self) -> &SnvId {
        if let WgtKind::Snv(snv_wgt) = self {
            return snv_wgt.snv_id();
        }
        //if let WgtKind::Snv(snv_wgt, _, _) = self {
        //    return snv_wgt;
        //}
        panic!("This wgtkind is not snv.");
    }

    pub fn snv_id_interaction(&self) -> (&SnvId, &SnvId) {
        if let WgtKind::SnvInteraction(snv_inter_wgt) = self {
            return snv_inter_wgt.snv_ids();
            //return (snv_wgt_1, snv_wgt_2);
        }
        panic!("This wgtkind is not snv interaction.");
    }

    // is_snv() is confusing
    pub fn is_snv_single(&self) -> bool {
        matches!(self, WgtKind::Snv(..))
    }

    pub fn is_snv_interaction(&self) -> bool {
        matches!(self, WgtKind::SnvInteraction(..))
    }

    pub fn is_snv_single_or_interaction(&self) -> bool {
        match self {
            WgtKind::Snv(..) => true,
            WgtKind::SnvInteraction(..) => true,
            _ => false,
        }
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

    // use set_snv_index_check()
    //pub fn set_snv_index(&mut self, mi: Option<usize>) {
    //    if let WgtKind::Snv(snv_wgt) = self {
    //        snv_wgt.set_index(mi);
    //        //*index = mi;
    //        return;
    //    }
    //    //if let WgtKind::Snv(_, _, ref mut index) = self {
    //    //    *index = mi;
    //    //    return;
    //    //}
    //    panic!("This wgtkind is not single snv.");
    //}

    // check consistecy if already set
    pub fn set_snv_index_check(&mut self, mi: Option<usize>) {
        if let WgtKind::Snv(snv_wgt) = self {
            snv_wgt.set_index_check_for_wgt(mi);
            return;
        }
        panic!("This wgtkind is not single snv.");
    }

    //pub fn set_snv_index_interaction(&mut self, mi_1: Option<usize>, mi_2: Option<usize>) {
    //    if let WgtKind::SnvInteraction(snv_inter_wgt) = self {
    //        snv_inter_wgt.set_index(mi_1, mi_2);
    //        return;
    //    }
    //    panic!("This wgtkind is not snv interaction.");
    //}

    pub fn set_snv_index_interaction_check(&mut self, mi_1: Option<usize>, mi_2: Option<usize>) {
        if let WgtKind::SnvInteraction(snv_inter_wgt) = self {
            snv_inter_wgt.set_index_check_for_wgt(mi_1, mi_2);
            return;
        }
        panic!("This wgtkind is not snv interaction.");
    }

    pub fn set_a1_freq(&mut self, a1_freq: Option<f64>) {
        if let WgtKind::Snv(snv_wgt) = self {
            snv_wgt.set_freq(a1_freq);
            return;
        }
        panic!("This wgtkind is not single snv.");
    }

    //pub fn set_snv_index_interaction(&mut self, mi_1: Option<usize>, mi_2: Option<usize>) {
}

/*
// moved to Cov
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

// TODO
#[derive(Debug, Clone)]
pub struct SnvWgt {
    snv_id: SnvId,
    index_inherit: Option<usize>, // := si, 0(dom) or 1(rec)
    index: Option<usize>,         // := mi index
    // TODO: rename to a1_frq
    maf: Option<f64>, // for missing values required for score
}

// index_predict: (index_inherit, index)
// inherit=0: dominant or 1: recessive
impl SnvWgt {
    pub fn snv_id(&self) -> &SnvId {
        &self.snv_id
    }

    pub fn index_predict(&self) -> (usize, usize) {
        (self.index_inherit.unwrap(), self.index.unwrap())
    }

    pub fn index(&self) -> Option<usize> {
        self.index
    }

    pub fn maf(&self) -> Option<f64> {
        self.maf
    }

    pub fn new(
        snv_id: SnvId,
        index_inherit: Option<usize>,
        index: Option<usize>,
        maf: Option<f64>,
    ) -> Self {
        Self {
            snv_id,
            index_inherit,
            index,
            maf,
        }
    }

    pub fn new_snv_id(snv_id: SnvId) -> Self {
        Self::new(snv_id, None, None, None)
    }

    pub fn new_score(snv_id: SnvId, maf: Option<f64>) -> Self {
        Self::new(snv_id, None, None, maf)
    }

    pub fn set_index(&mut self, index: Option<usize>) {
        self.index = index;
    }

    pub fn set_index_check_for_wgt(&mut self, index: Option<usize>) {
        check_set_snv_index(index, self.index());

        self.set_index(index);
    }

    pub fn set_freq(&mut self, a1_freq: Option<f64>) {
        self.maf = a1_freq;
    }
}

fn check_set_snv_index(index_new: Option<usize>, index_old: Option<usize>) {
    //fn check_set_snv_index(index_old: Option<usize>, index_new: Option<usize>) {
    match (index_new, index_old) {
        // Cannot overwrite to another index
        (Some(mi_new), Some(mi_old)) => {
            if mi_new != mi_old {
                panic!("Inconsistent mi. mi_now: {}, mi: {}", mi_old, mi_new);
            }
        }
        // Cannot overwrite to None
        (None, Some(mi_old)) => {
            panic!("Inconsistent mi. mi_now: {}, mi: {:?}", mi_old, index_new);
        }
        _ => {}
    }
}

#[derive(Debug, Clone)]
pub struct SnvInteractionWgt {
    snv_id_1: SnvId,
    snv_id_2: SnvId,
    index_1: Option<usize>, // := mi index
    index_2: Option<usize>, // := mi index
    maf_1: Option<f64>,     // for missing values required for score
    maf_2: Option<f64>,     // for missing values required for score
}

// index_predict: (index_inherit, index)
// inherit=0: dominant or 1: recessive
impl SnvInteractionWgt {
    pub fn snv_id_1(&self) -> &SnvId {
        &self.snv_id_1
    }

    pub fn snv_id_2(&self) -> &SnvId {
        &self.snv_id_2
    }

    pub fn snv_ids(&self) -> (&SnvId, &SnvId) {
        (self.snv_id_1(), &self.snv_id_2())
    }

    pub fn index_1(&self) -> Option<usize> {
        self.index_1
    }

    pub fn index_2(&self) -> Option<usize> {
        self.index_2
    }

    pub fn indexs(&self) -> (Option<usize>, Option<usize>) {
        (self.index_1(), self.index_2())
    }

    pub fn maf_1(&self) -> Option<f64> {
        self.maf_1
    }

    pub fn maf_2(&self) -> Option<f64> {
        self.maf_2
    }

    pub fn mafs(&self) -> (Option<f64>, Option<f64>) {
        (self.maf_1(), self.maf_2())
    }

    pub fn new(
        snv_id_1: SnvId,
        snv_id_2: SnvId,
        index_1: Option<usize>,
        index_2: Option<usize>,
        maf_1: Option<f64>,
        maf_2: Option<f64>,
    ) -> Self {
        Self {
            snv_id_1,
            snv_id_2,
            index_1,
            index_2,
            maf_1,
            maf_2,
        }
    }

    pub fn new_snv_id(snv_id_1: SnvId, snv_id_2: SnvId) -> Self {
        Self::new(snv_id_1, snv_id_2, None, None, None, None)
    }

    pub fn new_score(
        snv_id_1: SnvId,
        snv_id_2: SnvId,
        maf_1: Option<f64>,
        maf_2: Option<f64>,
    ) -> Self {
        Self::new(snv_id_1, snv_id_2, None, None, maf_1, maf_2)
    }

    pub fn set_index(&mut self, index_1: Option<usize>, index_2: Option<usize>) {
        self.index_1 = index_1;
        self.index_2 = index_2;
    }

    pub fn set_index_check_for_wgt(&mut self, index_1: Option<usize>, index_2: Option<usize>) {
        check_set_snv_index(index_1, self.index_1());
        check_set_snv_index(index_2, self.index_2());

        self.set_index(index_1, index_2);
        //self.set_index(index)
    }
}

/*
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
 */

#[cfg(test)]
mod tests {
    use super::*;

    // None->Some(1)
    #[test]
    fn test_wgt_kind_set_snv_index() {
        let snv_wgt = SnvWgt::new(SnvId::default(), None, None, None);
        let mut wgt_kind = WgtKind::new_snv(snv_wgt);
        wgt_kind.set_snv_index_check(Some(1));
        assert_eq!(wgt_kind.index_snv(), Some(1));
    }

    // None <- None
    // Some(1) <- None
    #[test]
    fn test_check_set_snv_index() {
        check_set_snv_index(None, None);
        check_set_snv_index(Some(1), None);
    }

    // None <- Some(1)
    #[should_panic]
    #[test]
    fn test_check_set_snv_index_2_panic() {
        check_set_snv_index(None, Some(1));
    }

    // Some(2) <- Some(1)
    #[should_panic]
    #[test]
    fn test_check_set_snv_index_3_panic() {
        check_set_snv_index(Some(2), Some(1));
    }
}
