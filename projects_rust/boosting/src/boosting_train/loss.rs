// pub for bench
pub mod calc;

use std::collections::HashSet;

//use super::coefficient;
//use super::compute_pred;
use super::epsilon;
use super::sample_weight::SampleWeight;
use super::table;
use super::BoostType;
use super::WgtBoost;
use crate::BoostParam;
use genetics::genot::prelude::*;
use genetics::genot::BaseGenot;
use genetics::Snvs;
//use genetics::vec::sort_float;
use genetics::vec;
use genetics::wgt::Wgt;
use genetics::Dataset;
use std::io::{BufWriter, Write};

#[derive(Debug, Clone)]
pub enum LossStruct {
    // loss, m, s=2
    // loss[mi, si]
    // (loss[0,0], loss[0,1], loss[1,0],loss[1,1]...)
    // TODO: create two Vec<f64>
    ConstAda(Vec<f64>, usize, usize),
    // ModelFree or Logit
    // loss, m
    LossOne(Vec<f64>, usize),
    // Interaction
    // vec for m and vec<vec<>> for m_inter
    // (loss, loss_inter)
    // m = loss.len()
    // m_inter = loss_inter.len()
    Interaction(Vec<f64>, Vec<f64>),
    //Interaction(Vec<f64>, usize, Vec<f64>, usize),
    //Interaction(Vec<f64>, usize, Vec<f64>, usize, usize),
    //Interaction(Vec<f64>, usize, Vec<Vec<f64>>, usize, usize),
}

impl LossStruct {
    pub fn inner_mut(&mut self) -> &mut [f64] {
        match self {
            LossStruct::ConstAda(ref mut v, _, _) => v,
            LossStruct::LossOne(ref mut v, _) => v,
            LossStruct::Interaction(ref mut v, ..) => v,
            //LossStruct::Interaction(..) => {
            //    unimplemented!("better use inner_single() or inner_interaction()")
            //}
        }
    }

    pub fn inner(&self) -> &[f64] {
        match self {
            LossStruct::ConstAda(v, _, _) => v,
            LossStruct::LossOne(v, _) => v,
            //LossStruct::ConstAda(ref v, _, _) => v,
            //LossStruct::LossOne(ref v, _) => v,
            // use inner_single() to avoid confusion;
            LossStruct::Interaction(..) => unimplemented!(),
        }
    }

    pub fn inner_single(&self) -> &[f64] {
        match self {
            LossStruct::Interaction(v, _, ..) => v,
            _ => unimplemented!("use innter_single()"),
        }
    }

    // mainly for interaction
    pub fn inner_single_mut(&mut self) -> &mut [f64] {
        match self {
            LossStruct::Interaction(ref mut v, _, ..) => v,
            _ => unimplemented!("use innter_mut()"),
        }
    }

    pub fn inner_interaction(&mut self) -> &[f64] {
        match self {
            // should be ok
            LossStruct::Interaction(_, v_inter) => v_inter,
            //LossStruct::Interaction(_, _, ref mut v_inter, _) => v_inter,
            _ => unimplemented!("use innter_mut()"),
        }
    }

    pub fn inner_interaction_mut(&mut self) -> &mut [f64] {
        match self {
            // should be ok
            LossStruct::Interaction(_, v_inter) => v_inter,
            //LossStruct::Interaction(_, _, ref mut v_inter, _) => v_inter,
            _ => unimplemented!("use innter_mut()"),
        }
    }

    /// for single snv
    pub fn access_loss(&self, mi: usize) -> f64 {
        match self {
            LossStruct::ConstAda(..) => unimplemented!(),
            LossStruct::LossOne(..) => self.inner()[mi],
            LossStruct::Interaction(..) => self.inner_single()[mi],
        }
    }

    pub fn new(boost_type: BoostType, m: usize) -> Self {
        match boost_type {
            BoostType::ConstAda => LossStruct::ConstAda(vec![f64::NAN; m * 2], m, 2),
            BoostType::FreeModelMissing
            | BoostType::Logit
            | BoostType::LogitNoMissing
            | BoostType::LogitAdd
            | BoostType::LogitMhcNoMissing
            | BoostType::LogitCommon => LossStruct::LossOne(vec![f64::NAN; m], m),
            BoostType::LogitAddInteraction => panic!("use new_interaction()"),
            //BoostType::LogitAddInteraction=> LossStruct::Interaction(vec![f64::NAN; m], m, vec![vec![f64::NAN; m2]; m1], m1, m2),
            BoostType::Ada => unimplemented!(),
        }
    }

    pub fn new_vec(boost_type: BoostType, v: Vec<f64>) -> Self {
        let m = v.len();
        match boost_type {
            BoostType::FreeModelMissing
            | BoostType::Logit
            | BoostType::LogitNoMissing
            | BoostType::LogitAdd
            | BoostType::LogitMhcNoMissing
            | BoostType::LogitCommon => LossStruct::LossOne(v, m),
            _ => unimplemented!(),
        }
    }

    pub fn new_interaction(boost_type: BoostType, m: usize, m_inter: usize) -> Self {
        match boost_type {
            BoostType::LogitAddInteraction => {
                LossStruct::Interaction(vec![f64::NAN; m], vec![f64::NAN; m_inter])
                //LossStruct::Interaction(vec![f64::NAN; m], m, vec![f64::NAN; m_inter], m_inter)
            }
            _ => panic!("use new()"),
        }
    }

    pub fn new_interaction_capacity(
        boost_type: BoostType,
        m: usize,
        m_inter: usize,
        m_inter_capacity: usize,
    ) -> Self {
        match boost_type {
            BoostType::LogitAddInteraction => {
                let mut v_inter: Vec<f64> = Vec::with_capacity(m_inter_capacity);
                v_inter.resize(m_inter, f64::NAN);

                assert_eq!(v_inter.len(), m_inter);
                assert_eq!(v_inter.capacity(), m_inter_capacity);

                LossStruct::Interaction(vec![f64::NAN; m], v_inter)
                //LossStruct::Interaction(vec![f64::NAN; m], vec![f64::NAN; m_inter])
                //LossStruct::Interaction(vec![f64::NAN; m], m, vec![f64::NAN; m_inter], m_inter)
            }
            _ => panic!("use new()"),
        }
    }

    pub fn new_interaction_vec(boost_type: BoostType, v: Vec<f64>, v_inter: Vec<f64>) -> Self {
        match boost_type {
            BoostType::LogitAddInteraction => LossStruct::Interaction(v, v_inter),
            _ => panic!("use new()"),
        }
    }

    pub fn resize_interaction(&mut self, m_inter: usize) {
        if let LossStruct::Interaction(_, loss_inter) = self {
            //if let LossStruct::Interaction(_, _, ref mut loss_inter) = self {
            // unnecessary
            //loss_inter.clear();
            loss_inter.resize(m_inter, f64::NAN);
            assert_eq!(loss_inter.len(), m_inter);
        } else {
            panic!("Sth wrong.")
        }
    }

    //pub fn new_vec(boost_type: BoostType, vec: Vec<f64>) -> Self {
    //    let m = vec.len();
    //    match boost_type {
    //        BoostType::ConstAda => LossStruct::ConstAda(vec, m / 2, 2),
    //        BoostType::FreeModelMissing
    //        | BoostType::Logit
    //        | BoostType::LogitNoMissing
    //        | BoostType::LogitAdd
    //        | BoostType::LogitMhcNoMissing => LossStruct::LossOne(vec, m),
    //        BoostType::Ada => unimplemented!(),
    //    }
    //}

    pub fn extend_interaction(&mut self, loss_2: &Self) {
        match self {
            LossStruct::Interaction(_, loss_inter) => {
                if let LossStruct::Interaction(_, loss_inter_2) = loss_2 {
                    loss_inter.extend_from_slice(loss_inter_2);
                } else {
                    panic!("loss_2 should have same type as self.");
                }
            }
            _ => panic!("Sth wrong."),
        }
    }

    pub fn extract_interaction(&mut self, use_snvs: &Vec<bool>) {
        match self {
            LossStruct::Interaction(_, ref mut loss_inter) => {
                let mut mi_ext = 0;
                for mi in 0..use_snvs.len() {
                    if use_snvs[mi] {
                        loss_inter[mi_ext] = loss_inter[mi];
                        mi_ext += 1;
                    }
                }
                // TODO: keep large vec (size=SIZE_MAX)
                loss_inter.resize(mi_ext, f64::NAN);
                log::debug!("loss extract {:?}", mi_ext);

                // TODO: do not create large vec again.
                //let loss_inter_ = loss_inter
                //    .iter()
                //    .zip(use_snvs.iter())
                //    .filter(|(_, &b)| b)
                //    .map(|(l, _)| *l)
                //    .collect::<Vec<f64>>();
                //*loss_inter = loss_inter_;
            }
            _ => panic!("Sth wrong."),
        }
    }

    //pub fn search_min(&self, skip_snv: &[usize]) -> (f64, Option<usize>, Option<usize>) {
    //pub fn search_min(&self, skip_snv: &HashSet<usize>) -> (f64, Option<usize>, Option<usize>) {
    pub fn search_min(
        &self,
        extract_snvs: Option<&HashSet<usize>>,
    ) -> (f64, Option<usize>, Option<usize>) {
        match self {
            LossStruct::ConstAda(..) => self.search_min_constada(extract_snvs),
            LossStruct::LossOne(..) => self.search_min_single(extract_snvs),
            LossStruct::Interaction(..) => self.search_min_single(extract_snvs),
            //LossStruct::Logit(..) => self.search_min_logit(skip_snv),
        }
    }

    pub fn search_min_interaction(
        &self,
        extract_interaction: &Vec<(usize, usize)>,
    ) -> (f64, Option<usize>, Option<usize>) {
        match self {
            LossStruct::ConstAda(..) | LossStruct::LossOne(..) => panic!("use search_min()"),
            LossStruct::Interaction(..) => {
                self.search_min_interaction_interaction(extract_interaction)
            } //LossStruct::Logit(..) => self.search_min_logit(skip_snv),
        }
    }

    fn get_loss_constada(&self, mi: usize, si: usize) -> f64 {
        if let LossStruct::ConstAda(loss, _, _) = self {
            let li = self.index_loss_constada(mi, si);
            loss[li]
        } else {
            panic!("Sth wrong.")
        }
    }

    /// from si,mi to li
    fn index_loss_constada(&self, mi: usize, si: usize) -> usize {
        2 * mi + si
        //(index_loss % 2, index_loss / 2)
    }

    /// return si,mi
    fn index_constada(&self, index_loss: usize) -> (usize, usize) {
        (index_loss % 2, index_loss / 2)
    }

    fn search_min_constada(
        &self,
        extract_snvs: Option<&HashSet<usize>>,
    ) -> (f64, Option<usize>, Option<usize>) {
        if let LossStruct::ConstAda(loss, _, _) = self {
            // since f64::INFINITY for not in extract_snvs in loss.rs
            let mut loss_min = f64::MAX;
            //let mut loss_min = f64::INFINITY;
            let mut var_si: Option<usize> = None;
            let mut var_mi: Option<usize> = None;

            for (li, loss) in loss.iter().enumerate() {
                let (_, var_mi_) = self.index_constada(li);
                if (*loss < loss_min)
                    && (extract_snvs.is_none() || extract_snvs.unwrap().contains(&var_mi_))
                {
                    loss_min = *loss;
                    let (var_si_, var_mi_) = self.index_constada(li);
                    var_si = Some(var_si_);
                    var_mi = Some(var_mi_);
                    /*
                    log::debug!(
                        "Temporary chosen index, si, mi, loss: {}, {}, {}, {:.6}",
                        i,
                        var_si.unwrap(),
                        var_mi.unwrap(),
                        *loss
                    );
                    */
                }
            }
            return (loss_min, var_mi, var_si);
        } else {
            panic!("Sth wrong.")
        }
    }

    pub fn flip_sign(&self) -> Self {
        match self {
            LossStruct::ConstAda(loss, m, s) => {
                let mut loss = loss.clone();
                for i in 0..loss.len() {
                    loss[i] = -loss[i];
                }
                LossStruct::ConstAda(loss, *m, *s)
            }
            LossStruct::LossOne(loss, m) => {
                let mut loss = loss.clone();
                for i in 0..loss.len() {
                    loss[i] = -loss[i];
                }
                LossStruct::LossOne(loss, *m)
            }
            LossStruct::Interaction(..) => unimplemented!(),
        }
    }

    fn search_min_single(
        &self,
        extract_snvs: Option<&HashSet<usize>>,
        //skip_snv: &HashSet<usize>,
    ) -> (f64, Option<usize>, Option<usize>) {
        //if let LossStruct::LossOne(loss, _) = self {
        let (loss_min, var_mi) = match self {
            LossStruct::LossOne(loss, _) => {
                // for Logit
                Self::search_min_single_vec(loss, extract_snvs)
                //let mut loss_min = f64::INFINITY;
                //let mut var_mi: Option<usize> = None;

                ////log::debug!("loss m {}", loss.len());
                //for (mi, loss) in loss.iter().enumerate() {
                //    //if mi % 100000 == 0 {
                //    //    log::debug!("min loss mi {} {}", mi, *loss);
                //    //}

                //    // if (*loss < loss_min) && !(skip_snv.contains(&mi)) {
                //    if (*loss < loss_min)
                //        && (extract_snvs.is_none() || extract_snvs.unwrap().contains(&mi))
                //    {
                //        //if *loss < loss_min {
                //        loss_min = *loss;
                //        let var_mi_ = mi;
                //        var_mi = Some(var_mi_);
                //        /*
                //        log::debug!(
                //            "Temporary chosen index, mi, loss: {}, {:.6}",
                //            //var_si.unwrap(),
                //            var_mi.unwrap(),
                //            *loss
                //        );
                //         */
                //    }
                //}
            }
            LossStruct::Interaction(loss, ..) => Self::search_min_single_vec(loss, extract_snvs),
            _ => unreachable!(),
        };

        if var_mi.is_none() {
            panic!("No SNVs were selected");
        }

        log::debug!(
            "Selected index, mi, loss: {}, {:.6}",
            var_mi.unwrap(),
            loss_min
        );
        return (loss_min, var_mi, None);
    }

    // TODO: test
    // check loss[mi] is nan if mi in extract_snvs -> unnecessary
    fn search_min_single_vec(
        loss: &Vec<f64>,
        extract_snvs: Option<&HashSet<usize>>,
    ) -> (f64, Option<usize>) {
        // for Logit
        let mut loss_min = f64::MAX;
        //let mut loss_min = f64::INFINITY;
        let mut var_mi: Option<usize> = None;

        //log::debug!("loss m {}", loss.len());
        for (mi, loss) in loss.iter().enumerate() {
            //if mi % 100000 == 0 {
            //    log::debug!("min loss mi {} {}", mi, *loss);
            //}

            // if (*loss < loss_min) && !(skip_snv.contains(&mi)) {
            if (*loss < loss_min) && (extract_snvs.is_none() || extract_snvs.unwrap().contains(&mi))
            {
                //if *loss < loss_min {
                loss_min = *loss;
                let var_mi_ = mi;
                var_mi = Some(var_mi_);
                /*
                log::debug!(
                    "Temporary chosen index, mi, loss: {}, {:.6}",
                    //var_si.unwrap(),
                    var_mi.unwrap(),
                    *loss
                );
                 */
            }
        }

        (loss_min, var_mi)
    }

    // check loss[mi] is nan if mi in extract_snvs
    fn search_min_interaction_interaction(
        &self,
        extract_interaction: &Vec<(usize, usize)>,
    ) -> (f64, Option<usize>, Option<usize>) {
        if let LossStruct::Interaction(_, loss_interaction) = self {
            if loss_interaction.len() == 0 {
                return (f64::MAX, None, None);
            }

            if loss_interaction.len() != extract_interaction.len() {
                panic!("loss_interaction.len()!=extract_interaction.len()");
            }

            // for Logit
            let mut loss_min = f64::MAX;
            //let mut loss_min = f64::INFINITY;
            let mut var_mi_pair: Option<(usize, usize)> = None;

            //log::debug!("loss m {}", loss.len());

            for (&snv_pair, loss) in extract_interaction.iter().zip(loss_interaction.iter()) {
                //if mi % 100000 == 0 {
                //    log::debug!("min loss mi {} {}", mi, *loss);
                //}

                // if (*loss < loss_min) && !(skip_snv.contains(&mi)) {
                if *loss < loss_min {
                    loss_min = *loss;
                    var_mi_pair = Some(snv_pair);
                    /*
                    log::debug!(
                        "Temporary chosen index, mi, loss: {}, {:.6}",
                        //var_si.unwrap(),
                        var_mi.unwrap(),
                        *loss
                    );
                     */
                }
            }

            if var_mi_pair.is_none() {
                panic!("No SNVs were selected");
            }

            log::debug!(
                "Selected index, mi, loss: {:?}, {:.6}",
                var_mi_pair.unwrap(),
                loss_min
            );

            return (
                loss_min,
                Some(var_mi_pair.unwrap().0),
                Some(var_mi_pair.unwrap().1),
            );
        } else {
            unreachable!();
        }
    }

    pub fn search_topprop(&self, prop: f64) -> Vec<bool> {
        match self {
            LossStruct::ConstAda(..) => unimplemented!(),
            LossStruct::LossOne(..) => self.search_topprop_modelfree(prop),
            LossStruct::Interaction(..) => unimplemented!(),
        }
    }

    fn search_topprop_modelfree(&self, prop: f64) -> Vec<bool> {
        if let LossStruct::LossOne(loss, _) = self {
            let m = loss.len();
            // takes at least 1 SNVs
            let topprop_n = ((m as f64 * prop) as usize).max(1);

            let (use_snvs, _) = self.search_topprop_modelfree_n(topprop_n, None);
            return use_snvs;
        } else {
            panic!("Sth wrong.")
        }
    }

    // for batch
    // return next loss
    pub fn search_topprop_n(
        &self,
        m: usize,
        extract_snvs: Option<&HashSet<usize>>,
    ) -> (Vec<bool>, f64) {
        match self {
            LossStruct::ConstAda(..) => unimplemented!(),
            LossStruct::LossOne(..) => self.search_topprop_modelfree_n(m, extract_snvs),
            LossStruct::Interaction(..) => self.search_topprop_modelfree_n(m, extract_snvs),
        }
    }

    // for batch interactions
    // return next loss
    pub fn search_topprop_interaction_n(
        &self,
        m: usize,
        //extract_interaction: &Vec<(usize, usize)>,
    ) -> (Vec<bool>, f64) {
        match self {
            LossStruct::ConstAda(..) | LossStruct::LossOne(..) => panic!("use search_topprop_n()"),
            LossStruct::Interaction(..) => {
                self.search_topprop_n_interaction(m)
                //self.search_topprop_n_interaction(m, extract_interaction)
            }
        }
    }

    // TODO: test
    /// also return next loss
    /// if extract snvs < m_top, return # extract snvs
    fn search_topprop_modelfree_n(
        &self,
        m_top: usize,
        extract_snvs: Option<&HashSet<usize>>,
    ) -> (Vec<bool>, f64) {
        //if let LossStruct::LossOne(loss, _)  = self {
        match self {
            LossStruct::LossOne(loss, _) | LossStruct::Interaction(loss, ..) => {
                //let mut loss_topprop = 10.0;

                let m = loss.len();

                //log::debug!("loss {:?}", loss);
                if loss.iter().any(|x| x.is_nan()) {
                    log::debug!("loss any nan {:?}", loss.iter().any(|x| x.is_nan()));
                    log::debug!(
                        "loss nan {:?}",
                        loss.iter()
                            .filter(|x| x.is_nan())
                            .map(|x| *x)
                            .collect::<Vec<f64>>()
                    );
                }

                let mut loss_sort = loss.clone();

                loss_sort.iter_mut().for_each(|x| {
                    if x.is_nan() {
                        *x = f64::MAX
                    }
                });

                // make loss of skip_snv f64::MAX
                if let Some(extract_snvs) = extract_snvs {
                    loss_sort.iter_mut().enumerate().for_each(|(i, x)| {
                        if !extract_snvs.contains(&i) {
                            *x = f64::MAX
                        }
                    });
                    //skip_snv.iter().for_each(|i| loss_sort[*i] = f64::MAX);
                }
                vec::sort_float(&mut loss_sort);

                // takes at least 1 SNVs
                //let topprop_n = ((m as f64 * prop) as usize).max(1);

                let top_n = m_top.min(m);
                if top_n >= loss_sort.len() {
                    return (vec![true; m], f64::MAX);
                }
                // FIXED: BUG
                // -2 for loss_topprop_top_outside
                //let topprop_n = m_top.min(m) - 2;

                // the smaller, the better
                let loss_topprop = loss_sort[top_n - 1];
                let loss_topprop_top_outside = loss_sort[top_n];

                let mut use_snvs = vec![false; m];

                for (mi, loss) in loss.iter().enumerate() {
                    if (*loss <= loss_topprop)
                        && (extract_snvs.is_none() || extract_snvs.unwrap().contains(&mi))
                    {
                        //if (*loss <= loss_topprop) && !(skip_snv.contains(&mi)) {
                        use_snvs[mi] = true;
                    }
                }

                //if let Some(skip_snv) = skip_snv {
                //    for (mi, loss) in loss.iter().enumerate() {
                //        //if *loss <= loss_topprop {
                //        if (*loss <= loss_topprop) && !(skip_snv.contains(&mi)) {
                //            use_snvs[mi] = true;
                //        }
                //    }
                //} else {
                //    for (mi, loss) in loss.iter().enumerate() {
                //        if *loss <= loss_topprop {
                //            use_snvs[mi] = true;
                //        }
                //    }
                //}

                log::debug!(
                    "#SNVs in use_snvs by loss {} in {}",
                    use_snvs.iter().filter(|b| **b).count(),
                    m
                );

                return (use_snvs, loss_topprop_top_outside);
                //return (loss_min, var_mi, var_si);
            }
            _ => panic!("Sth wrong."),
        }
    }

    // TODO: commo with search_topprop_modelfree_n()
    /// also return next loss
    fn search_topprop_n_interaction(
        &self,
        m_top: usize,
        //extract_interaction: &Vec<(usize, usize)>,
    ) -> (Vec<bool>, f64) {
        //if let LossStruct::LossOne(loss, _)  = self {
        match self {
            // why error?
            //LossStruct::Interaction(_, _, loss, &m) => {
            LossStruct::Interaction(_, loss) => {
                let m = loss.len();
                //log::debug!("loss len {:?}",m);

                //log::debug!("loss {:?}", loss);
                if loss.iter().any(|x| x.is_nan()) {
                    log::debug!("loss any nan {:?}", loss.iter().any(|x| x.is_nan()));
                    log::debug!(
                        "loss nan {} / {}",
                        loss.iter().filter(|x| x.is_nan()).map(|x| *x).count(),
                        loss.len() //.collect::<Vec<f64>>().len()
                    );
                }

                let mut loss_sort = loss.clone();

                loss_sort.iter_mut().for_each(|x| {
                    if x.is_nan() {
                        *x = f64::MAX
                    }
                });

                // no skip for interaction
                //// make loss of skip_snv f64::MAX
                //loss_sort.iter_mut().enumerate().for_each(|(i, x)| {
                //    if !extract_snvs.contains(&i) {
                //        *x = f64::MAX
                //    }
                //});
                //skip_snv.iter().for_each(|i| loss_sort[*i] = f64::MAX);
                vec::sort_float(&mut loss_sort);

                // takes at least 1 SNVs
                //let topprop_n = ((m as f64 * prop) as usize).max(1);

                let top_n = m_top.min(m);
                if top_n >= loss_sort.len() {
                    return (vec![true; m], f64::MAX);
                }

                //let topprop_n = top_n - 1;
                // FIXED: BUG
                // -2 for loss_topprop_top_outside
                //let topprop_n = m_top.min(m) - 2;

                // the smaller, the better
                let loss_topprop = loss_sort[top_n - 1];
                let loss_topprop_top_outside = loss_sort[top_n];
                //let loss_topprop = loss_sort[topprop_n];
                //let loss_topprop_top_outside = loss_sort[topprop_n + 1];

                let mut use_snvs = vec![false; m];

                // # of true SNVs might be larger than m when multiple SNVs have the same loss
                for (mi, loss) in loss.iter().enumerate() {
                    if *loss <= loss_topprop {
                        use_snvs[mi] = true;
                    }
                }

                //if let Some(skip_snv) = skip_snv {
                //    for (mi, loss) in loss.iter().enumerate() {
                //        //if *loss <= loss_topprop {
                //        if (*loss <= loss_topprop) && !(skip_snv.contains(&mi)) {
                //            use_snvs[mi] = true;
                //        }
                //    }
                //} else {
                //    for (mi, loss) in loss.iter().enumerate() {
                //        if *loss <= loss_topprop {
                //            use_snvs[mi] = true;
                //        }
                //    }
                //}

                log::debug!(
                    "#SNVs in use_snvs_interaction by loss {} in {}",
                    use_snvs.iter().filter(|b| **b).count(),
                    m
                );

                return (use_snvs, loss_topprop_top_outside);
                //return (loss_min, var_mi, var_si);
            }
            _ => panic!("Sth wrong."),
        }
    }

    pub fn decrease_by(&self, loss_2: &Self) -> Self {
        match self {
            LossStruct::ConstAda(..) => {
                unimplemented!()
                //let mut loss_2 = loss_2.inner();
                //let mut loss = loss.inner_mut();
                //for i in 0..loss.len() {
                //    loss[i] -= loss_2[i];
                //}
                //LossStruct::ConstAda(loss.to_vec(), *m, *s)
            }
            LossStruct::LossOne(loss, m) => {
                if let LossStruct::LossOne(loss_2, _) = loss_2 {
                    let mut loss_by = loss.clone();
                    //let mut loss_2 = loss_2.clone();
                    //let mut loss = loss.clone();
                    for i in 0..loss.len() {
                        loss_by[i] = loss_2[i] - loss[i];
                        //loss[i] -= loss_2[i];
                    }
                    LossStruct::LossOne(loss_by, *m)
                } else {
                    panic!("loss_2 should have same type as self.");
                }
            }
            LossStruct::Interaction(..) => unimplemented!(),
        }
    }

    // FIXME constada
    // modelfree only now
    fn string_write(&self, snvs: &Snvs) -> String {
        if let LossStruct::LossOne(loss, _) = self {
            log::debug!("loss[0], [1]: {}, {}", loss[0], loss[1]);
            let str_header = "m\trs\tchrom\tpos\tA1\tA2\tloss".to_owned();
            let strings = loss
                .iter()
                .zip(snvs.snv_ids().iter())
                .enumerate()
                .map(|(mi, (l, snv_id))| {
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        mi.to_string(),
                        snv_id.id(),
                        snv_id.chrom().to_string(),
                        snv_id.pos().to_string(),
                        snv_id.a1(),
                        snv_id.a2(),
                        l.to_string()
                    )
                    //mi.to_string()
                    //    + "\t"
                    //    + snv_id.id()
                    //    + "\t"
                    //    + &snv_id.chrom().to_string()
                    //    + "\t"
                    //    + &snv_id.pos().to_string()
                    //    + "\t"
                    //    + snv_id.a1()
                    //    + "\t"
                    //    + snv_id.a2()
                    //    + "\t"
                    //    + &l.to_string()
                })
                .collect::<Vec<String>>();
            str_header + "\n" + &strings.join("\n")
        } else {
            unimplemented!();
        }
    }

    // only write :.3
    fn string_write_initial(
        &self,
        snvs_interaction: &Vec<(usize, usize)>,
        alphas: &Vec<f64>,
    ) -> String {
        if let LossStruct::Interaction(_loss_single, loss_inter) = self {
            //log::debug!(
            //    "loss_single[0], [1]: {}, {}",
            //    loss_single[0],
            //    loss_single[1]
            //);
            log::debug!("loss_inter[0], [1]: {}, {}", loss_inter[0], loss_inter[1]);
            let str_header = "m1\tm2\talpha\tloss".to_owned();
            // use string_write_initial_single()
            //// for single
            //let strings = loss_single
            //    .iter()
            //    .zip(alphas.iter())
            //    .enumerate()
            //    .map(|(mi, (l, alpha))| format!("{}\t{}\t{:.5e}\t{:.3}", mi, "NaN", alpha, l))
            //    .collect::<Vec<String>>();
            // for interaction
            let strings_inter = loss_inter
                .iter()
                .zip(alphas.iter())
                .zip(snvs_interaction.iter())
                .map(|((l, alpha), (mi_1, mi_2))| {
                    format!("{}\t{}\t{:.4e}\t{:.4e}", mi_1, mi_2, alpha, l)
                })
                .collect::<Vec<String>>();
            str_header + "\n" + &strings_inter.join("\n")
            //str_header + "\n" + &strings.join("\n") + "\n" + &strings_inter.join("\n")
        } else {
            unimplemented!();
        }
    }

    fn string_write_initial_single(&self, alphas: &Vec<f64>) -> String {
        if let LossStruct::Interaction(loss_single, _) = self {
            log::debug!(
                "loss_single[0], [1]: {}, {}",
                loss_single[0],
                loss_single[1]
            );
            let str_header = "m\talpha\tloss".to_owned();
            let strings = loss_single
                .iter()
                .zip(alphas.iter())
                .enumerate()
                .map(|(mi, (l, alpha))| format!("{}\t{:.4e}\t{:.4e}", mi, alpha, l))
                .collect::<Vec<String>>();
            str_header + "\n" + &strings.join("\n")
        } else {
            unimplemented!();
        }
    }

    //pub fn write_writer<W: std::io::Write>(&self, writer: &mut BufWriter<W>, iteration: usize) {
    pub fn write_writer<W: std::io::Write>(&self, writer: &mut BufWriter<W>, snvs: &Snvs) {
        /*
        let fwgt = io::get_fname_wgt(fout);
        let file = OpenOptions::new().append(true).open(fwgt).unwrap();
        //let mut writer = BufWriter::new(File::create(fout).unwrap());
        let mut writer = BufWriter::new(file);
         */
        let str = self.string_write(snvs);
        //log::debug!("str wgt {}", &str);
        //let str = "a\tb\n".to_owned();
        writer.write_all(str.as_bytes()).unwrap();
        //writer.write(str.as_bytes()).unwrap();
        // capture error
        writer.flush().unwrap();
    }

    // initial has large file size so only write snv index and loss
    // save both single and interaction
    pub fn write_initial_writer<W: std::io::Write>(
        &self,
        writer: &mut BufWriter<W>,
        snvs_interaction: &Vec<(usize, usize)>,
        alphas: &Vec<f64>,
    ) {
        let str = self.string_write_initial(snvs_interaction, alphas);
        writer.write_all(str.as_bytes()).unwrap();
        // This stops at 2GB
        //writer.write(str.as_bytes()).unwrap();
        // capture error
        writer.flush().unwrap();
    }

    pub fn write_initial_single_writer<W: std::io::Write>(
        &self,
        writer: &mut BufWriter<W>,
        alphas: &Vec<f64>,
        //snvs_interaction: &Vec<(usize, usize)>,
    ) {
        let str = self.string_write_initial_single(alphas);
        writer.write_all(str.as_bytes()).unwrap();
        writer.flush().unwrap();
    }
}

pub fn search_min_loss_gt(
    loss: &mut LossStruct,
    iteration: usize,
    dataset: &Dataset,
    //genot: &Genot,
    //sample_weight: &SampleWeight,
    //phe: &Phe,
    //snvs: &Snvs,
    boost_param: &BoostParam,
    //boost_param: BoostParam,
    extract_snvs: Option<&HashSet<usize>>,
    //skip_snv: &HashSet<usize>,
    //use_adjloss: bool,
) -> WgtBoost {
    let genot = dataset.genot();
    let phe = dataset.samples().phe_unwrap();
    let snvs = dataset.snvs();

    // calc loss: moved up
    //calculate_loss_gt(loss, dataset, sample_weight, boost_param, skip_snv);
    //calculate_loss_gt(loss, genot, sample_weight, phe, boost_param, skip_snv);

    // search min loss
    let (loss_min, var_mi, var_si) = loss.search_min(extract_snvs);
    let var_mi: usize = var_mi.expect("No SNVs were loss<1.0");

    log::debug!(
        "contingency table count {:?}",
        genot
            .to_genot_snv(var_mi)
            .stat_contingency_table_nosimd(phe)
    );

    // create wgt
    match boost_param.boost_type() {
        BoostType::ConstAda => {
            let var_si: usize = var_si.expect("No SNVs were loss<1.0");
            let wgt_min = Wgt::new_snv_id(
                snvs.snv_ids()[var_mi].clone(),
                genetics::THRESHOLD_SNV[var_si],
                (var_si, var_mi),
                snvs.mafs().unwrap()[var_mi],
            );

            let wgt_boost_min = WgtBoost::new_wgt_loss(wgt_min, iteration, Some(loss_min));
            wgt_boost_min
        }

        BoostType::FreeModelMissing
        | BoostType::Logit
        | BoostType::LogitNoMissing
        | BoostType::LogitAdd
        | BoostType::LogitMhcNoMissing
        | BoostType::LogitCommon => {
            let wgt_min = Wgt::new_snv_id_freemodel(
                snvs.snv_ids()[var_mi].clone(),
                var_mi,
                snvs.mafs().unwrap()[var_mi],
            );
            let wgt_boost_min = WgtBoost::new_wgt_loss(wgt_min, iteration, Some(loss_min));
            wgt_boost_min
        }
        BoostType::LogitAddInteraction => unimplemented!(),

        BoostType::Ada => {
            unimplemented!()
        }
    }
}

pub fn search_min_loss_gt_interaction(
    loss: &mut LossStruct,
    iteration: usize,
    dataset: &Dataset,
    _boost_param: &BoostParam,
    extract_snvs: Option<&HashSet<usize>>,
    extract_interaction: &Vec<(usize, usize)>,
) -> WgtBoost {
    let genot = dataset.genot();
    let phe = dataset.samples().phe_unwrap();
    let snvs = dataset.snvs();

    // calc loss: moved up
    //calculate_loss_gt(loss, dataset, sample_weight, boost_param, skip_snv);
    //calculate_loss_gt(loss, genot, sample_weight, phe, boost_param, skip_snv);

    // search min loss for single snv
    let (loss_min_single, var_mi_single, _) = loss.search_min(extract_snvs);

    // search min loss for interaction
    let (loss_min_interaction, var_mi_1, var_mi_2) = if extract_interaction.len() == 0 {
        (f64::MAX, None, None)
        //(f64::INFINITY, None, None)
    } else {
        //let (loss_min_interaction, var_mi_1, var_mi_2) =
        loss.search_min_interaction(extract_interaction)
    };

    // directly compare loss_min and loss_min_interaction
    if loss_min_single < loss_min_interaction {
        // adopt single
        let var_mi_single: usize = var_mi_single.expect("No SNVs were selected.");

        log::debug!(
            "contingency table count {:?}",
            genot
                .to_genot_snv(var_mi_single)
                .stat_contingency_table_nosimd(phe)
        );

        let wgt_min = Wgt::new_snv_id_freemodel(
            snvs.snv_ids()[var_mi_single].clone(),
            var_mi_single,
            snvs.mafs().unwrap()[var_mi_single],
        );
        let wgt_boost_min = WgtBoost::new_wgt_loss(wgt_min, iteration, Some(loss_min_single));
        wgt_boost_min
    } else {
        // adopt intearction
        let var_mi_1: usize = var_mi_1.expect("No SNVs were selected.");
        let var_mi_2: usize = var_mi_2.expect("No SNVs were selected.");

        //log::debug!(
        //    "contingency table count {:?}",
        //    genot
        //        .to_genot_snv(var_mi_single)
        //        .stat_contingency_table_nosimd(phe)
        //);

        // adopt single
        let wgt_min = Wgt::new_snv_id_interaction(
            snvs.snv_ids()[var_mi_1].clone(),
            snvs.snv_ids()[var_mi_2].clone(),
            var_mi_1,
            var_mi_2,
            snvs.mafs().unwrap()[var_mi_1],
            snvs.mafs().unwrap()[var_mi_2],
        );
        let wgt_boost_min = WgtBoost::new_wgt_loss(wgt_min, iteration, Some(loss_min_single));
        wgt_boost_min
    }

    // create wgt
    //match boost_param.boost_type() {
    //    BoostType::FreeModelMissing
    //    | BoostType::Logit
    //    | BoostType::LogitNoMissing
    //    | BoostType::LogitAdd
    //    | BoostType::LogitMhcNoMissing
    //    | BoostType::LogitCommon => {
    //        let wgt_min = Wgt::new_snv_index_freemodel(snvs.snv_indexs()[var_mi].clone(), var_mi);
    //        let wgt_boost_min = WgtBoost::new_wgt_loss(wgt_min, iteration, Some(loss_min));
    //        wgt_boost_min
    //    }
    //    BoostType::LogitAddInteraction => unimplemented!(),

    //    BoostType::Ada => {
    //        unimplemented!()
    //    }
    //}
}

pub fn calculate_loss_gt(
    loss: &mut LossStruct,
    dataset: &Dataset,
    sample_weight: &mut SampleWeight,
    boost_param: &BoostParam,
    extract_snvs: Option<&HashSet<usize>>,
    //skip_snv: &HashSet<usize>,
    //use_adjloss: bool,
    alphas_save: Option<&mut Vec<f64>>,
) {
    // should be done in score but rerun for safety
    sample_weight.clear_pad();

    let func = match boost_param.boost_type() {
        BoostType::Logit | BoostType::LogitNoMissing => calc::calc_loss_logit,
        BoostType::LogitAdd => calc::calc_loss_logit_add,
        BoostType::LogitMhcNoMissing => calc::calc_loss_logit_mhcnomissing,
        BoostType::LogitCommon => calc::calc_loss_logit_common,
        BoostType::ConstAda => calc::calc_loss_constada,
        BoostType::FreeModelMissing => calc::calc_loss_freemodelmissing,
        BoostType::Ada => {
            unimplemented!()
        }
        // for single SNVs
        BoostType::LogitAddInteraction => calc::calc_loss_logit_add,
    };

    func(
        loss,
        dataset,
        sample_weight,
        boost_param,
        extract_snvs,
        alphas_save,
    );

    //#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    //{
    //    // TOFIX: This would raise error when no avx2 but on x86??,
    //    // since it is not compiled with avx2
    //    if is_x86_feature_detected!("avx2") {
    //        log::debug!("Use SIMD.");
    //        unsafe {
    //            let func = match boost_param.boost_type() {
    //                BoostType::Logit | BoostType::LogitNoMissing => {
    //                    calc::calculate_loss_gt_logit_simd
    //                }
    //                BoostType::LogitAdd => calc::calculate_loss_gt_logit_add_simd,
    //                BoostType::LogitMhcNoMissing => calc::calculate_loss_gt_logit_mhcnomissing_simd,
    //                BoostType::LogitCommon => calc::calculate_loss_gt_logit_common_simd,
    //                BoostType::ConstAda => calc::calculate_loss_gt_constada_simd,
    //                BoostType::FreeModelMissing => calc::calculate_loss_gt_freemodelmissing_simd,
    //                BoostType::Ada => {
    //                    unimplemented!()
    //                }
    //                // for single SNVs
    //                BoostType::LogitAddInteraction => calc::calculate_loss_gt_logit_add_simd,
    //            };

    //            func(loss, dataset, sample_weight, boost_param, extract_snvs);
    //        }
    //        return;
    //    }
    //}

    //log::debug!("Do not use SIMD.");

    //match boost_param.boost_type() {
    //    BoostType::ConstAda => {
    //        calc::calculate_loss_gt_constada_nosimd(
    //            loss,
    //            dataset,
    //            //genot,
    //            sample_weight,
    //            //sample_weight.ps().unwrap(),
    //            //phe,
    //            boost_param,
    //            extract_snvs,
    //        );
    //    }
    //    BoostType::FreeModelMissing => {
    //        calc::calculate_loss_gt_freemodelmissing_nosimd(
    //            loss,
    //            dataset,
    //            //genot,
    //            sample_weight,
    //            //sample_weight.ps().unwrap(),
    //            //phe,
    //            boost_param,
    //            extract_snvs,
    //        );
    //    }
    //    // FIXME: BoostType::LogitNoMissing
    //    _ => unimplemented!(),
    //}
}

pub fn calculate_loss_gt_interaction(
    loss: &mut LossStruct,
    dataset: &Dataset,
    sample_weight: &SampleWeight,
    boost_param: &BoostParam,
    extract_snvs: Option<&HashSet<usize>>,
    extract_interaction: &Vec<(usize, usize)>,
    //extract_interaction: Option<&HashSet<(usize, usize)>>,
) {
    match boost_param.boost_type() {
        BoostType::LogitAddInteraction => {
            calc::calc_loss_logitaddinteraction(
                loss,
                dataset,
                sample_weight,
                boost_param,
                extract_snvs,
                extract_interaction,
            );
        }
        _ => {
            panic!("");
        }
    }

    //#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    //{
    //    if is_x86_feature_detected!("avx2") {
    //        log::debug!("Use SIMD.");
    //        unsafe {
    //            match boost_param.boost_type() {
    //                BoostType::LogitAddInteraction => {
    //                    calc::calculate_loss_gt_logitaddinteraction_simd(
    //                        loss,
    //                        dataset,
    //                        sample_weight,
    //                        boost_param,
    //                        extract_snvs,
    //                        extract_interaction,
    //                    );
    //                }
    //                _ => {
    //                    panic!("");
    //                }
    //            };
    //        }
    //        return;
    //    }
    //}
}

/* pub fn search_min_loss_gt_pruned(
    loss_gt: &[f64],
    iteration: usize,
    predictions: &GenotBi<Vec<u8>>,
    //predictions: &[B8],
    ps: &[f64],
    ys: &[B8],
    //m: usize,
    n: usize,
    snvs: &[Snv],
    is_pruned: &[bool],
) -> WgtBoost {
    //calculate::calculate_loss_gt(loss_gt, predictions, ps, ys, n);

    //log::debug!("loss_gt {:?}", loss_gt);

    let mut loss_min = 1.0;
    let mut var_si: Option<usize> = None;
    let mut var_mi: Option<usize> = None;
    /*
    let mut var_mi = usize::MAX; //NaN
    let mut var_si = usize::MAX;
    */

    for (i, loss) in loss_gt.iter().enumerate() {
        let var_mi_ = index_from_loss(i);
        if is_pruned[var_mi_] {
            continue;
        }
        if *loss < loss_min {
            loss_min = *loss;
            let (var_si_, var_mi_) = predict_index_from_loss(i);
            var_si = Some(var_si_);
            var_mi = Some(var_mi_);
            log::debug!(
                "Temporary chosen index, si, mi, loss: {}, {}, {}, {:.6}",
                i,
                var_si.unwrap(),
                var_mi.unwrap(),
                *loss
            );
        }
    }

    log::debug!(
        "Chosen si, mi, loss: {}, {}, {:.6}",
        var_si.unwrap(),
        var_mi.unwrap(),
        loss_min
    );

    let var_si: usize = var_si.expect("No SNVs were loss<1.0");
    let var_mi: usize = var_mi.expect("No SNVs were loss<1.0");

    //assert_ne!(var_mi, usize::MAX, "No SNVs were loss<1.0");

    //pub fn calculate_abcd_sum(mis_sm: &[B8], ps: &[f64], ys: &[B8], n: usize) -> (f64, f64, f64, f64) {
    let abcd_sum = calculate::calculate_abcd_sum(
        predictions.predictions_snv(var_si, var_mi),
        //predict::predictions_snv(predictions, var_si, var_mi, n),
        //mistake::get_predictions_snv(predictions, var_si, var_mi, m, n),
        ps,
        ys,
        n,
    );

    let snv_wgt = SnvWgt::construct(snvs[var_mi].clone(), (var_si, var_mi));
    //let kind = WgtKind::Snv(snv_wgt);
    let wgt_min = Wgt::construct_wgt(
        WgtKind::Snv(snv_wgt),
        Model::Binary(BinaryModel::construct_threshold(
            super::THRESHOLD_SNV[var_si],
        )),
        //iteration,
        //loss_min,
        //abcd_sum,
    );

    let wgt_boost_min = WgtBoost::construct_wgt(wgt_min, iteration, Some(loss_min), Some(abcd_sum));

    /*
    let wgt_min = Wgt::construct_wgt(
        //WgtKind::Snv(snvs[var_mi].clone()),
        WgtKind::Snv(SnvWgt::construct_snv_wgt(
            snvs[var_mi].clone(),
            (var_si, var_mi),
        )),
        iter,
        //var_si,
        //var_mi,
        loss_min,
        Some(super::THRESHOLD_SNV[var_si]),
        abcd_sum,
    );
     */
    wgt_boost_min
} */

pub fn search_min_loss_gt_second(
    wgt_first: &WgtBoost,
    loss: &mut LossStruct,
    iteration: usize,
    snvs: &Snvs,
    boost_param: &BoostParam,
) -> WgtBoost {
    let (var_si_first, var_mi) = wgt_first.wgt().kind().index_predict_snv();
    let var_si = if var_si_first == 0 {
        1
    } else if var_si_first == 1 {
        0
    } else {
        panic!("Impossible var_si_first: {}", var_si_first)
    };

    let loss_second = loss.get_loss_constada(var_mi, var_si);

    // create wgt
    match boost_param.boost_type() {
        BoostType::ConstAda => {
            let wgt_min = Wgt::new_snv_id(
                snvs.snv_ids()[var_mi].clone(),
                genetics::THRESHOLD_SNV[var_si],
                (var_si, var_mi),
                snvs.mafs().unwrap()[var_mi],
            );

            let wgt_boost_min = WgtBoost::new_wgt_loss(wgt_min, iteration, Some(loss_second));
            wgt_boost_min
        }
        _ => {
            panic!()
        }
    }
}

/*
/// set alpha too
pub fn search_min_loss_gt_pruned_loss_ss(
    sum_stats: &[SumStat],
    iteration: usize,
    is_pruned: &[bool],
) -> WgtBoost {
    let mut loss_min = 1.0;
    //let mut var_si: Option<usize> = None;
    let mut var_mi: Option<usize> = None;

    for (mi, snv) in sum_stats.iter().enumerate() {
        if is_pruned[mi] {
            continue;
        }
        // TODO; None
        let loss = snv.loss().unwrap();
        if loss < loss_min {
            loss_min = loss;
            var_mi = Some(mi);
            log::debug!(
                "Temporary chosen mi, loss: {}, {:.6}",
                var_mi.unwrap(),
                loss
            );
        }
    }

    log::debug!("Chosen mi, loss: {}, {:.6}", var_mi.unwrap(), loss_min);

    let var_mi: usize = var_mi.expect("No SNVs were loss<1.0");

    // TODO: for None
    let alpha = sum_stats[var_mi].alpha().unwrap();

    let snv_wgt = SnvWgt::construct_ss(sum_stats[var_mi].snv_index().clone(), var_mi);
    //let kind = WgtKind::Snv(snv_wgt);
    let wgt_min = Wgt::construct_wgt(
        WgtKind::Snv(snv_wgt),
        Model::Linear(LinearModel::new(Coef::Linear(alpha))),
    );

    let wgt_loss_min = WgtBoost::construct_wgt(wgt_min, iteration, Some(loss_min), None);

    wgt_loss_min
}
 */

/*
/// set alpha too
pub fn search_min_loss_gt_pruned_loss_ss_second(
    wgt_first: &WgtBoost,
    sum_stats: &[SumStat],
    iteration: usize,
) -> WgtBoost {
    let var_mi = wgt_first.wgt().kind().index_snv().unwrap();

    let loss_min = sum_stats[var_mi].loss().unwrap();

    log::debug!("Chosen mi, loss: {}, {:.6}", var_mi, loss_min);

    // TODO: for None
    let alpha = sum_stats[var_mi].alpha().unwrap();

    let snv_wgt = SnvWgt::construct_ss(sum_stats[var_mi].snv_index().clone(), var_mi);
    //let kind = WgtKind::Snv(snv_wgt);
    let wgt_min = Wgt::construct_wgt(
        WgtKind::Snv(snv_wgt),
        Model::Linear(LinearModel::new(Coef::Linear(alpha))),
        //iteration,
        //loss_min,
        //(f64::NAN, f64::NAN, f64::NAN, f64::NAN),
    );

    let wgt_boost_min = WgtBoost::construct_wgt(wgt_min, iteration, Some(loss_min), None);

    wgt_boost_min
}
*/

// // not using now
// pub fn create_loss_const(
//     iter: usize,
//     ps: &[f64],
//     pred: &mut [u8],
//     //pred: &mut [B8],
//     phe: &Phe,
//     //_n: usize,
//     //_boosting_param: BoostParam,
// ) -> WgtBoost {
//     /*
//     let mut const_var = Var::construct_var(CovKind::Const, "const".to_owned());
//     const_var.set_vals_const(n);
//     let const_wgt = CovWgt::construct(const_var);
//      */
//     compute_pred::compute_pred_const(pred);

//     //let abcd_sum = calculate::calculate_abcd_sum(
//     let (table2_sum, is_eps) = table::calculate_table2(pred, ps, phe);

//     let loss = calc::calculate_loss_ab(table2_sum);

//     // since all predict=1
//     let wgt = Wgt::construct_const_threshold(0.5);

//     /*
//     let wgt = Wgt::construct_wgt(
//         WgtKind::Cov(const_wgt),
//         Model::Binary(BinaryModel::construct_threshold(0.5)),
//     );
//      */
//     let mut wgt_boost = WgtBoost::construct_wgt(wgt, iter, loss, table2_sum, is_eps, None);

//     unimplemented!("Why using Eps::Med?->introduce eff_eps");
//     // no learning rate
//     let coef_ti = coefficient::calculate_coef_ada_eps(
//         wgt_boost.contingency_table().unwrap(),
//         BoostType::Ada,
//         1.0,
//         None, //Some(crate::Eps::Med), // will not be used
//     );

//     wgt_boost.wgt_mut().model_mut().set_coef(coef_ti);
//     //wgt.model_mut().set_coef_binary(coef_ti);

//     wgt_boost

//     /*
//     if let WgtKind::Cov(ref mut cov_wgt) = &mut (wgt.get_kind_mut()) {
//         if let Model::Binary(ref mut binary_model) = &mut (cov_wgt.get_model_mut()) {
//             binary_model.set_coef(coef_ti);
//         }
//     }
//     */
//     /*
//     if let WgtKind::Snv(ref mut snv_wgt) = &mut (wgt.get_kind_mut()) {
//         snv_wgt.get_model_mut().set_coef(coef_ti);
//         //let model = snv_wgt.get_model_mut();
//         //model.set_coef(coef_ti);
//         //(*(*snv_wgt).get_model()).set_coef(coef_ti);
//         //(*(*snv_wgt).get_model()).set_coef(coef_ti);
//     } else {
//         panic!("wgt should be Snv.")
//     }
//     */
//     //// TODO: make these func
//     //compute::compute_mis(&mut mis, &wgt, ys, mistakes, n);
//     //let abcd_sum = calculate::calculate_abcd_sum();

//     //wgt
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resize_interaction() {
        let m = 3;
        let m_inter = 5;
        let mut loss_inter =
            LossStruct::new_interaction(BoostType::LogitAddInteraction, m, m_inter);

        loss_inter.resize_interaction(10);

        if let LossStruct::Interaction(_, loss) = loss_inter {
            assert_eq!(loss.len(), 10);
        }
    }

    #[test]
    fn test_search_min() {
        let loss = vec![0.1, 0.5, 0.2, 0.4, 0.3];
        let loss_inter = LossStruct::new_vec(BoostType::Logit, loss);

        // [0.5, 0.2, 0.3]
        let extract = vec![1usize, 2, 4];
        let extract: HashSet<usize> = HashSet::from_iter(extract.into_iter());

        let (loss_min, mi_ans, si_ans) = loss_inter.search_min(Some(&extract));
        assert_eq!(loss_min, 0.2);
        assert_eq!(mi_ans, Some(2));
        assert_eq!(si_ans, None);
    }

    #[test]
    fn test_search_min_interaction() {
        let loss = vec![0.2, 0.5, 0.1, 0.4, 0.3];
        let loss_inter =
            LossStruct::new_interaction_vec(BoostType::LogitAddInteraction, vec![], loss);

        let extract = vec![(0usize, 1usize), (1, 2), (2, 3), (3, 4), (4, 5)];

        let (loss_min, pair_1, pair_2) = loss_inter.search_min_interaction(&extract);
        assert_eq!(loss_min, 0.1);
        assert_eq!(pair_1, Some(2));
        assert_eq!(pair_2, Some(3));
    }

    #[test]
    fn test_search_topprop_n() {
        let loss = vec![0.1, 0.5, 0.2, 0.4, 0.3];
        let loss_inter = LossStruct::new_vec(BoostType::Logit, loss);

        // [0.5, 0.2,0.4, 0.3]
        let extract = vec![1usize, 2, 3, 4];
        let extract: HashSet<usize> = HashSet::from_iter(extract.into_iter());

        let top_n = 2;

        let (use_snvs, loss_next) = loss_inter.search_topprop_n(top_n, Some(&extract));
        assert_eq!(use_snvs, vec![false, false, true, false, true]);
        assert_eq!(loss_next, 0.4);
    }

    #[test]
    fn test_search_topprop_n_interaction() {
        let loss = vec![0.1, 0.5, 0.2, 0.4, 0.3];
        let loss_inter =
            LossStruct::new_interaction_vec(BoostType::LogitAddInteraction, vec![], loss);

        // [0.5, 0.2,0.4, 0.3]
        let top_n = 2;

        let (use_snvs, loss_next) = loss_inter.search_topprop_n_interaction(top_n);
        assert_eq!(use_snvs, vec![true, false, true, false, false]);
        assert_eq!(loss_next, 0.3);
    }
}
