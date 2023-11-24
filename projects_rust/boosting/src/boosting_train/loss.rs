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
use genetics::samples::prelude::*;
use genetics::vec::sort_float;
use genetics::wgt::Wgt;
use genetics::{Genot, Snvs};
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
    // loss, m
    //Logit(Vec<f64>, usize),
}

impl LossStruct {
    pub fn new(boost_type: BoostType, m: usize) -> Self {
        match boost_type {
            BoostType::ConstAda => LossStruct::ConstAda(vec![f64::NAN; m * 2], m, 2),
            BoostType::FreeModelMissing
            | BoostType::Logit
            | BoostType::LogitNoMissing
            | BoostType::LogitAdd => LossStruct::LossOne(vec![f64::NAN; m], m),
            BoostType::Ada => unimplemented!(),
        }
    }

    pub fn new_vec(boost_type: BoostType, vec: Vec<f64>) -> Self {
        let m = vec.len();
        match boost_type {
            BoostType::ConstAda => LossStruct::ConstAda(vec, m / 2, 2),
            BoostType::FreeModelMissing
            | BoostType::Logit
            | BoostType::LogitNoMissing
            | BoostType::LogitAdd => LossStruct::LossOne(vec, m),
            BoostType::Ada => unimplemented!(),
        }
    }

    //pub fn search_min(&self, skip_snv: &[usize]) -> (f64, Option<usize>, Option<usize>) {
    pub fn search_min(&self, skip_snv: &HashSet<usize>) -> (f64, Option<usize>, Option<usize>) {
        match self {
            // FIXME: implement skip_snv for constada
            LossStruct::ConstAda(..) => self.search_min_constada(),
            LossStruct::LossOne(..) => self.search_min_modelfree(skip_snv),
            //LossStruct::Logit(..) => self.search_min_logit(skip_snv),
        }
    }
    pub fn inner_mut(&mut self) -> &mut [f64] {
        // TODO: all right?
        match self {
            LossStruct::ConstAda(ref mut v, _, _) => v,
            LossStruct::LossOne(ref mut v, _) => v,
        }
    }

    pub fn inner(&self) -> &[f64] {
        // TODO: all right?
        match self {
            LossStruct::ConstAda(ref v, _, _) => v,
            LossStruct::LossOne(ref v, _) => v,
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
    fn search_min_constada(&self) -> (f64, Option<usize>, Option<usize>) {
        if let LossStruct::ConstAda(loss, _, _) = self {
            let mut loss_min = f64::INFINITY;
            let mut var_si: Option<usize> = None;
            let mut var_mi: Option<usize> = None;

            for (i, loss) in loss.iter().enumerate() {
                if *loss < loss_min {
                    loss_min = *loss;
                    let (var_si_, var_mi_) = self.index_constada(i);
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

    fn search_min_modelfree(
        &self,
        skip_snv: &HashSet<usize>,
    ) -> (f64, Option<usize>, Option<usize>) {
        if let LossStruct::LossOne(loss, _) = self {
            // for Logit
            let mut loss_min = f64::INFINITY;
            let mut var_mi: Option<usize> = None;

            //log::debug!("loss m {}", loss.len());
            for (mi, loss) in loss.iter().enumerate() {
                //if mi % 100000 == 0 {
                //    log::debug!("min loss mi {} {}", mi, *loss);
                //}

                if (*loss < loss_min) && !(skip_snv.contains(&mi)) {
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

            log::debug!(
                "Selected index, mi, loss: {}, {:.6}",
                var_mi.unwrap(),
                loss_min
            );

            return (loss_min, var_mi, None);
        } else {
            panic!("Sth wrong.")
        }
    }

    pub fn search_topprop(&self, prop: f64) -> Vec<bool> {
        match self {
            LossStruct::ConstAda(..) => unimplemented!(),
            LossStruct::LossOne(..) => self.search_topprop_modelfree(prop),
        }
    }

    // FIXME: make search_topn() like search_topprop()
    /// also return next f64
    pub fn search_topprop_modelfree_n(
        &self,
        m_top: usize,
        skip_snv: Option<&HashSet<usize>>,
    ) -> (Vec<bool>, f64) {
        if let LossStruct::LossOne(loss, _) = self {
            //let mut loss_topprop = 10.0;

            let m = loss.len();

            //log::debug!("loss {:?}", loss);
            log::debug!("loss any {:?}", loss.iter().any(|x| x.is_nan()));
            log::debug!(
                "loss no_nan {:?}",
                loss.iter()
                    .filter(|x| x.is_nan())
                    .map(|x| *x)
                    .collect::<Vec<f64>>()
            );

            let mut loss_sort = loss.clone();

            loss_sort.iter_mut().for_each(|x| {
                if x.is_nan() {
                    *x = f64::MAX
                }
            });

            // make loss of skip_snv f64::MAX
            if let Some(skip_snv) = skip_snv {
                skip_snv.iter().for_each(|i| loss_sort[*i] = f64::MAX);
            }
            sort_float(&mut loss_sort);

            // takes at least 1 SNVs
            //let topprop_n = ((m as f64 * prop) as usize).max(1);

            // -2 for loss_topprop_top_outside
            let topprop_n = m_top.min(m) - 2;

            // the smaller, the better
            let loss_topprop = loss_sort[topprop_n];
            let loss_topprop_top_outside = loss_sort[topprop_n + 1];

            // FIXME: if loss_topprop, loss_topprop_ == f64::MAX

            let mut use_snvs = vec![false; m];

            if let Some(skip_snv) = skip_snv {
                for (mi, loss) in loss.iter().enumerate() {
                    //if *loss <= loss_topprop {
                    if (*loss <= loss_topprop) && !(skip_snv.contains(&mi)) {
                        use_snvs[mi] = true;
                    }
                }
            } else {
                for (mi, loss) in loss.iter().enumerate() {
                    if *loss <= loss_topprop {
                        use_snvs[mi] = true;
                    }
                }
            }

            log::debug!(
                "#SNVs in use_snvs by loss {} in {}",
                use_snvs.iter().filter(|b| **b).count(),
                m
            );

            return (use_snvs, loss_topprop_top_outside);
            //return (loss_min, var_mi, var_si);
        } else {
            panic!("Sth wrong.")
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

    // FIXME constada
    // modelfree only now
    fn string_write(&self, snvs: &Snvs) -> String {
        if let LossStruct::LossOne(loss, _) = self {
            log::debug!("loss");
            log::debug!("{}", loss[0]);
            log::debug!("{}", loss[1]);
            let strings = loss
                .iter()
                .zip(snvs.snv_indexs().iter())
                .enumerate()
                .map(|(mi, (l, snv_id))| {
                    mi.to_string()
                        + "\t"
                        + snv_id.rs()
                        + "\t"
                        + &snv_id.chrom().to_string()
                        + "\t"
                        + &snv_id.pos().to_string()
                        + "\t"
                        + snv_id.a1()
                        + "\t"
                        + snv_id.a2()
                        + "\t"
                        + &l.to_string()
                })
                .collect::<Vec<String>>();
            let str = "m\trs\tchrom\tpos\tA1\tA2\tloss".to_owned();
            str + "\n" + &strings.join("\n")
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
        writer.write(str.as_bytes()).unwrap();
        // capture error
        writer.flush().unwrap();
    }
}

pub fn search_min_loss_gt(
    loss: &mut LossStruct,
    iteration: usize,
    genot: &Genot,
    sample_weight: &SampleWeight,
    phe: &Phe,
    snvs: &Snvs,
    boost_param: BoostParam,
    skip_snv: &HashSet<usize>,
    //skip_snv: &[usize],
    use_adjloss: bool,
) -> WgtBoost {
    // calc loss
    calculate_loss_gt(
        loss,
        genot,
        sample_weight,
        phe,
        boost_param,
        skip_snv,
        use_adjloss,
    );

    // search min loss
    let (loss_min, var_mi, var_si) = loss.search_min(skip_snv);
    let var_mi: usize = var_mi.unwrap_or_else(|| panic!("No SNVs were loss<1.0"));

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
            let wgt_min = Wgt::construct_snv_index(
                snvs.snv_indexs()[var_mi].clone(),
                genetics::THRESHOLD_SNV[var_si],
                (var_si, var_mi),
            );

            let wgt_boost_min = WgtBoost::construct_wgt_loss(wgt_min, iteration, Some(loss_min));
            wgt_boost_min
        }

        BoostType::FreeModelMissing
        | BoostType::Logit
        | BoostType::LogitNoMissing
        | BoostType::LogitAdd => {
            let wgt_min =
                Wgt::construct_snv_index_freemodel(snvs.snv_indexs()[var_mi].clone(), var_mi);
            let wgt_boost_min = WgtBoost::construct_wgt_loss(wgt_min, iteration, Some(loss_min));
            wgt_boost_min
        }

        BoostType::Ada => {
            unimplemented!()
        }
    }
}

pub fn calculate_loss_gt(
    loss_gt: &mut LossStruct,
    genot: &Genot,
    sample_weight: &SampleWeight,
    phe: &Phe,
    boost_param: BoostParam,
    skip_snv: &HashSet<usize>,
    use_adjloss: bool,
) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            log::debug!("Use SIMD.");
            unsafe {
                match boost_param.boost_type() {
                    BoostType::Logit | BoostType::LogitNoMissing => {
                        calc::calculate_loss_gt_logit_simd(
                            loss_gt,
                            genot,
                            sample_weight,
                            //sample_weight.wzs_pad().unwrap(),
                            //sample_weight.wls_pad().unwrap(),
                            //sample_weight.zs_pad().unwrap(),
                            phe,
                            boost_param,
                            skip_snv,
                            use_adjloss,
                        );
                    }
                    BoostType::LogitAdd => {
                        calc::calculate_loss_gt_logit_add_simd(
                            loss_gt,
                            genot,
                            sample_weight,
                            //sample_weight.wzs_pad().unwrap(),
                            //sample_weight.wls_pad().unwrap(),
                            //sample_weight.zs_pad().unwrap(),
                            //phe,
                            boost_param,
                            skip_snv,
                            use_adjloss,
                        );
                    }
                    BoostType::ConstAda => {
                        calc::calculate_loss_gt_constada_simd(
                            loss_gt,
                            genot,
                            sample_weight,
                            //sample_weight.ps_pad().unwrap(),
                            phe,
                            boost_param,
                            skip_snv,
                        );
                    }
                    BoostType::FreeModelMissing => {
                        calc::calculate_loss_gt_freemodelmissing_simd(
                            loss_gt,
                            genot,
                            sample_weight,
                            //sample_weight.ps_pad().unwrap(),
                            phe,
                            boost_param,
                            skip_snv,
                        );
                    }
                    BoostType::Ada => {
                        panic!("");
                    }
                };
            }
            return;
        }
    }

    log::debug!("Do not use SIMD.");

    match boost_param.boost_type() {
        BoostType::ConstAda => {
            calc::calculate_loss_gt_constada_nosimd(
                loss_gt,
                genot,
                sample_weight,
                //sample_weight.ps().unwrap(),
                phe,
                boost_param,
            );
        }
        BoostType::FreeModelMissing => {
            calc::calculate_loss_gt_freemodelmissing_nosimd(
                loss_gt,
                genot,
                sample_weight,
                //sample_weight.ps().unwrap(),
                phe,
                boost_param,
            );
        }
        // FIXME: BoostType::LogitNoMissing
        _ => unimplemented!("Use SIMD."),
    }
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
    boost_param: BoostParam,
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
            let wgt_min = Wgt::construct_snv_index(
                snvs.snv_indexs()[var_mi].clone(),
                genetics::THRESHOLD_SNV[var_si],
                (var_si, var_mi),
            );

            let wgt_boost_min = WgtBoost::construct_wgt_loss(wgt_min, iteration, Some(loss_second));
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

/*
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_min_loss_gt() {}
}
*/
