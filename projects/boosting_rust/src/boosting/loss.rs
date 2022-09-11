// pub for bench
pub mod calc;

use super::coefficient;
use super::compute_pred;
use super::epsilon;
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

#[derive(Debug)]
pub enum LossStruct {
    // loss, m, s=2
    // loss[mi, si]
    // (loss[0,0], loss[0,1], loss[1,0],loss[1,1]...)
    ConstAda(Vec<f64>, usize, usize),
    // loss, m
    ModelFree(Vec<f64>, usize),
}

impl LossStruct {
    pub fn new(boost_type: BoostType, m: usize) -> Self {
        match boost_type {
            BoostType::ConstAda => LossStruct::ConstAda(vec![f64::NAN; m * 2], m, 2),
            BoostType::FreeModelMissing => LossStruct::ModelFree(vec![f64::NAN; m], m),
            _ => unimplemented!(),
        }
    }

    pub fn search_min(&self) -> (f64, Option<usize>, Option<usize>) {
        match self {
            LossStruct::ConstAda(..) => self.search_min_constada(),
            LossStruct::ModelFree(..) => self.search_min_modelfree(),
        }
    }
    pub fn inner(&mut self) -> &mut [f64] {
        // TODO: all right?
        match self {
            LossStruct::ConstAda(ref mut v, _, _) => v,
            LossStruct::ModelFree(ref mut v, _) => v,
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
            // make much larger than 1.0 since could be >1.0 by eps
            let mut loss_min = 10.0;
            //let mut loss_min = 1.0;
            let mut var_si: Option<usize> = None;
            let mut var_mi: Option<usize> = None;

            for (i, loss) in loss.iter().enumerate() {
                if *loss < loss_min {
                    loss_min = *loss;
                    let (var_si_, var_mi_) = self.index_constada(i);
                    var_si = Some(var_si_);
                    var_mi = Some(var_mi_);
                    println!(
                        "Temporary chosen index, si, mi, loss: {}, {}, {}, {:.6}",
                        i,
                        var_si.unwrap(),
                        var_mi.unwrap(),
                        *loss
                    );
                }
            }
            return (loss_min, var_mi, var_si);
        } else {
            panic!("Sth wrong.")
        }
    }

    fn search_min_modelfree(&self) -> (f64, Option<usize>, Option<usize>) {
        if let LossStruct::ModelFree(loss, _) = self {
            // make much larger than 1.0 since could be >1.0 by eps
            let mut loss_min = 10.0;
            //let var_si: Option<usize> = None;
            let mut var_mi: Option<usize> = None;

            //println!("loss m {}", loss.len());
            for (mi, loss) in loss.iter().enumerate() {
                //if mi % 100000 == 0 {
                //    println!("min loss mi {} {}", mi, *loss);
                //}
                if *loss < loss_min {
                    loss_min = *loss;
                    let var_mi_ = mi;
                    //let (var_si_, var_mi_) = predict_index_from_loss(i);
                    //var_si = Some(var_si_);
                    var_mi = Some(var_mi_);
                    /*
                    var_si = Some(i % 2);
                    var_mi = Some(i / 2);
                     */
                    //var_si = i / m;
                    //var_mi = i % m;
                    println!(
                        "Temporary chosen index, mi, loss: {}, {:.6}",
                        //var_si.unwrap(),
                        var_mi.unwrap(),
                        *loss
                    );
                }
            }

            return (loss_min, var_mi, None);
            //return (loss_min, var_mi, var_si);
        } else {
            panic!("Sth wrong.")
        }
    }

    pub fn search_topprop(&self, prop: f64) -> Vec<bool> {
        match self {
            LossStruct::ConstAda(..) => unimplemented!(),
            LossStruct::ModelFree(..) => self.search_topprop_modelfree(prop),
        }
    }

    fn search_topprop_modelfree(&self, prop: f64) -> Vec<bool> {
        if let LossStruct::ModelFree(loss, _) = self {
            //let mut loss_topprop = 10.0;

            let m = loss.len();

            let mut loss_sort = loss.clone();
            sort_float(&mut loss_sort);

            // takes at least 1 SNVs
            let topprop_n = ((m as f64 * prop) as usize).max(1);
            // the smaller, the better
            let loss_topprop = loss_sort[topprop_n];

            let mut use_snvs = vec![false; m];

            for (mi, loss) in loss.iter().enumerate() {
                if *loss <= loss_topprop {
                    use_snvs[mi] = true;
                }
            }

            println!(
                "use_snvs by loss {} in {}",
                use_snvs.iter().filter(|b| **b).count(),
                m
            );

            return use_snvs;
            //return (loss_min, var_mi, var_si);
        } else {
            panic!("Sth wrong.")
        }
    }

    // FIXME constada
    // modelfree only now
    fn string_write(&self, snvs: &Snvs) -> String {
        if let LossStruct::ModelFree(loss, _) = self {
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
        //println!("str wgt {}", &str);
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
    ps: &[f64],
    phe: &Phe,
    snvs: &Snvs,
    boost_param: BoostParam,
) -> WgtBoost {
    // calc loss
    calculate_loss_gt(loss, genot, ps, phe, boost_param);

    // search min loss
    let (loss_min, var_mi, var_si) = loss.search_min();
    let var_mi: usize = var_mi.expect("No SNVs were loss<1.0");

    println!(
        "contingency table count {:?}",
        genot.to_genot_snv(var_mi).stat_contingency_table(phe)
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
        BoostType::FreeModelMissing => {
            let wgt_min =
                Wgt::construct_snv_index_freemodel(snvs.snv_indexs()[var_mi].clone(), var_mi);
            let wgt_boost_min = WgtBoost::construct_wgt_loss(wgt_min, iteration, Some(loss_min));
            wgt_boost_min
        }
        _ => {
            panic!()
        }
    }
}

pub fn calculate_loss_gt(
    loss_gt: &mut LossStruct,
    genot: &Genot,
    ps: &[f64],
    phe: &Phe,
    boost_param: BoostParam,
) {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            println!("Use SIMD.");
            //panic!("simd results in bug");
            unsafe {
                match boost_param.boost_type() {
                    BoostType::ConstAda => {
                        calc::calculate_loss_gt_constada_simd(loss_gt, genot, ps, phe, boost_param);
                    }
                    BoostType::FreeModelMissing => {
                        calc::calculate_loss_gt_freemodelmissing_simd(
                            loss_gt,
                            genot,
                            ps,
                            phe,
                            boost_param,
                        );
                    }
                    _ => {}
                };
            }
            return;
        }
    }

    println!("Do not use SIMD.");

    match boost_param.boost_type() {
        BoostType::ConstAda => {
            calc::calculate_loss_gt_constada_nosimd(loss_gt, genot, ps, phe, boost_param);
        }
        BoostType::FreeModelMissing => {
            calc::calculate_loss_gt_freemodelmissing_nosimd(loss_gt, genot, ps, phe, boost_param);
        }
        _ => {}
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

    //println!("loss_gt {:?}", loss_gt);

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
            println!(
                "Temporary chosen index, si, mi, loss: {}, {}, {}, {:.6}",
                i,
                var_si.unwrap(),
                var_mi.unwrap(),
                *loss
            );
        }
    }

    println!(
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
            println!(
                "Temporary chosen mi, loss: {}, {:.6}",
                var_mi.unwrap(),
                loss
            );
        }
    }

    println!("Chosen mi, loss: {}, {:.6}", var_mi.unwrap(), loss_min);

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

    println!("Chosen mi, loss: {}, {:.6}", var_mi, loss_min);

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

pub fn create_loss_const(
    iter: usize,
    ps: &[f64],
    pred: &mut [u8],
    //pred: &mut [B8],
    ys: &Phe,
    _n: usize,
    _boosting_param: BoostParam,
) -> WgtBoost {
    /*
    let mut const_var = Var::construct_var(CovKind::Const, "const".to_owned());
    const_var.set_vals_const(n);
    let const_wgt = CovWgt::construct(const_var);
     */

    compute_pred::compute_pred_const(pred);

    //let abcd_sum = calculate::calculate_abcd_sum(
    let (ab_sum, is_eps) = table::calculate_table2_sum(pred, ps, ys);

    let loss = calc::calculate_loss_ab(ab_sum);

    // since all predict=1
    let wgt = Wgt::construct_const_threshold(0.5);

    /*
    let wgt = Wgt::construct_wgt(
        WgtKind::Cov(const_wgt),
        Model::Binary(BinaryModel::construct_threshold(0.5)),
    );
     */
    let mut wgt_boost = WgtBoost::construct_wgt(wgt, iter, loss, ab_sum, is_eps);

    // no learning rate
    let coef_ti = coefficient::calculate_coefficients(
        wgt_boost.contingency_table().unwrap(),
        BoostType::Ada,
        None,
    );

    wgt_boost.wgt_mut().model_mut().set_coef(coef_ti);
    //wgt.model_mut().set_coef_binary(coef_ti);

    wgt_boost

    /*
    if let WgtKind::Cov(ref mut cov_wgt) = &mut (wgt.get_kind_mut()) {
        if let Model::Binary(ref mut binary_model) = &mut (cov_wgt.get_model_mut()) {
            binary_model.set_coef(coef_ti);
        }
    }
    */

    /*
    if let WgtKind::Snv(ref mut snv_wgt) = &mut (wgt.get_kind_mut()) {
        snv_wgt.get_model_mut().set_coef(coef_ti);
        //let model = snv_wgt.get_model_mut();
        //model.set_coef(coef_ti);
        //(*(*snv_wgt).get_model()).set_coef(coef_ti);
        //(*(*snv_wgt).get_model()).set_coef(coef_ti);
    } else {
        panic!("wgt should be Snv.")
    }
    */

    //// TODO: make these func
    //compute::compute_mis(&mut mis, &wgt, ys, mistakes, n);
    //let abcd_sum = calculate::calculate_abcd_sum();

    //wgt
}

/*
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_min_loss_gt() {}
}
*/
