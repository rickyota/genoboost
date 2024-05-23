use super::sample_weight::SampleWeight;
use crate::{boosting_train::coefficient, wgt_boost::WgtBoost, BoostParam};
use genetics::{samples::prelude::Phe, BaseGenot, Dataset, Genot, WgtKind};

pub fn update_coef(
    wgt: &mut WgtBoost,
    boost_param: &BoostParam,
    //pred: &mut [u8],
    _phe: &Phe,
    genot: &Genot,
    sample_weight: &SampleWeight,
    dataset: &Dataset,
) {
    if boost_param.boost_type().is_type_ada() {
        //let ps_pad = sample_weight.ps_pad().unwrap();
        /*
        // create table
        let (table_sum, is_eps) =
            table::calculate_table_eps(&pred, ps, phe, boost_param.eps(), boost_param.boost_type());
        log::debug!("table {:?}", table_sum);
        wgt.set_contingency_table(table_sum, is_eps);

        // TODO: just end without error when loss==1.0

        let coef_ti = coefficient::calculate_coefficients(
            wgt.contingency_table().unwrap(),
            boost_param.boost_type(),
            boost_param.learning_rate(),
            boost_param.eps(),
        );
        log::debug!("coef {:?}", coef_ti);
        */

        unimplemented!()

        //let mi: usize = match wgt.wgt().kind() {
        //    WgtKind::Snv(_, _, mi) => mi.unwrap(),
        //    _ => panic!(),
        //};

        // TO implement not using pred
        //// not set cont table
        ////wgt.set_contingency_table(table_sum, is_eps);
        //let (coef_ti, is_eps, is_eff_eps) = coefficient::calculate_coef_ada_update(
        //    &pred,
        //    &genot.to_genot_snv(mi),
        //    sample_weight,
        //    //ps_pad,
        //    phe,
        //    boost_param.learning_rate(),
        //    boost_param.eps(),
        //    boost_param.eff_eps(),
        //    boost_param.boost_type(),
        //);
        //log::debug!("coef {:?}", coef_ti);
        //wgt.set_coef(coef_ti);
        //wgt.set_is_eps(is_eps);
        //wgt.set_is_eff_eps(is_eff_eps);
        //log::debug!("wgt {:?}", wgt);
    } else if boost_param.boost_type().is_type_logit() {
        //let mi: usize = match wgt.wgt().kind() {
        //    WgtKind::Snv(_, _, mi) => mi.unwrap(),
        //    _ => panic!(),
        //};

        match wgt.wgt().kind() {
            //WgtKind::Snv(_, _, mi) => {
            WgtKind::Snv(snv_wgt) => {
                let mi = snv_wgt.index().unwrap();
                //let mi = mi.unwrap();

                // TOFIX: pass Coef:: type here, but now wgt here does not have Coef
                // add coef in loss?
                // TODO: cleaner
                let (coef_ti, is_eps, is_eff_eps) = coefficient::calculate_coef_logit_update(
                    &genot.to_genot_snv(mi),
                    mi,
                    dataset,
                    sample_weight,
                    boost_param.learning_rate(),
                    //boost_param.eps(),
                    //boost_param.eff_eps(),
                    boost_param.boost_type(),
                    boost_param,
                );

                log::debug!("coef {:?}", coef_ti);
                wgt.set_coef(coef_ti);
                wgt.set_is_eps(is_eps);
                wgt.set_is_eff_eps(is_eff_eps);

                log::debug!("wgt {:?}", wgt);
            }
            //WgtKind::SnvInteraction(_, mi_1, _, mi_2) => {
            WgtKind::SnvInteraction(snv_inter_wgt) => {
                let (mi_1, mi_2) = snv_inter_wgt.indexs();
                let mi_1 = mi_1.unwrap();
                let mi_2 = mi_2.unwrap();

                let mafs = dataset.snvs().mafs().unwrap();

                let (coef_ti, is_eps, is_eff_eps) =
                    coefficient::calculate_coef_logit_interaction_update(
                        &genot.to_genot_snv(mi_1),
                        &genot.to_genot_snv(mi_2),
                        mi_1,
                        mi_2,
                        dataset,
                        sample_weight,
                        boost_param.learning_rate(),
                        boost_param.boost_type(),
                        boost_param,
                        mafs[mi_1],
                        mafs[mi_2],
                    );

                log::debug!("coef {:?}", coef_ti);
                wgt.set_coef(coef_ti);
                wgt.set_is_eps(is_eps);
                wgt.set_is_eff_eps(is_eff_eps);

                log::debug!("wgt {:?}", wgt);
            }
            _ => panic!(),
        }
    };
}
