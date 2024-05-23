use super::LossFunc;
//use super::WgtBoost;
use crate::boosting_param::{SampleWeightClip, SampleWeightWlsClip};
use crate::BoostType;
use genetics::samples::prelude::*;
use genetics::vec;
//use genetics::wgt::Coef;
use genetics::{alloc, SampleScore};

// add score or not? -> no, since score can be used without wgt
// NEVER use clone() since vec is not aligned -> use clone_align()
#[derive(Default, Debug)]
pub struct SampleWeight {
    n: usize,
    len_n: usize,
    //score: Vec<f64>,
    // Ada
    ws: Option<Vec<f64>>, //len=n
    ps_pad: Option<Vec<f64>>,
    // Logit
    probs: Option<Vec<f64>>,   //len=n
    zs_pad: Option<Vec<f64>>,  // z
    wzs_pad: Option<Vec<f64>>, // w*z
    wls_pad: Option<Vec<f64>>, // w for LogitBoost
    // TOFIX: remove config??
    // config
    boost_type: BoostType,
    loss_func: LossFunc,
    // for ws and wzs // TODO: split
    sample_weight_clip: Option<SampleWeightClip>,
    sample_weight_wls_clip: Option<SampleWeightWlsClip>,
}

// padding <-> supress
fn supress_len(v: Option<&[f64]>, n: usize) -> Option<&[f64]> {
    v.map(|x| &x[..n])
}
fn supress_len_mut(v: Option<&mut [f64]>, n: usize) -> Option<&mut [f64]> {
    //v.map(|x| &x[..self.n])
    //let ws_mut=self.ws.as_deref_mut();
    v.map(|x| &mut x[..n])
}

impl SampleWeight {
    pub fn n(&self) -> usize {
        self.n
    }
    // len_n >= n
    pub fn len_n(&self) -> usize {
        self.len_n
    }

    // idealy, this should be ws() but usually supressed len is used
    //pub fn ws_pad(&self) -> Option<&[f64]> {
    //    self.ws_pad.as_deref()
    //}
    pub fn ws(&self) -> Option<&[f64]> {
        //ws is n
        self.ws.as_deref()
        //self.ws().map(|x| &x[..self.n])
        //supress_len(self.ws_pad(), self.n)
    }
    //pub fn ws_pad_mut(&mut self) -> Option<&mut [f64]> {
    //    self.ws_pad.as_deref_mut()
    //}
    pub fn ws_mut(&mut self) -> Option<&mut [f64]> {
        self.ws.as_deref_mut()
        //let n = self.n;
        //supress_len_mut(self.ws_pad_mut(), n)
    }

    pub fn ps_pad(&self) -> Option<&[f64]> {
        self.ps_pad.as_deref()
    }
    pub fn ps_pad_mut(&mut self) -> Option<&mut [f64]> {
        self.ps_pad.as_deref_mut()
    }
    pub fn ps(&self) -> Option<&[f64]> {
        supress_len(self.ps_pad(), self.n)
    }
    pub fn ps_mut(&mut self) -> Option<&mut [f64]> {
        let n = self.n;
        supress_len_mut(self.ps_pad_mut(), n)
    }

    //pub fn probs_pad(&self) -> Option<&[f64]> {
    //    self.probs_pad.as_deref()
    //}
    //pub fn probs_pad_mut(&mut self) -> Option<&mut [f64]> {
    //    self.probs.as_deref_mut()
    //}
    pub fn probs(&self) -> Option<&[f64]> {
        self.probs.as_deref()
        //supress_len(self.probs_pad(), self.n)
    }
    pub fn probs_mut(&mut self) -> Option<&mut [f64]> {
        self.probs.as_deref_mut()
        //let n = self.n;
        //supress_len_mut(self.probs_pad_mut(), n)
    }

    pub fn zs_pad(&self) -> Option<&[f64]> {
        self.zs_pad.as_deref()
    }
    pub fn zs_pad_mut(&mut self) -> Option<&mut [f64]> {
        self.zs_pad.as_deref_mut()
    }
    pub fn zs(&self) -> Option<&[f64]> {
        supress_len(self.zs_pad(), self.n)
    }
    pub fn zs_mut(&mut self) -> Option<&mut [f64]> {
        let n = self.n;
        supress_len_mut(self.zs_pad_mut(), n)
    }

    pub fn wzs_pad(&self) -> Option<&[f64]> {
        self.wzs_pad.as_deref()
    }
    pub fn wzs_pad_mut(&mut self) -> Option<&mut [f64]> {
        self.wzs_pad.as_deref_mut()
    }
    pub fn wzs(&self) -> Option<&[f64]> {
        supress_len(self.wzs_pad(), self.n)
    }
    pub fn wzs_mut(&mut self) -> Option<&mut [f64]> {
        let n = self.n;
        supress_len_mut(self.wzs_pad_mut(), n)
    }

    pub fn wls_pad(&self) -> Option<&[f64]> {
        self.wls_pad.as_deref()
    }
    pub fn wls_pad_mut(&mut self) -> Option<&mut [f64]> {
        self.wls_pad.as_deref_mut()
    }
    pub fn wls(&self) -> Option<&[f64]> {
        supress_len(self.wls_pad(), self.n)
    }
    pub fn wls_mut(&mut self) -> Option<&mut [f64]> {
        let n = self.n;
        supress_len_mut(self.wls_pad_mut(), n)
    }

    pub fn clear_pad(&mut self) {
        let n = self.n();

        if let Some(x) = self.ps_pad_mut() {
            x[n..].iter_mut().for_each(|x| *x = 0.0f64)
        };
        if let Some(x) = self.zs_pad_mut() {
            x[n..].iter_mut().for_each(|x| *x = 0.0f64)
        };
        if let Some(x) = self.wzs_pad_mut() {
            x[n..].iter_mut().for_each(|x| *x = 0.0f64)
        };
        if let Some(x) = self.wls_pad_mut() {
            x[n..].iter_mut().for_each(|x| *x = 0.0f64)
        };

        // self.ps_pad_mut().unwrap()[n..]
        //     .iter_mut()
        //     .for_each(|x| *x = 0.0f64);
        // self.zs_pad_mut().unwrap()[n..]
        //     .iter_mut()
        //     .for_each(|x| *x = 0.0f64);
        // self.wzs_pad_mut().unwrap()[n..]
        //     .iter_mut()
        //     .for_each(|x| *x = 0.0f64);
        // self.wls_pad_mut().unwrap()[n..]
        //     .iter_mut()
        //     .for_each(|x| *x = 0.0f64);
    }

    // TODO: add check len
    pub fn new(
        n: usize,
        boost_type: BoostType,
        loss_func: LossFunc,
        sample_weight_clip: Option<SampleWeightClip>,
        sample_weight_wls_clip: Option<SampleWeightWlsClip>,
    ) -> SampleWeight {
        if boost_type.is_type_ada() {
            let ws: Vec<f64> = vec![1.0 / (n as f64); n];
            let len_n = n + 32;
            let mut ps_pad: Vec<f64> = alloc::with_capacity_align::<f64>(len_n);
            ps_pad.resize(n + 32, 0.0f64);
            SampleWeight {
                n,
                len_n,
                ws: Some(ws),
                ps_pad: Some(ps_pad),
                boost_type,
                loss_func,
                sample_weight_clip,
                sample_weight_wls_clip,
                ..Default::default() // might cause bug
            }
        } else if boost_type.is_type_logit() {
            let probs: Vec<f64> = vec![0.0f64; n];
            //let mut zs: Vec<f64> = vec![0.0f64; n];
            // len != n only for wzs, wls
            // last 32=0.0
            let len_n = n + 32;
            let mut zs_pad: Vec<f64> = alloc::with_capacity_align::<f64>(len_n);
            zs_pad.resize(len_n, 0.0f64);
            let mut wzs_pad: Vec<f64> = alloc::with_capacity_align::<f64>(len_n);
            wzs_pad.resize(len_n, 0.0f64);
            let mut wls_pad: Vec<f64> = alloc::with_capacity_align::<f64>(len_n);
            wls_pad.resize(len_n, 0.0f64);
            SampleWeight {
                n,
                len_n,
                probs: Some(probs),
                zs_pad: Some(zs_pad),
                wzs_pad: Some(wzs_pad),
                wls_pad: Some(wls_pad),
                boost_type,
                loss_func,
                sample_weight_clip,
                sample_weight_wls_clip,
                ..Default::default()
            }
        } else {
            panic!();
        }
    }

    // TODO: implement as Clone trait?
    pub fn clone_align(&self) -> Self {
        //let len_n = self.len_n;
        let ps_pad = self.ps_pad().map(|x| {
            let mut y: Vec<f64> = alloc::with_capacity_align::<f64>(x.len());
            y.resize(x.len(), 0.0f64);
            y.copy_from_slice(x);
            y
        });
        let zs_pad = self.zs_pad().map(|x| {
            // len=n samples
            //log::debug!("zs_pad len {}", x.len());
            let mut y: Vec<f64> = alloc::with_capacity_align::<f64>(x.len());
            // len=0 and capacity=n samples
            //log::debug!("y len {}", y.len());
            y.resize(x.len(), 0.0f64);
            y.copy_from_slice(x);
            y
        });
        let wzs_pad = self.wzs_pad().map(|x| {
            let mut y: Vec<f64> = alloc::with_capacity_align::<f64>(x.len());
            y.resize(x.len(), 0.0f64);
            y.copy_from_slice(x);
            y
        });
        let wls_pad = self.wls_pad().map(|x| {
            let mut y: Vec<f64> = alloc::with_capacity_align::<f64>(x.len());
            y.resize(x.len(), 0.0f64);
            y.copy_from_slice(x);
            y
        });

        SampleWeight {
            n: self.n,
            len_n: self.len_n,
            ws: self.ws.clone(),
            ps_pad,
            probs: self.probs.clone(),
            zs_pad,
            wzs_pad,
            wls_pad,
            boost_type: self.boost_type,
            loss_func: self.loss_func,
            sample_weight_clip: self.sample_weight_clip,
            sample_weight_wls_clip: self.sample_weight_wls_clip,
        }
    }

    // for debug
    pub fn _new_test(n: usize, ps_pad: Vec<f64>) -> Self {
        let len_n = ps_pad.len();
        SampleWeight {
            n,
            len_n,
            ps_pad: Some(ps_pad),
            ..Default::default()
        }
    }

    /*
    // use renew_sample_weight() with score=[0.0]
    pub fn renew_sample_weight_init(
        &mut self,
        phe: &Phe,
    ) {
    }
    */

    pub fn renew_sample_weight(
        &mut self,
        scores: &SampleScore,
        phe: &Phe,
        //boost_type: BoostType,
        //loss_func: LossFunc,
        //sample_weight_clip: SampleWeightClip,
        //sample_weight_wls_clip: SampleWeightWlsClip,
    ) {
        if self.boost_type.is_type_ada() {
            let loss_func = self.loss_func;
            let sample_weight_clip = self.sample_weight_clip;
            renew_ws(
                self.ws_mut().unwrap(),
                //self.ws_pad_mut().unwrap(),
                scores.scores(),
                //&scores,
                phe,
                loss_func,
                sample_weight_clip,
            );
            // TODO: how to avoid clone?
            // simultaneously calling ws() and ps_mut() is prohibited.
            // directly call self.ps_pad here
            // https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=e1ae68ad90ad083665a4fb157bcf3b8f
            //
            // error
            //let ws = self.ws().unwrap().to_vec();
            //let ws = self.ws().unwrap();
            //let ps = self.ps_mut().unwrap();

            // cannot call ps_mut() here
            let n = self.n();
            let ps = &mut self.ps_pad.as_deref_mut().unwrap()[..n];
            let ws = self.ws.as_deref().unwrap();

            renew_ps(ps, ws);
            //renew_ps(self.ps_mut().unwrap(), &ws);
        } else if self.boost_type.is_type_logit() {
            renew_logit_prob(self.probs_mut().unwrap(), &scores.scores());

            let n = self.n();
            let probs = self.probs.as_deref().unwrap();
            let zs = &mut self.zs_pad.as_deref_mut().unwrap()[..n];
            //let probs = self.probs().unwrap().to_vec();
            //let sample_weight_clip = self.sample_weight_clip;
            renew_logit_zs(zs, probs, phe, self.sample_weight_clip);
            //renew_logit_zs(self.zs_mut().unwrap(), &probs, phe, self.sample_weight_clip);

            let wls = &mut self.wls_pad.as_deref_mut().unwrap()[..n];
            //let sample_weight_wls_clip = self.sample_weight_wls_clip;
            renew_logit_wls(wls, probs, phe, self.sample_weight_wls_clip);

            let wzs = &mut self.wzs_pad.as_deref_mut().unwrap()[..n];
            //let wls = self.wls().unwrap().to_vec();
            //let zs = self.zs().unwrap().to_vec();
            renew_logit_wzs(wzs, wls, zs);
            //renew_logit_wzs(self.wzs_mut().unwrap(), &wls, &zs);
        } else {
            panic!();
        }

        //pub fn renew_sample_weight(
        //    &mut self,
        //    scores: &[f64],
        //    phe: &Phe,
        //    //boost_type: BoostType,
        //    //loss_func: LossFunc,
        //    //sample_weight_clip: SampleWeightClip,
        //    //sample_weight_wls_clip: SampleWeightWlsClip,
        //) {
        //    if self.boost_type.is_type_ada() {
        //        let loss_func = self.loss_func;
        //        let sample_weight_clip = self.sample_weight_clip;
        //        renew_ws(
        //            self.ws_mut().unwrap(),
        //            //self.ws_pad_mut().unwrap(),
        //            &scores,
        //            phe,
        //            loss_func,
        //            sample_weight_clip,
        //        );
        //        // TODO: how to avoid clone?
        //        let ws = self.ws().unwrap().to_vec();
        //        renew_ps(self.ps_mut().unwrap(), &ws);
        //    } else if self.boost_type.is_type_logit() {
        //        renew_logit_prob(self.probs_mut().unwrap(), &scores);

        //        let probs = self.probs().unwrap().to_vec();
        //        //let probs = self.probs_pad.clone().unwrap();
        //        let sample_weight_clip = self.sample_weight_clip;
        //        renew_logit_zs(self.zs_mut().unwrap(), &probs, phe, sample_weight_clip);

        //        let sample_weight_wls_clip = self.sample_weight_wls_clip;
        //        renew_logit_wls(self.wls_mut().unwrap(), &probs, phe, sample_weight_wls_clip);

        //        let wls = self.wls().unwrap().to_vec();
        //        let zs = self.zs().unwrap().to_vec();
        //        renew_logit_wzs(self.wzs_mut().unwrap(), &wls, &zs);
        //    } else {
        //        panic!();
        //    }
    }

    pub fn print_stat(&self, phe: &Phe) {
        if let Some(x) = self.ws() {
            self._print_stat_each(x, phe, "ws");
        }
        if let Some(x) = self.ps() {
            self._print_stat_each(x, phe, "ps");
        }
        if let Some(x) = self.probs() {
            self._print_stat_each(x, phe, "probs");
        }
        if let Some(x) = self.zs() {
            self._print_stat_each(x, phe, "zs");
        }
        if let Some(x) = self.wzs() {
            self._print_stat_each(x, phe, "wzs");
        }
        if let Some(x) = self.wls() {
            self._print_stat_each(x, phe, "wls");
        }
        //if let Some(x) = &self.ws_pad {
        //    self._print_stat_each(x, phe, "ws");
        //}
        //if let Some(x) = &self.ps_pad {
        //    self._print_stat_each(x, phe, "ps");
        //}
        //if let Some(x) = &self.probs_pad {
        //    self._print_stat_each(x, phe, "probs");
        //}
        //if let Some(x) = &self.zs_pad {
        //    self._print_stat_each(x, phe, "zs");
        //}
        //if let Some(x) = &self.wzs_pad {
        //    self._print_stat_each(x, phe, "wzs");
        //}
        //if let Some(x) = &self.wls_pad {
        //    self._print_stat_each(x, phe, "wls");
        //}
    }

    fn _print_stat_each(&self, ps: &[f64], phe: &Phe, name: &str) {
        log::debug!("{:?} stat case/control", name);
        log::debug!(
            "{:?} sum {:.4e}, {:.4e}",
            name,
            sum_ps(ps, phe, true),
            sum_ps(ps, phe, false)
        );
        log::debug!(
            "{:?} mean {:.4e}, {:.4e}",
            name,
            mean_ps(ps, phe, true),
            mean_ps(ps, phe, false)
        );
        log::debug!(
            "{:?} med {:.4e}, {:.4e}",
            name,
            med_ps(ps, phe, true),
            med_ps(ps, phe, false)
        );
        log::debug!(
            "{:?} max {:.4e}, {:.4e}",
            name,
            max_ps(ps, phe, true),
            max_ps(ps, phe, false)
        );
        log::debug!(
            "{:?} min {:.4e}, {:.4e}",
            name,
            min_ps(ps, phe, true),
            min_ps(ps, phe, false)
        );
    }
}

// use crate::score::add_score()
//// TODO: integrate to score.rs::add_score
//// TODO: implement constonly
//pub fn renew_score(score_ti: &mut [f64], pred: &[u8], wgt: &WgtBoost) {
//    match wgt.wgt().model().coef() {
//        Coef::Single(alpha) => {
//            for score in score_ti.iter_mut() {
//                *score += alpha;
//            }
//        }
//        Coef::Binary((const_ti, alpha_ti)) => {
//            let score_add_high = const_ti + alpha_ti;
//            let score_add_low = const_ti - alpha_ti;
//            for (ni, score) in score_ti.iter_mut().enumerate() {
//                let score_add = match pred[ni] {
//                    1 => score_add_high,
//                    0 => score_add_low,
//                    _ => panic!("Sth wrong"),
//                };
//
//                *score += score_add;
//            }
//        }
//        Coef::Score3((s0, s1, s2)) => {
//            for (ni, score) in score_ti.iter_mut().enumerate() {
//                let score_add = match pred[ni] {
//                    0 => s0,
//                    1 => s1,
//                    2 => s2,
//                    _ => panic!("Sth wrong"),
//                };
//
//                *score += score_add;
//            }
//        }
//        Coef::Score4((s0, s1, s2, m)) => {
//            for (ni, score) in score_ti.iter_mut().enumerate() {
//                let score_add = match pred[ni] {
//                    0 => s0,
//                    1 => s1,
//                    2 => s2,
//                    3 => m,
//                    _ => panic!("Sth wrong"),
//                };
//
//                *score += score_add;
//            }
//        }
//        Coef::Linear(_) => {
//            unimplemented!();
//            // -> these are Single()
//            /*
//            let score_add_high = alpha_ti;
//            let score_add_low = -alpha_ti;
//            for (ni, score) in score_ti.iter_mut().enumerate() {
//                let score_add = match pred[ni] {
//                    1 => score_add_high,
//                    0 => score_add_low,
//                    _ => panic!("Sth wrong"),
//                };
//                /*                 let score_add = match bit::bget(&pred, ni) {
//                    true => score_add_high,
//                    false => score_add_low,
//                }; */
//
//                *score += score_add;
//            }
//             */
//        }
//        Coef::LinearConst((const_ti, alpha_ti)) => {
//            let s2 = const_ti + 2.0 * alpha_ti;
//            let s1 = const_ti + alpha_ti;
//            let s0 = const_ti;
//            // TODO: use zip for pred
//            for (ni, score) in score_ti.iter_mut().enumerate() {
//                // TODO: *score += alpha_ti * val;
//                let score_add = match pred[ni] {
//                    0 => s0,
//                    1 => s1,
//                    2 => s2,
//                    //3 => m,
//                    _ => panic!("Sth wrong"),
//                };
//
//                *score += score_add;
//            }
//        }
//        _ => unimplemented!(),
//    }
//}

// moved to boosting::score
/// for linear model
//pub fn renew_score_var(score_ti: &mut [f64], coef: Coef, vals: &[f64]) {
//    //if let WgtKind::Snv(snv_wgt) = wgt.get_kind() {
//
//    //match wgt.wgt().model().coef() {
//    match coef {
//        Coef::Linear(alpha_ti) => {
//            //let vals = cov_wgt.var().vals();
//            for (score, val) in score_ti.iter_mut().zip(vals.iter()) {
//                *score += alpha_ti * val;
//                //log::debug!("score {}", score);
//            }
//        }
//        _ => unimplemented!(),
//    }
//}

pub fn renew_ps(ps: &mut [f64], ws: &[f64]) -> f64 {
    assert_eq!(ps.len(), ws.len());

    let ws_sum = ws.iter().sum();

    ps.iter_mut()
        .zip(ws.iter())
        .for_each(|(p, w)| *p = *w / ws_sum);

    ws_sum
}

pub fn extract_ws_label(ws: &[f64], phe: &Phe, is_case: bool) -> Vec<f64> {
    ws.iter()
        .zip(phe.iter())
        .filter(|(_, p)| *p == is_case)
        .map(|(w, _)| *w)
        .collect()
}

fn clip_ws_top(ws: &mut [f64], phe: &Phe, prop: f64) {
    for is_case in [true, false] {
        let mut ws_label = extract_ws_label(ws, phe, is_case);

        vec::sort_float(&mut ws_label);
        let n = ws_label.len();
        let top_n = ((n as f64) * (1.0 - prop)) as usize;
        let ws_lim = ws_label[top_n];
        log::debug!("cliptop ws_lim {:.6e}", ws_lim);
        ws.iter_mut()
            .zip(phe.iter())
            .filter(|(w, b)| (**w > ws_lim) && (*b == is_case))
            .for_each(|(w, _)| *w = ws_lim);
    }
}

fn clip_ws_both(ws: &mut [f64], phe: &Phe, prop: f64) {
    for is_case in [true, false] {
        let mut ws_label = extract_ws_label(ws, phe, is_case);

        vec::sort_float(&mut ws_label);
        let n = ws_label.len();
        let top_n = ((n as f64) * (1.0 - prop)) as usize;
        let ws_top_lim = ws_label[top_n];
        log::debug!("cliptop ws_top_lim {:.6e}", ws_top_lim);
        ws.iter_mut()
            .zip(phe.iter())
            .filter(|(w, b)| (**w > ws_top_lim) && (*b == is_case))
            .for_each(|(w, _)| *w = ws_top_lim);

        let bot_n = ((n as f64) * prop) as usize;
        let ws_bot_lim = ws_label[bot_n];
        log::debug!("clipbot ws_bot_lim {:.6e}", ws_bot_lim);
        ws.iter_mut()
            .zip(phe.iter())
            .filter(|(w, b)| (**w < ws_bot_lim) && (*b == is_case))
            .for_each(|(w, _)| *w = ws_bot_lim);
    }
}

fn renew_ws_f(ws: &mut [f64], score: &[f64], phe: &Phe, f: fn(f64, bool) -> f64) {
    ws.iter_mut()
        .zip(score.iter())
        .zip(phe.iter())
        .for_each(|((w, s), p)| {
            *w = f(*s, p);
        });
}

pub fn renew_ws(
    ws: &mut [f64],
    score: &[f64],
    phe: &Phe,
    loss_func: LossFunc,
    clip_weight: Option<SampleWeightClip>,
) {
    let renew_f = match loss_func {
        // AdaOri can be the same as Ada?
        LossFunc::Exp => renew_w_exp,
        LossFunc::Logistic => renew_w_logistic,
    };
    renew_ws_f(ws, score, phe, renew_f);

    if clip_weight.is_none() {
        return;
    }
    match clip_weight.unwrap() {
        SampleWeightClip::Top(prop) => clip_ws_top(ws, phe, prop),
        SampleWeightClip::Both(prop) => clip_ws_both(ws, phe, prop),
    }
}

//pub fn renew_ws(
//    ws: &mut [f64],
//    score: &[f64],
//    phe: &Phe,
//    loss_func: LossFunc,
//    clip_weight: Option<SampleWeightClip>,
//) {
//    let renew_f = match loss_func {
//        // AdaOri can be the same as Ada?
//        LossFunc::Exp => renew_w_exp,
//        LossFunc::Logistic => renew_w_logistic,
//    };
//    renew_ws_f(ws, score, phe, renew_f);
//
//    if clip_weight.is_none() {
//        return;
//    }
//    match clip_weight.unwrap() {
//        SampleWeightClip::Top(prop) => clip_ws_top(ws, phe, prop),
//        SampleWeightClip::Both(prop) => clip_ws_both(ws, phe, prop),
//    }
//}

fn renew_w_exp(s: f64, b: bool) -> f64 {
    let pn = (2 * b as i32 - 1) as f64;
    (-pn * s).exp()
}
fn renew_w_logistic(s: f64, b: bool) -> f64 {
    // pn: +1 or -1
    let pn = (2 * b as i32 - 1) as f64;
    1.0f64 / (1.0 + (pn * s).exp())
}

pub fn extract_ps_label(ps: &[f64], phe: &Phe, is_case: bool) -> Vec<f64> {
    //let n = phe.n();
    //ps[..n]
    ps.iter()
        .zip(phe.iter())
        .filter(|(_, ph)| *ph == is_case)
        .map(|(p, _)| *p)
        .collect()
}

/// Median of ps
/// Not cosidering when n is odd.
pub fn med_ps(ps: &[f64], phe: &Phe, is_case: bool) -> f64 {
    let mut ps_label: Vec<f64> = extract_ps_label(ps, phe, is_case);
    vec::sort_float(&mut ps_label);
    let med_n = ps_label.len() / 2;
    ps_label[med_n]
}

/// Mean of ps
pub fn mean_ps(ps: &[f64], phe: &Phe, is_case: bool) -> f64 {
    let ps_label: Vec<f64> = extract_ps_label(ps, phe, is_case);

    ps_label.iter().sum::<f64>() / (ps_label.len() as f64)
}

/// Min
pub fn min_ps(ps: &[f64], phe: &Phe, is_case: bool) -> f64 {
    let ps_label: Vec<f64> = extract_ps_label(ps, phe, is_case);

    ps_label.iter().fold(0.0 / 0.0, |v, v1| v1.min(v))
}

/// Max
pub fn max_ps(ps: &[f64], phe: &Phe, is_case: bool) -> f64 {
    let ps_label: Vec<f64> = extract_ps_label(ps, phe, is_case);

    ps_label.iter().fold(0.0 / 0.0, |v, v1| v1.max(v))
}

/// Sum
pub fn sum_ps(ps: &[f64], phe: &Phe, is_case: bool) -> f64 {
    let ps_label: Vec<f64> = extract_ps_label(ps, phe, is_case);

    ps_label.iter().sum::<f64>()
}

fn renew_prob_logit(s: f64) -> f64 {
    1.0f64 / (1.0 + (-s).exp())
}

pub fn renew_logit_prob(probs: &mut [f64], scores: &[f64]) {
    assert_eq!(probs.len(), scores.len());

    probs
        .iter_mut()
        .zip(scores.iter())
        .for_each(|(p, s)| *p = renew_prob_logit(*s));
    /*     let n = probs.len();
    probs[..n]
        .iter_mut()
        .zip(scores.iter())
        .for_each(|(p, s)| *p = renew_prob_logit(*s)); */
}

fn renew_z_logit(pr: f64, b: bool) -> f64 {
    if b {
        1.0 / pr
    } else {
        -1.0 / (1.0 - pr)
    }
}

pub fn renew_logit_zs(
    zs: &mut [f64],
    probs: &[f64],
    phe: &Phe,
    //loss_func: LossFunc,
    clip_weight: Option<SampleWeightClip>,
) {
    assert_eq!(zs.len(), probs.len());
    assert_eq!(zs.len(), phe.n());

    zs.iter_mut()
        .zip(probs.iter())
        .zip(phe.iter())
        .for_each(|((z, pr), p)| {
            *z = renew_z_logit(*pr, p);
        });

    //let renew_f = match loss_func {
    //    // AdaOri can be the same as Ada?
    //    LossFunc::Exp => renew_w_exp,
    //    LossFunc::Logistic => renew_w_logistic,
    //    _ => panic!("wrong")
    //};
    //renew_ws_f(ws, score, phe, renew_f);

    if clip_weight.is_none() {
        return;
    }
    match clip_weight.unwrap() {
        SampleWeightClip::Top(prop) => clip_logit_zs_top(zs, phe, prop),
        SampleWeightClip::Both(prop) => clip_logit_zs_both(zs, phe, prop),
    }
}

fn clip_logit_zs_top(ws: &mut [f64], phe: &Phe, prop: f64) {
    // clip zs
    // for phe=+1, clip top zs
    // for phe=-1, clip bottom zs (since zs is all negative)
    for is_case in [true, false] {
        let mut ws_label = extract_ws_label(ws, phe, is_case);

        vec::sort_float(&mut ws_label);
        let n = ws_label.len();
        let top_n = ((n as f64) * (1.0 - prop)) as usize;
        let ws_lim = ws_label[top_n];
        log::debug!("cliptop ws_lim {:.6e}", ws_lim);

        if is_case {
            ws.iter_mut()
                .zip(phe.iter())
                .filter(|(w, b)| (**w > ws_lim) && (*b == is_case))
                .for_each(|(w, _)| *w = ws_lim);
        } else {
            ws.iter_mut()
                .zip(phe.iter())
                .filter(|(w, b)| (**w < ws_lim) && (*b == is_case))
                .for_each(|(w, _)| *w = ws_lim);
        }
    }
}

// clipboth not suitable for logit
fn clip_logit_zs_both(zs: &mut [f64], phe: &Phe, prop: f64) {
    clip_ws_both(zs, phe, prop)
}

fn renew_wls_logit(pr: f64) -> f64 {
    pr * (1.0 - pr)
}

pub fn renew_logit_wls(
    wls: &mut [f64],
    probs: &[f64],
    phe: &Phe,
    clip_weight: Option<SampleWeightWlsClip>,
) {
    assert_eq!(wls.len(), probs.len());
    assert_eq!(wls.len(), phe.n());

    wls.iter_mut()
        .zip(probs.iter())
        .for_each(|(w, pr)| *w = renew_wls_logit(*pr));

    /*     let n = probs.len();
    wls[..n]
        .iter_mut()
        .zip(probs.iter())
        .for_each(|(w, pr)| *w = renew_wls_logit(*pr)); */

    if clip_weight.is_none() {
        return;
    }
    match clip_weight.unwrap() {
        // same as ws in modelfree
        SampleWeightWlsClip::Top(prop) => clip_ws_top(wls, phe, prop),
        SampleWeightWlsClip::Both(prop) => clip_ws_both(wls, phe, prop),
    }
}

fn renew_wz_logit(wl: f64, z: f64) -> f64 {
    wl * z
}

pub fn renew_logit_wzs(wzs: &mut [f64], wls: &[f64], zs: &[f64]) {
    assert_eq!(wzs.len(), wls.len());
    assert_eq!(wzs.len(), zs.len());

    wzs.iter_mut()
        .zip(wls.iter())
        .zip(zs.iter())
        .for_each(|((wz, wl), z)| {
            *wz = renew_wz_logit(*wl, *z);
        });

    /*     let n = zs.len();
    wzs[..n]
        .iter_mut()
        .zip(wls.iter())
        .zip(zs.iter())
        .for_each(|((wz, wl), z)| {
            *wz = renew_wz_logit(*wl, *z);
        }); */
}

#[cfg(test)]
mod tests {
    use super::*;
    use genetics::alloc;

    fn setup() -> SampleWeight {
        let ps_pad = Some(alloc::vec_align_f64(vec![0.1f64, 0.2, 0.3], 5));
        //let ps_pad = Some(vec![0.1f64, 0.2, 0.3, 0.0, 0.0]);
        SampleWeight {
            n: 3,
            len_n: 5,
            ps_pad,
            ..Default::default()
        }
    }

    #[test]
    fn test_sample_weight_len() {
        let mut sw = setup();

        assert_eq!(sw.ps_pad().unwrap().len(), 5);
        assert_eq!(sw.ps_pad_mut().unwrap().len(), 5);
        assert_eq!(sw.ps().unwrap().len(), 3);
        assert_eq!(sw.ps_mut().unwrap().len(), 3);
    }

    /// test padding will not changed
    #[test]
    fn test_sample_weight_mut() {
        let mut sw = setup();

        let ps_mut = sw.ps_mut().unwrap();
        for p in ps_mut.iter_mut() {
            *p = *p + 1.0
        }

        assert_eq!(sw.ps_pad().unwrap(), vec![1.1f64, 1.2, 1.3, 0.0, 0.0]);
    }

    #[test]
    fn test_sample_weight_clear_pad() {
        let mut sw = setup();
        sw.ps_pad_mut().unwrap()[4] = 1.0;
        assert_eq!(sw.ps_pad_mut().unwrap()[4], 1.0);
        sw.clear_pad();
        assert_eq!(sw.ps_pad().unwrap()[4], 0.0);
    }

    #[test]
    fn test_renew_ps() {
        let ws = vec![1.0, 2.0, 3.0, 4.0];
        let mut ps = vec![0.0; 4];
        //let mut ps = vec![0.0; 6];
        renew_ps(&mut ps, &ws);
        assert_eq!(ps, vec![0.1, 0.2, 0.3, 0.4]);
        //assert_eq!(ps, vec![0.1, 0.2, 0.3, 0.4, 0.0, 0.0]);
    }

    #[test]
    fn test_clip_ws_top() {
        let mut ws = vec![1.0, 2.0, 3.0, 10.0, 1.0];
        let phe_v = vec![true, true, true, true, false];
        let phe = Phe::new(&phe_v);
        clip_ws_top(&mut ws, &phe, 0.3);
        assert_eq!(ws, vec![1.0, 2.0, 3.0, 3.0, 1.0]);
    }

    #[test]
    fn test_clip_ws_top_control() {
        let mut ws = vec![1.0, 2.0, 3.0, 10.0, 1.0];
        let phe_v = vec![false, false, false, false, true];
        let phe = Phe::new(&phe_v);
        clip_ws_top(&mut ws, &phe, 0.3);
        assert_eq!(ws, vec![1.0, 2.0, 3.0, 3.0, 1.0]);
    }

    #[test]
    fn test_extract_ps_label() {
        let mut ps = vec![0.1, 0.2, 0.3, 0.4, 1.0, 1.0, 1.0];
        let phe_v = vec![true, true, false, false];
        let phe = Phe::new(&phe_v);
        let ps_case = extract_ps_label(&mut ps, &phe, true);
        assert_eq!(ps_case, vec![0.1, 0.2]);
        let ps_cont = extract_ps_label(&mut ps, &phe, false);
        assert_eq!(ps_cont, vec![0.3, 0.4]);
    }
}
