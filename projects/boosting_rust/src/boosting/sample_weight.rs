use super::LossFunc;
use super::WgtBoost;
use crate::boosting_param::SampleWeightClip;
use genetics::samples::prelude::*;
use genetics::vec;
use genetics::wgt::Coef;

// TODO: integrate to score.rs::add_score
// TODO: implement constonly
pub fn renew_score(score_ti: &mut [f64], pred: &[u8], wgt: &WgtBoost) {
    match wgt.wgt().model().coef() {
        Coef::Binary((const_ti, alpha_ti)) => {
            let score_add_high = const_ti + alpha_ti;
            let score_add_low = const_ti - alpha_ti;
            for (ni, score) in score_ti.iter_mut().enumerate() {
                let score_add = match pred[ni] {
                    1 => score_add_high,
                    0 => score_add_low,
                    _ => panic!("Sth wrong"),
                };

                *score += score_add;
            }
        }
        Coef::Score4((s0, s1, s2, m)) => {
            for (ni, score) in score_ti.iter_mut().enumerate() {
                let score_add = match pred[ni] {
                    0 => s0,
                    1 => s1,
                    2 => s2,
                    3 => m,
                    _ => panic!("Sth wrong"),
                };

                *score += score_add;
            }
        }
        Coef::Linear(alpha_ti) => {
            let score_add_high = alpha_ti;
            let score_add_low = -alpha_ti;
            for (ni, score) in score_ti.iter_mut().enumerate() {
                let score_add = match pred[ni] {
                    1 => score_add_high,
                    0 => score_add_low,
                    _ => panic!("Sth wrong"),
                };
                /*                 let score_add = match bit::bget(&pred, ni) {
                    true => score_add_high,
                    false => score_add_low,
                }; */

                *score += score_add;
            }
        }
        _ => unimplemented!(),
    }
}

// TODO: combine into above?
/// for linearmodel
pub fn renew_score_var(score_ti: &mut [f64], coef: Coef, vals: &[f64]) {
    //if let WgtKind::Snv(snv_wgt) = wgt.get_kind() {

    //match wgt.wgt().model().coef() {
    match coef {
        Coef::Linear(alpha_ti) => {
            //let vals = cov_wgt.var().vals();
            for (score, val) in score_ti.iter_mut().zip(vals.iter()) {
                *score += alpha_ti * val;
                //println!("score {}", score);
            }
        }
        _ => unimplemented!(),
    }
}

pub fn renew_ps(ps: &mut [f64], ws: &[f64]) -> f64 {
    let ws_sum = ws.iter().sum();
    let n = ws.len();

    ps[..n]
        .iter_mut()
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
        println!("cliptop ws_lim {}", ws_lim);
        ws.iter_mut()
            .zip(phe.iter())
            .filter(|(w, b)| (**w > ws_lim) & (*b == is_case))
            .for_each(|(w, _)| *w = ws_lim);
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
        _ => unimplemented!(),
    }
}

fn renew_w_exp(s: f64, b: bool) -> f64 {
    let pn = (2 * b as i32 - 1) as f64;
    (-pn * s).exp()
}
fn renew_w_logistic(s: f64, b: bool) -> f64 {
    // phe: +1 or -1
    let pn = (2 * b as i32 - 1) as f64;
    1.0f64 / (1.0 + (pn * s).exp())
}

pub fn extract_ps_label(ps: &[f64], phe: &Phe, is_case: bool) -> Vec<f64> {
    let n = phe.n();
    ps[..n]
        .iter()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_renew_ps() {
        let ws = vec![1.0, 2.0, 3.0, 4.0];
        let mut ps = vec![0.0; 6];
        renew_ps(&mut ps, &ws);
        assert_eq!(ps, vec![0.1, 0.2, 0.3, 0.4, 0.0, 0.0]);
    }

    #[test]
    fn test_clip_ws_top() {
        let mut ws = vec![1.0, 2.0, 3.0, 10.0, 1.0];
        let phe_v = vec![true, true, true, true, false];
        let phe = Phe::new(phe_v);
        clip_ws_top(&mut ws, &phe, 0.3);
        assert_eq!(ws, vec![1.0, 2.0, 3.0, 3.0, 1.0]);
    }

    #[test]
    fn test_clip_ws_top_control() {
        let mut ws = vec![1.0, 2.0, 3.0, 10.0, 1.0];
        let phe_v = vec![false, false, false, false, true];
        let phe = Phe::new(phe_v);
        clip_ws_top(&mut ws, &phe, 0.3);
        assert_eq!(ws, vec![1.0, 2.0, 3.0, 3.0, 1.0]);
    }

    #[test]
    fn test_extract_ps_label() {
        let mut ps = vec![0.1, 0.2, 0.3, 0.4, 1.0, 1.0, 1.0];
        let phe_v = vec![true, true, false, false];
        let phe = Phe::new(phe_v);
        let ps_case = extract_ps_label(&mut ps, &phe, true);
        assert_eq!(ps_case, vec![0.1, 0.2]);
        let ps_cont = extract_ps_label(&mut ps, &phe, false);
        assert_eq!(ps_cont, vec![0.3, 0.4]);
    }
}
