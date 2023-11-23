use crate::EffEps;
//use genetics::wgt::Coef;

///
///  scores: // (s0, s1, s2)
///  return scores: // (s0, s1, s2)
///
pub fn adjust_eff_logit_no_missing(
    scores: (f64, f64, f64), // (s0, s1, s2)
    table8_count: (usize, usize, usize, usize, usize, usize, usize, usize),
    eff_eps: Option<EffEps>,
    verbose: bool,
) -> ((f64, f64, f64), bool) {
    // since will be used to Modelfree()
    //) -> (Coef, bool) {
    let (s0, s1, s2) = scores;

    if let Some(eff_eps) = eff_eps {
        match eff_eps {
            EffEps::LimScore(clim, sratio) => {
                let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
                let (s2_new, is_eff_eps_2) =
                    adjust_eff_score_lims2(s2, d2, n2, (s1, s0), clim, sratio);
                let (s1_new, is_eff_eps_1) =
                    adjust_eff_score_lims2(s1, d1, n1, (s2, s0), clim, sratio);
                let (s0_new, is_eff_eps_0) =
                    adjust_eff_score_lims2(s0, d0, n0, (s2, s1), clim, sratio);
                let is_eff_eps = is_eff_eps_2 | is_eff_eps_1 | is_eff_eps_0;
                ((s0_new, s1_new, s2_new), is_eff_eps)
                //(Coef::Score3((s0_new, s1_new, s2_new)), is_eff_eps)
            }
            EffEps::LimScoreProp(plim, sratio) => {
                //let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
                let (s2_new, is_eff_eps_2) = adjust_eff_score_lim_score_prop(
                    s2,
                    table8_prop(table8_count, 2),
                    (s1, s0),
                    plim,
                    sratio,
                );
                let (s1_new, is_eff_eps_1) = adjust_eff_score_lim_score_prop(
                    s1,
                    table8_prop(table8_count, 1),
                    (s2, s0),
                    plim,
                    sratio,
                );
                let (s0_new, is_eff_eps_0) = adjust_eff_score_lim_score_prop(
                    s0,
                    table8_prop(table8_count, 0),
                    (s2, s1),
                    plim,
                    sratio,
                );
                let is_eff_eps = is_eff_eps_2 | is_eff_eps_1 | is_eff_eps_0;
                ((s0_new, s1_new, s2_new), is_eff_eps)
                //(Coef::Score3((s0_new, s1_new, s2_new)), is_eff_eps)
            }
            EffEps::LimS2GmodelProp(plim, rec_max_ratio, rec_min_ratio) => {
                //let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
                let (s2_new, is_eff_eps_2) = adjust_eff_score_lim_s2_gmodel_prop(
                    s2,
                    table8_prop(table8_count, 2),
                    s0,
                    s1,
                    plim,
                    rec_max_ratio,
                    rec_min_ratio,
                );
                let is_eff_eps = is_eff_eps_2;
                ((s0, s1, s2_new), is_eff_eps)
                //(Coef::Score3((s0, s1, s2_new)), is_eff_eps)
            }
            EffEps::LimS12GmodelProp(
                plim,
                rec_max_ratio,
                rec_min_ratio,
                het_max_ratio,
                het_min_ratio,
            )
            | EffEps::LimS12GmodelPropOnUpdate(
                plim,
                rec_max_ratio,
                rec_min_ratio,
                het_max_ratio,
                het_min_ratio,
            ) => {
                //let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
                let (s2_new, s1_new, is_eff_eps_2) = adjust_eff_score_lim_s12_gmodel_prop(
                    s2,
                    s1,
                    table8_prop(table8_count, 2),
                    s0,
                    plim,
                    rec_max_ratio,
                    rec_min_ratio,
                    het_max_ratio,
                    het_min_ratio,
                );
                let is_eff_eps = is_eff_eps_2;
                ((s0, s1_new, s2_new), is_eff_eps)
                //(Coef::Score3((s0, s1_new, s2_new)), is_eff_eps)
            }
            EffEps::LimS2GmodelOverProp(
                plim,
                rec_max_ratio,
                rec_min_ratio,
                het_max_ratio,
                //het_min_ratio,
            )
            | EffEps::LimS2GmodelOverPropOnUpdate(
                plim,
                rec_max_ratio,
                rec_min_ratio,
                het_max_ratio,
            ) => {
                //let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
                let (s2_new, is_eff_eps_2) = adjust_eff_score_lim_s2_gmodel_over_prop(
                    s2,
                    s1,
                    table8_prop(table8_count, 2),
                    s0,
                    plim,
                    rec_max_ratio,
                    rec_min_ratio,
                    het_max_ratio,
                    //het_min_ratio,
                );
                let is_eff_eps = is_eff_eps_2;
                ((s0, s1, s2_new), is_eff_eps)
                //(Coef::Score3((s0, s1, s2_new)), is_eff_eps)
            }
            EffEps::LimS2GmodelOverKeepSignProp(
                plim,
                rec_max_ratio,
                rec_min_ratio,
                het_max_ratio,
            ) => {
                //let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
                let (s2_new, is_eff_eps_2) = adjust_eff_score_lim_s2_gmodel_over_keep_sign_prop(
                    s2,
                    s1,
                    table8_prop(table8_count, 2),
                    s0,
                    plim,
                    rec_max_ratio,
                    rec_min_ratio,
                    het_max_ratio,
                    //het_min_ratio,
                );
                let is_eff_eps = is_eff_eps_2;
                ((s0, s1, s2_new), is_eff_eps)
                //(Coef::Score3((s0, s1, s2_new)), is_eff_eps)
            }

            EffEps::LimS2GmodelBorderProp(
                plim,
                rec_add_ratio,
                add_dom_ratio,
                dom_het_ratio,
                het_rec_ratio,
            ) => {
                //let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
                let (s2_new, is_eff_eps_2) = adjust_eff_score_lim_s2_gmodel_border_prop(
                    s2,
                    table8_prop(table8_count, 2),
                    s0,
                    s1,
                    plim,
                    rec_add_ratio,
                    add_dom_ratio,
                    dom_het_ratio,
                    het_rec_ratio,
                    verbose,
                );
                let is_eff_eps = is_eff_eps_2;
                ((s0, s1, s2_new), is_eff_eps)
                //(Coef::Score3((s0, s1, s2_new)), is_eff_eps)
            }
            EffEps::LimS2AddProp(plim) => {
                //let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
                let (s2_new, is_eff_eps_2) = adjust_eff_score_lim_s2_add_prop(
                    s2,
                    table8_prop(table8_count, 2),
                    s0,
                    s1,
                    plim,
                );
                let is_eff_eps = is_eff_eps_2;
                ((s0, s1, s2_new), is_eff_eps)
                //(Coef::Score3((s0, s1, s2_new)), is_eff_eps)
            }
        }
    } else {
        ((s0, s1, s2), false)
        //(Coef::Score3((s0, s1, s2)), false)
    }
}

fn table8_prop(
    table8_count: (usize, usize, usize, usize, usize, usize, usize, usize),
    genot: usize,
) -> f64 {
    let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;
    let sum = (d2 + n2 + d1 + n1 + d0 + n0) as f64;
    match genot {
        2 => ((d2 + n2) as f64) / sum,
        1 => ((d1 + n1) as f64) / sum,
        0 => ((d0 + n0) as f64) / sum,
        _ => panic!(),
    }
}

/*
// TODO: test
fn table8_maf(table8_count: (usize, usize, usize, usize, usize, usize, usize, usize)) -> f64 {
    let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;

    let maf = (2.0 * (d2 + n2) + (d1 + n1)) / (2.0 * (d2 + n2 + d1 + n1 + d0 + n0));
    if maf > 0.5 {
        1.0 - maf
    } else {
        maf
    }
}
 */

fn lim_score(s: f64, s_others: (f64, f64), sratio: f64) -> (f64, bool) {
    // positive
    let s_others_diff = s_others.0.abs().max(s_others.1.abs());

    if s.abs() > sratio * s_others_diff {
        (s.signum() * (sratio * s_others_diff), true)
    } else {
        (s, false)
    }
}

fn lim_s2_add(_s2: f64, s1: f64, s0: f64) -> (f64, bool) {
    let diff_s1_s0 = s1 - s0;

    (s0 + 2.0 * diff_s1_s0, true)
}

fn adjust_eff_score_lims2(
    s: f64,
    d: usize,
    n: usize,
    s_others: (f64, f64),
    //eff_eps: EffEps,
    clim: usize,
    sratio: f64,
) -> (f64, bool) {
    if (d < clim) | (n < clim) {
        lim_score(s, s_others, sratio)
    } else {
        (s, false)
    }
}

fn adjust_eff_score_lim_score_prop(
    s: f64,
    p: f64,
    s_others: (f64, f64),
    //eff_eps: EffEps,
    plim: f64,
    sratio: f64,
) -> (f64, bool) {
    if p < plim {
        lim_score(s, s_others, sratio)
    } else {
        (s, false)
    }
}

fn adjust_eff_score_lim_s2_add_prop(
    s2: f64,
    p: f64,
    s0: f64,
    s1: f64,
    //s_others: (f64, f64),
    //eff_eps: EffEps,
    plim: f64,
    //sratio: f64,
) -> (f64, bool) {
    if p < plim {
        lim_s2_add(s2, s1, s0)
    } else {
        (s2, false)
    }
}

fn lim_s2_gmodel(s2: f64, s1: f64, s0: f64, rec_max_ratio: f64, rec_min_ratio: f64) -> (f64, bool) {
    let rs1 = s1 - s0;
    let rs2 = s2 - s0;

    let ratio = rs2.abs() / rs1.abs();

    if (rs2 * rs1).is_sign_positive() {
        // sign are same
        if ratio > rec_max_ratio {
            return (s0 + rec_max_ratio * rs1, true);
        }
    } else {
        // sign are different
        if ratio > rec_min_ratio {
            return (s0 + -rec_min_ratio * rs1, true);
        }
    }

    (s2, false)
}

fn lim_s12_gmodel(
    s2: f64,
    s1: f64,
    s0: f64,
    rec_max_ratio: f64,
    rec_min_ratio: f64,
    het_max_ratio: f64,
    het_min_ratio: f64,
) -> (f64, f64, bool) {
    let rs1 = s1 - s0;
    let rs2 = s2 - s0;

    let ratio2 = rs2.abs() / rs1.abs();

    // for s2
    if (rs2 * rs1).is_sign_positive() {
        // sign are same
        // TODO: this should be upper
        if ratio2 > rec_max_ratio {
            return (s0 + rec_max_ratio * rs1, s1, true);
        }
    } else {
        // sign are different
        if ratio2 > rec_min_ratio {
            return (s0 + -rec_min_ratio * rs1, s1, true);
        }
    }

    let ratio1 = rs1.abs() / rs2.abs();

    // for s1
    if (rs2 * rs1).is_sign_positive() {
        // sign are same
        if ratio1 > het_max_ratio {
            return (s2, s0 + het_max_ratio * rs2, true);
        }
    } else {
        // sign are different
        if ratio1 > het_min_ratio {
            return (s2, s0 + -het_min_ratio * rs2, true);
        }
    }

    (s2, s1, false)
}

// TOOD: when s2 is not minor homo but major homo
fn lim_s2_gmodel_over(
    s2: f64,
    s1: f64,
    s0: f64,
    rec_max_ratio: f64,
    rec_min_ratio: f64,
    het_max_ratio: f64,
) -> (f64, bool) {
    let rs1 = s1 - s0;
    let rs2 = s2 - s0;

    let ratio2 = rs2.abs() / rs1.abs();

    // for rec
    if (rs2 * rs1).is_sign_positive() {
        // sign are same
        if ratio2 > rec_max_ratio {
            return (s0 + rec_max_ratio * rs1, true);
        }
    } else {
        // sign are different
        if ratio2 > rec_min_ratio {
            return (s0 + -rec_min_ratio * rs1, true);
        }
    }

    let sratio2 = rs2 / rs1;

    // for overdom
    if sratio2 < (1.0 / het_max_ratio) {
        return (s0 + (1.0 / het_max_ratio) * rs1, true);
    }

    (s2, false)
}

fn lim_s2_gmodel_over_keep_sign(
    s2: f64,
    s1: f64,
    s0: f64,
    rec_max_ratio: f64,
    rec_min_ratio: f64,
    het_max_ratio: f64,
) -> (f64, bool) {
    let rs1 = s1 - s0;
    let rs2 = s2 - s0;

    let ratio2 = rs2.abs() / rs1.abs();

    // for rec
    if (rs2 * rs1).is_sign_positive() {
        // sign are same
        if ratio2 > rec_max_ratio {
            return (s0 + rec_max_ratio * rs1, true);
        }
    } else {
        // sign are different
        if ratio2 > rec_min_ratio {
            return (s0 + -rec_min_ratio * rs1, true);
        }
    }

    // for overdom
    if (rs2 * rs1).is_sign_negative() {
        return (s0, true);
    }

    // different het_max_ratio from gmodel_over()
    if ratio2 < het_max_ratio {
        return (s0 + het_max_ratio * rs1, true);
    }

    (s2, false)
}

fn adjust_eff_score_lim_s2_gmodel_prop(
    s2: f64,
    p: f64,
    s0: f64,
    s1: f64,
    //s_others: (f64, f64),
    //eff_eps: EffEps,
    plim: f64,
    rec_max_ratio: f64,
    rec_min_ratio: f64,
) -> (f64, bool) {
    if p < plim {
        lim_s2_gmodel(s2, s1, s0, rec_max_ratio, rec_min_ratio)
    } else {
        (s2, false)
    }
}

fn adjust_eff_score_lim_s12_gmodel_prop(
    s2: f64,
    s1: f64,
    p: f64,
    s0: f64,
    //s_others: (f64, f64),
    //eff_eps: EffEps,
    plim: f64,
    rec_max_ratio: f64,
    rec_min_ratio: f64,
    het_max_ratio: f64,
    het_min_ratio: f64,
) -> (f64, f64, bool) {
    if p < plim {
        lim_s12_gmodel(
            s2,
            s1,
            s0,
            rec_max_ratio,
            rec_min_ratio,
            het_max_ratio,
            het_min_ratio,
        )
    } else {
        (s2, s1, false)
    }
}

fn adjust_eff_score_lim_s2_gmodel_over_prop(
    s2: f64,
    s1: f64,
    p: f64,
    s0: f64,
    //s_others: (f64, f64),
    //eff_eps: EffEps,
    plim: f64,
    rec_max_ratio: f64,
    rec_min_ratio: f64,
    het_max_ratio: f64,
    //het_min_ratio: f64,
) -> (f64, bool) {
    if p < plim {
        lim_s2_gmodel_over(
            s2,
            s1,
            s0,
            rec_max_ratio,
            rec_min_ratio,
            het_max_ratio,
            //het_min_ratio,
        )
    } else {
        (s2, false)
    }
}

fn adjust_eff_score_lim_s2_gmodel_over_keep_sign_prop(
    s2: f64,
    s1: f64,
    p: f64,
    s0: f64,
    //s_others: (f64, f64),
    //eff_eps: EffEps,
    plim: f64,
    rec_max_ratio: f64,
    rec_min_ratio: f64,
    het_max_ratio: f64,
    //het_min_ratio: f64,
) -> (f64, bool) {
    if p < plim {
        lim_s2_gmodel_over_keep_sign(
            s2,
            s1,
            s0,
            rec_max_ratio,
            rec_min_ratio,
            het_max_ratio,
            //het_min_ratio,
        )
    } else {
        (s2, false)
    }
}

fn lim_s2_gmodel_border(
    s2: f64,
    s1: f64,
    s0: f64,
    rec_add_ratio: f64,
    add_dom_ratio: f64,
    dom_het_ratio: f64,
    het_rec_ratio: f64,
    verbose: bool,
) -> (f64, bool) {
    let diff_s1_s0 = s1 - s0;
    let diff_s2_s0 = s2 - s0;

    let ratio = diff_s1_s0 / diff_s2_s0;

    if verbose {
        log::debug!("gmodel_border rs1/rs2 {}", ratio);
    }

    if 0.0 <= ratio && ratio < rec_add_ratio {
        // ex. rec_add_ratio=0.25
        if verbose {
            log::debug!("gmodel_border: rec");
        }
        (s0 + diff_s1_s0 / rec_add_ratio, true)
    } else if rec_add_ratio <= ratio && ratio <= add_dom_ratio {
        if verbose {
            log::debug!("gmodel_border: add");
        }
        (s0 + diff_s1_s0 / 0.5, true)
    } else if add_dom_ratio < ratio && ratio <= dom_het_ratio {
        if verbose {
            log::debug!("gmodel_border: dom");
        }
        (s1, true)
    } else if dom_het_ratio < ratio || ratio < het_rec_ratio {
        if verbose {
            log::debug!("gmodel_border: het");
        }
        (s0, true)
    } else if het_rec_ratio <= ratio && ratio < 0.0 {
        if verbose {
            log::debug!("gmodel_border: rec");
        }
        (s0 + diff_s1_s0 / het_rec_ratio, true)
    } else {
        panic!("wrong");
    }

    /*
    // ex. (rratio, hratio)=(5.0,0.0), (5.0, -1.0)
    if ratio > rec_ratio {
        (s0 + rec_ratio * diff_s1_s0, true)
    } else if ratio < het_ratio {
        (s0 + het_ratio * diff_s1_s0, true)
    } else {
        (s2, false)
    }
    */

    //(s0 + 2.0 * diff_s1_s0, true)
}

fn adjust_eff_score_lim_s2_gmodel_border_prop(
    s2: f64,
    p: f64,
    s0: f64,
    s1: f64,
    //s_others: (f64, f64),
    //eff_eps: EffEps,
    plim: f64,
    rec_add_ratio: f64,
    add_dom_ratio: f64,
    dom_het_ratio: f64,
    het_rec_ratio: f64,
    verbose: bool,
) -> (f64, bool) {
    if p < plim {
        lim_s2_gmodel_border(
            s2,
            s1,
            s0,
            rec_add_ratio,
            add_dom_ratio,
            dom_het_ratio,
            het_rec_ratio,
            verbose,
        )
    } else {
        (s2, false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lim_score() {
        // no change
        let s = 0.1;
        let s_others = (0.1, 0.1);
        let sratio = 5.0;

        assert_float_absolute_eq!(s, lim_score(s, s_others, sratio).0);
        //assert!(is_eq_f64(s, lim_score(s, s_others, sratio).0, 1e-7));
    }

    #[test]
    fn test_lim_score_2() {
        // no change
        let s = -10.0;
        let s_others = (0.1, 0.1);
        let sratio = 5.0;

        assert_float_absolute_eq!(-0.5, lim_score(s, s_others, sratio).0);
    }

    #[test]
    fn test_lim_s12_gmodel() {
        // all patterns

        let (rec_max_ratio, rec_min_ratio) = (4.0, 4.0);
        let (het_max_ratio, het_min_ratio) = (1.25, 0.25);

        let s0 = 0.1;
        // pairs of (s1, s2) of original and answer
        let pairs = [
            ((0.3, 0.5), (0.3, 0.5, false)),     // no change
            ((0.3, 1.0), (0.3, 0.9, true)),      // s2 changes; (rs1, rs2)=(+,+)
            ((-0.1, 1.0), (-0.1, 0.9, true)),    // (-,+)
            ((0.3, -0.8), (0.3, -0.7, true)),    // (+,-)
            ((-0.1, -0.8), (-0.1, -0.7, true)),  // (-,-)
            ((0.4, 0.3), (0.35, 0.3, true)),     // s1 changes; (rs1, rs2)=(+,+)
            ((0.0, 0.3), (0.05, 0.3, true)),     // (-,+)
            ((0.2, -0.1), (0.15, -0.1, true)),   // (+,-)
            ((-0.2, -0.1), (-0.15, -0.1, true)), // (-,-)
        ];

        for ((s1, s2), (s1_ans, s2_ans, is_changed_ans)) in pairs.iter() {
            let (s2_ret, s1_ret, is_changed) = lim_s12_gmodel(
                *s2,
                *s1,
                s0,
                rec_max_ratio,
                rec_min_ratio,
                het_max_ratio,
                het_min_ratio,
            );
            println!("{}, {}, {}", s2, *s2_ans, s2_ret);

            assert_float_absolute_eq!(s2_ret, *s2_ans);
            assert_float_absolute_eq!(s1_ret, *s1_ans);
            assert_eq!(is_changed, *is_changed_ans);
        }
    }

    #[test]
    fn test_lim_s2_gmodel_over() {
        // all patterns

        let (rec_max_ratio, rec_min_ratio) = (4.0, 4.0);
        let het_max_ratio = 1.25;

        let s0 = 0.1;
        // pairs of (s1, s2) of original and answer
        let pairs = [
            ((0.3, 0.5), (0.5, false)),    // no change
            ((0.3, 1.0), (0.9, true)),     // s2 changes; (rs1, rs2)=(+,+)
            ((-0.1, 1.0), (0.9, true)),    // (-,+)
            ((0.3, -0.8), (-0.7, true)),   // (+,-)
            ((-0.1, -0.8), (-0.7, true)),  // (-,-)
            ((0.4, 0.3), (0.34, true)),    // s1 changes; (rs1, rs2)=(+,+)
            ((0.0, 0.3), (0.02, true)),    // (-,+)
            ((0.2, -0.1), (0.18, true)),   // (+,-)
            ((-0.2, -0.1), (-0.14, true)), // (-,-)
        ];

        for ((s1, s2), (s2_ans, is_changed_ans)) in pairs.iter() {
            let (s2_ret, is_changed) =
                lim_s2_gmodel_over(*s2, *s1, s0, rec_max_ratio, rec_min_ratio, het_max_ratio);
            println!("{}, {}, {}", s2, *s2_ans, s2_ret);

            assert_float_absolute_eq!(s2_ret, *s2_ans);
            //assert!(is_eq_f64(s1_ret, *s1_ans, 1e-7));
            assert_eq!(is_changed, *is_changed_ans);
        }
    }
}
