use super::epsilon;
use crate::{BoostType, Eps};
use genetics::genot::prelude::*;
use genetics::samples::prelude::*;

pub const CONTINGENCY_TABLE_FILL: f64 = 1e-12;

pub const MACHINE_EPS: f64 = 1e-12;

//pub type WgtColumn = Coef;

// TODO: make Five
// TODO:remove Six (same as Seven)
#[derive(Copy, Clone, Debug)]
pub enum ContingencyTable {
    // AdaBoost
    // (T, F) = (True, False)
    Two((f64, f64)),
    // ConstAda
    // (TP, FP, FN, TN) = (D1, N1, D0, N0); Diseased, Not diseased
    Four((f64, f64, f64, f64)),
    // (TP, FP, FN, TN, M) = (D1, N1, D0, N0, M); Diseased, Not diseased
    Five((f64, f64, f64, f64, f64)),
    // FreeModel
    // (D2, N2, D1, N1, D0, N0)
    Six((f64, f64, f64, f64, f64, f64)),
    // FreeModelMissing
    // (D2, N2, D1, N1, D0, N0, M); Diseased, Not diseased, Missing
    Seven((f64, f64, f64, f64, f64, f64, f64)),
}

impl ContingencyTable {
    pub fn new_four(x: (f64, f64, f64, f64)) -> ContingencyTable {
        let x = ContingencyTable::Four(x);
        x.check();
        x
    }
    pub fn new_five(x: (f64, f64, f64, f64, f64)) -> ContingencyTable {
        let x = ContingencyTable::Five(x);
        x.check();
        x
    }

    pub fn new_seven(x: (f64, f64, f64, f64, f64, f64, f64)) -> ContingencyTable {
        let x = ContingencyTable::Seven(x);
        x.check();
        x
    }
    pub fn two(self) -> (f64, f64) {
        match self {
            ContingencyTable::Two(x) => x,
            _ => panic!("Not four."),
        }
    }
    pub fn four(self) -> (f64, f64, f64, f64) {
        match self {
            ContingencyTable::Four(x) => x,
            _ => panic!("Not four."),
        }
    }
    pub fn five(self) -> (f64, f64, f64, f64, f64) {
        match self {
            ContingencyTable::Five(x) => x,
            _ => panic!("Not four."),
        }
    }
    pub fn seven(self) -> (f64, f64, f64, f64, f64, f64, f64) {
        match self {
            ContingencyTable::Seven(x) => x,
            _ => panic!("Not seven."),
        }
    }

    pub fn is_five(self) -> bool {
        if let ContingencyTable::Seven(_) = self {
            true
        } else {
            false
        }
    }
    pub fn is_seven(self) -> bool {
        if let ContingencyTable::Seven(_) = self {
            true
        } else {
            false
        }
    }

    pub fn check(self) {
        self.check_all_positive();
    }

    fn check_all_positive(self) {
        let f = match self {
            Self::Two((t, f)) => (t > -MACHINE_EPS) && (f > -MACHINE_EPS),
            Self::Four((d1, n1, d0, n0)) => {
                (d1 > -MACHINE_EPS)
                    && (n1 > -MACHINE_EPS)
                    && (d0 > -MACHINE_EPS)
                    && (n0 > -MACHINE_EPS)
            }
            Self::Five((d1, n1, d0, n0, m)) => {
                (d1 > -MACHINE_EPS)
                    && (n1 > -MACHINE_EPS)
                    && (d0 > -MACHINE_EPS)
                    && (n0 > -MACHINE_EPS)
                    && (m > -MACHINE_EPS)
            }
            Self::Six((d2, n2, d1, n1, d0, n0)) => {
                (d2 > -MACHINE_EPS)
                    && (n2 > -MACHINE_EPS)
                    && (d1 > -MACHINE_EPS)
                    && (n1 > -MACHINE_EPS)
                    && (d0 > -MACHINE_EPS)
                    && (n0 > -MACHINE_EPS)
            }
            Self::Seven((d2, n2, d1, n1, d0, n0, m)) => {
                (d2 > -MACHINE_EPS)
                    && (n2 > -MACHINE_EPS)
                    && (d1 > -MACHINE_EPS)
                    && (n1 > -MACHINE_EPS)
                    && (d0 > -MACHINE_EPS)
                    && (n0 > -MACHINE_EPS)
                    && (m >= -MACHINE_EPS)
            }
        };
        if !f {
            panic!(
                "Some of ContingencyTable are negative or non-positive: {:?}",
                self
            );
        }
    }

    /// check sum ~1.0
    pub fn check_sum_table4(self) {
        let sum = match self {
            ContingencyTable::Four((x0, x1, x2, x3)) => x0 + x1 + x2 + x3,
            _ => panic!("Not four."),
        };
        if (sum - 1.0).abs() > 1e-10 {
            panic!("Sum of table4 is not 1.0, diff is: {}", (sum - 1.0).abs())
        }
    }
}

// pred_s -> GenotSnv seems better?
//      for this, Ada, ConstAda should also give si
// -> ok: now pred is simple enough
pub fn calculate_table_eps(
    pred_s: &[u8],
    //pred_s: &GenotTwinSnvRef,
    //pred_s: &[B8],
    ps: &[f64],
    phe: &Phe,
    eps: Option<Eps>,
    boost_type: BoostType,
) -> (ContingencyTable, bool) {
    match boost_type {
        BoostType::Ada => calculate_table2(pred_s, ps, phe),
        BoostType::ConstAda => calculate_table4_eps(pred_s, ps, phe, eps),
        BoostType::FreeModelMissing => calculate_table7_eps(pred_s, ps, phe, eps),
        _ => panic!(),
    }
}

pub fn calculate_table2(pred_sm: &[u8], ps: &[f64], phe: &Phe) -> (ContingencyTable, bool) {
    calculate_table2_sm(pred_sm, ps, phe)
}

/// assumed to be called by outside
pub fn calculate_table4_eps(
    pred_sm: &[u8],
    ps: &[f64],
    phe: &Phe,
    eps: Option<Eps>,
) -> (ContingencyTable, bool) {
    let epsilons = epsilon::calculate_epsilons(ps, phe, eps);

    // first calculate E:=case_ws_sum and F:=control_ws_sum
    let ef_ = calculate_ef_sum(ps, phe);

    calculate_table4_epsilons(pred_sm, ps, phe, ef_, epsilons, eps)
}

pub fn calculate_table7_eps(
    pred_s: &[u8],
    //pred_s: &GenotTwinSnvRef,
    ps: &[f64],
    phe: &Phe,
    eps: Option<Eps>,
) -> (ContingencyTable, bool) {
    let epsilons = epsilon::calculate_epsilons(ps, phe, eps);

    calculate_table7_epsilons(pred_s, ps, phe, epsilons, eps)
}

pub fn calculate_table7_epsilons(
    pred_s: &[u8],
    //pred_s: &GenotTwinSnvRef,
    ps: &[f64],
    phe: &Phe,
    epsilons: (f64, f64), // (epsilon_case: f64, epsilon_cont: f64,)
    eps: Option<Eps>,
) -> (ContingencyTable, bool) {
    //fixed_seven_sum((d2_, n2_, d1_, n1_, d0_, n0_), epsilons, eps)

    let table7_ori = calculate_table7_pred(pred_s, ps, phe);

    if eps.is_none() {
        return (table7_ori, false);
    }

    if eps.unwrap().dom() {
        // TODO: cleaner: how about basically hold original table7_ori and convert them on using ex. calculate loss and eff
        //
        // do not convert 7 to 5 since do this in wgt.calculate_cofficients()
        // otherwise, cannot distinguish concat 2/1 or 1/0
        // table7_or_5
        let (_, is_eps) = adjust_eps_table7_dom(table7_ori, eps.unwrap());
        // RETURN original table7
        (table7_ori, is_eps)
    } else {
        //table7
        adjust_eps_table7_nondom(table7_ori, epsilons, eps.unwrap())
    }
}

pub fn calculate_table7_pred(
    pred_s: &[u8],
    //pred_s: &GenotTwinSnvRef,
    ps: &[f64],
    phe: &Phe,
) -> ContingencyTable {
    //let n = phe.n();

    let mut d2_ = 0.0f64;
    let mut n2_ = 0.0f64;
    let mut d1_ = 0.0f64;
    let mut n1_ = 0.0f64;
    let mut d0_ = 0.0f64;
    let mut n0_ = 0.0f64;

    //ps[..n]
    ps.iter()
        .zip(pred_s.iter())
        .zip(phe.iter())
        .for_each(|((p, count), y)| {
            if y {
                match count {
                    2 => d2_ += *p,
                    1 => d1_ += *p,
                    0 => d0_ += *p,
                    _ => (),
                }
            } else {
                match count {
                    2 => n2_ += *p,
                    1 => n1_ += *p,
                    0 => n0_ += *p,
                    _ => (),
                }
            }
        });

    let m_sum: f64 = 1.0 - (d2_ + n2_ + d1_ + n1_ + d0_ + n0_);
    // TODO: m_sum>MACHNE_EPSILON
    let m_sum = m_sum.max(0.0);

    ContingencyTable::new_seven((d2_, n2_, d1_, n1_, d0_, n0_, m_sum))
}

/// TODO: dom ver.
pub fn adjust_eps_logit(
    wzs_sum: (f64, f64, f64),
    wls_sum: (f64, f64, f64),
    epsilons_wzs: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    epsilons_wls: (f64, f64), //(epsilon_case: f64, epsilon_cont: f64,)
    eps: Option<Eps>,
    table8_count: (usize, usize, usize, usize, usize, usize, usize, usize),
) -> ((f64, f64, f64), (f64, f64, f64), bool) {
    if eps.is_some() && eps.unwrap().dom() {
        panic!("wrong eps")
    }

    //if eps != Eps::MedLarge2 {
    //    panic!("wrong eps")
    //}

    let (wzs_sum_adj, wzs_is_eps) = match eps {
        None => (wzs_sum, false),
        Some(x) => {
            match x {
                Eps::MedLarge2 => {
                    adjust_eps_logit_nondom_add_cat(wzs_sum, epsilons_wzs, table8_count, true)
                }
                Eps::MedLargeAllCellAllSnvs => adjust_eps_logit_nondom_add_cat_allcellallsnvs(
                    wzs_sum,
                    epsilons_wzs,
                    //table8_count,
                    true,
                ),
                _ => unimplemented!(),
            }
        }
    };
    let (wls_sum_adj, wls_is_eps) = match eps {
        None => (wls_sum, false),
        Some(x) => {
            match x {
                Eps::MedLarge2 => {
                    adjust_eps_logit_nondom_add_cat(wls_sum, epsilons_wls, table8_count, false)
                }
                Eps::MedLargeAllCellAllSnvs => adjust_eps_logit_nondom_add_cat_allcellallsnvs(
                    wls_sum,
                    epsilons_wls,
                    //table8_count,
                    false,
                ),
                //Eps::None => (wls_sum, false),
                _ => unimplemented!(),
            }
        }
    };

    //let (wzs_sum_adj, wzs_is_eps) =
    //    adjust_eps_logit_nondom_add_cat(wzs_sum, epsilons_wzs, table8_count, true);
    //let (wls_sum_adj, wls_is_eps) =
    //    adjust_eps_logit_nondom_add_cat(wls_sum, epsilons_wls, table8_count, false);

    assert_eq!(wzs_is_eps, wls_is_eps);

    (wzs_sum_adj, wls_sum_adj, wzs_is_eps)
}

fn adjust_eps_logit_nondom_add_cat(
    wzs_sum: (f64, f64, f64),
    epsilons: (f64, f64),
    table8_count: (usize, usize, usize, usize, usize, usize, usize, usize),
    negative_for_cont: bool,
) -> ((f64, f64, f64), bool) {
    let (wzs_sum2, wzs_sum1, wzs_sum0) = wzs_sum;
    //let (mut wzs_sum2, mut wzs_sum1, mut wzs_sum0)=wzs_sum;
    //let (epsilon_case, epsilon_cont) = epsilons_wzs;

    let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;

    fn add_eps_sum(
        wzs_sum: f64,
        d: usize,
        n: usize,
        epsilons: (f64, f64),
        negative_for_cont: bool,
    ) -> (f64, bool) {
        let mut wzs_sum_ = wzs_sum;
        let (epsilon_case, epsilon_cont) = epsilons;

        let mut is_eps = false;

        // if d==n, wzs_sum_=0.0 for negative
        if d == 0 {
            if negative_for_cont {
                // wzs: wzs<0
                wzs_sum_ = (wzs_sum_ + epsilon_case).min(0.0);
                //wzs_sum_ = wzs_sum_ + epsilon_case;
            } else {
                // wls: wls>0
                wzs_sum_ = wzs_sum_ + epsilon_case;
            }
            is_eps = true
        }
        if n == 0 {
            if negative_for_cont {
                // wzs: wzs>0
                wzs_sum_ = (wzs_sum_ - epsilon_cont).max(0.0);
                //wzs_sum_ = wzs_sum_ - epsilon_cont;
            } else {
                // wls: wls>0
                wzs_sum_ = wzs_sum_ + epsilon_cont;
            }
            is_eps = true
        }

        (wzs_sum_, is_eps)
    }

    let (wzs_sum2_, is_eps2) = add_eps_sum(wzs_sum2, d2, n2, epsilons, negative_for_cont);
    let (wzs_sum1_, is_eps1) = add_eps_sum(wzs_sum1, d1, n1, epsilons, negative_for_cont);
    let (wzs_sum0_, is_eps0) = add_eps_sum(wzs_sum0, d0, n0, epsilons, negative_for_cont);

    let is_eps = is_eps2 | is_eps1 | is_eps0;

    ((wzs_sum2_, wzs_sum1_, wzs_sum0_), is_eps)
}

fn adjust_eps_logit_nondom_add_cat_allcellallsnvs(
    wzs_sum: (f64, f64, f64),
    epsilons: (f64, f64),
    //table8_count: (usize, usize, usize, usize, usize, usize, usize, usize),
    negative_for_cont: bool,
) -> ((f64, f64, f64), bool) {
    let (wzs_sum2, wzs_sum1, wzs_sum0) = wzs_sum;
    //let (mut wzs_sum2, mut wzs_sum1, mut wzs_sum0)=wzs_sum;
    //let (epsilon_case, epsilon_cont) = epsilons_wzs;

    //let (d2, n2, d1, n1, d0, n0, _dm, _nm) = table8_count;

    fn add_eps_sum(wzs_sum: f64, epsilons: (f64, f64), negative_for_cont: bool) -> (f64, bool) {
        let mut wzs_sum_ = wzs_sum;
        let (epsilon_case, _) = epsilons;

        //let mut is_eps = false;

        if negative_for_cont {
            // for wzs
            if wzs_sum_ > 0.0 {
                wzs_sum_ = (wzs_sum_ - epsilon_case).max(0.0);
            } else {
                wzs_sum_ = (wzs_sum_ + epsilon_case).min(0.0);
            }
        } else {
            // for wls
            // assume eps_case=eps_cont
            wzs_sum_ = wzs_sum_ + epsilon_case;
        }

        /*
        // if d==n, wzs_sum_=0.0 for negative
        if d == 0 {
            // wzs<0.0
            wzs_sum_ = wzs_sum_ + epsilon_case;
            is_eps = true
        }
        if n == 0 {
            // negative for wzs
            // positive for wls
            if negative_for_cont {
                wzs_sum_ = wzs_sum_ - epsilon_cont;
            } else {
                wzs_sum_ = wzs_sum_ + epsilon_cont;
            }
            is_eps = true
        }
         */

        (wzs_sum_, true)
    }

    let (wzs_sum2_, is_eps2) = add_eps_sum(wzs_sum2, epsilons, negative_for_cont);
    let (wzs_sum1_, is_eps1) = add_eps_sum(wzs_sum1, epsilons, negative_for_cont);
    let (wzs_sum0_, is_eps0) = add_eps_sum(wzs_sum0, epsilons, negative_for_cont);

    let is_eps = is_eps2 | is_eps1 | is_eps0;

    ((wzs_sum2_, wzs_sum1_, wzs_sum0_), is_eps)
}

pub fn adjust_eps_table7_dom(table7: ContingencyTable, eps: Eps) -> (ContingencyTable, bool) {
    if !eps.dom() {
        panic!("wrong eps")
    }

    //let mut is_eps = false;

    if let ContingencyTable::Seven(ws_sum) = table7 {
        let (d2, n2, d1, n1, d0, n0, m) = ws_sum;

        if (d2 < CONTINGENCY_TABLE_FILL) || (n2 < CONTINGENCY_TABLE_FILL) {
            // concat d2/d1 and n2/n1
            let d1_new = d2 + d1;
            let n1_new = n2 + n1;
            return (
                ContingencyTable::new_five((d1_new, n1_new, d0, n0, m)),
                true,
            );
        } else if (d0 < CONTINGENCY_TABLE_FILL) || (n0 < CONTINGENCY_TABLE_FILL) {
            // concat d0/d1 and n0/n1
            let d1_new = d0 + d1;
            let n1_new = n0 + n1;
            return (
                ContingencyTable::new_five((d2, n2, d1_new, n1_new, m)),
                true,
            );
        }

        // assume d1, n1 are never 0.0

        return (table7, false);
    } else {
        panic!("Not Seven.");
    }
}

/// make missing >=0.0
/// should run before adjust_eps
fn adjust_missing_table7(table7: ContingencyTable) -> ContingencyTable {
    let (d2, n2, d1, n1, d0, n0, _m) = table7.seven();

    let mut m = 1.0 - (d2 + n2 + d1 + n1 + d0 + n0);
    m = m.max(0.0);

    ContingencyTable::new_seven((d2, n2, d1, n1, d0, n0, m))
}

// TODO:  integrate with dom
// TODO:  integrate with table4
pub fn adjust_eps_table7_nondom(
    //ws_sum: (f64, f64, f64, f64, f64, f64),
    table7: ContingencyTable,
    epsilons: (f64, f64),
    eps: Eps,
) -> (ContingencyTable, bool) {
    if eps.dom() {
        panic!("wrong eps")
    }

    if !table7.is_seven() {
        panic!("Not Seven.");
    }

    // TODO: cleaner
    //if let ContingencyTable::Seven(ws_sum) = table7 {
    //let ws_sum_6=(ws_sum.0,ws_sum.1,ws_sum.2,ws_sum.3,ws_sum.4,ws_sum.5);
    //fixed_seven_sum(ws_sum_6,epsilons,eps)
    //ContingencyTable::Seven(ws_sum_6) => fixed_sven_sum(ws_sum,epsilons,eps),

    // make m>=0.0
    let table7 = adjust_missing_table7(table7);

    if eps.allsnv() {
        if eps.allcell() {
            // add for all snvs and to allcell
            return adjust_eps_table7_nondom_allsnv_allcell(table7, epsilons);
        } else {
            unimplemented!()
        }
    } else {
        if eps.allcell() {
            // add to d2~n0
            return adjust_eps_table7_nondom_allcell(table7, epsilons);
        } else {
            return adjust_eps_table7_nondom_add_cat(table7, epsilons);
        }
    }
}

fn adjust_eps_table7_nondom_allsnv_allcell(
    table7: ContingencyTable,
    epsilons: (f64, f64),
) -> (ContingencyTable, bool) {
    let (mut d2, mut n2, mut d1, mut n1, mut d0, mut n0, m) = table7.seven();
    let (epsilon_case, epsilon_cont) = epsilons;

    d2 = d2 + epsilon_case;
    d1 = d1 + epsilon_case;
    d0 = d0 + epsilon_case;
    n2 = n2 + epsilon_cont;
    n1 = n1 + epsilon_cont;
    n0 = n0 + epsilon_cont;
    let is_eps = true;

    (
        ContingencyTable::new_seven((d2, n2, d1, n1, d0, n0, m)),
        is_eps,
    )
}

fn adjust_eps_table7_nondom_allcell(
    table7: ContingencyTable,
    epsilons: (f64, f64),
) -> (ContingencyTable, bool) {
    let (mut d2, mut n2, mut d1, mut n1, mut d0, mut n0, m) = table7.seven();
    let (epsilon_case, epsilon_cont) = epsilons;

    let mut is_eps = false;

    if (d2 < CONTINGENCY_TABLE_FILL)
        || (n2 < CONTINGENCY_TABLE_FILL)
        || (d1 < CONTINGENCY_TABLE_FILL)
        || (n1 < CONTINGENCY_TABLE_FILL)
        || (d0 < CONTINGENCY_TABLE_FILL)
        || (n0 < CONTINGENCY_TABLE_FILL)
    {
        d2 = d2 + epsilon_case;
        d1 = d1 + epsilon_case;
        d0 = d0 + epsilon_case;
        n2 = n2 + epsilon_cont;
        n1 = n1 + epsilon_cont;
        n0 = n0 + epsilon_cont;
        is_eps = true;
    }

    (
        ContingencyTable::new_seven((d2, n2, d1, n1, d0, n0, m)),
        is_eps,
    )
}

fn adjust_eps_table7_nondom_add_cat(
    table7: ContingencyTable,
    epsilons: (f64, f64),
) -> (ContingencyTable, bool) {
    let (mut d2, mut n2, mut d1, mut n1, mut d0, mut n0, m) = table7.seven();
    let (epsilon_case, epsilon_cont) = epsilons;
    let mut is_eps = false;

    if (d2 < CONTINGENCY_TABLE_FILL) || (n2 < CONTINGENCY_TABLE_FILL) {
        d2 = d2 + epsilon_case;
        n2 = n2 + epsilon_cont;
        is_eps = true
    }
    if (d1 < CONTINGENCY_TABLE_FILL) || (n1 < CONTINGENCY_TABLE_FILL) {
        d1 = d1 + epsilon_case;
        n1 = n1 + epsilon_cont;
        is_eps = true
    }
    if (d0 < CONTINGENCY_TABLE_FILL) || (n0 < CONTINGENCY_TABLE_FILL) {
        d0 = d0 + epsilon_case;
        n0 = n0 + epsilon_cont;
        is_eps = true
    }

    (
        ContingencyTable::new_seven((d2, n2, d1, n1, d0, n0, m)),
        is_eps,
    )
}

/*
pub fn fixed_seven_sum(
    ws_sum: (f64, f64, f64, f64, f64, f64),
    epsilons: (f64, f64),
    eps: Eps,
) -> (ContingencyTable, bool) {
    // input a,b and return appropriate abcd

    let (epsilon_case, epsilon_cont) = epsilons;
    let (mut d2, mut n2, mut d1, mut n1, mut d0, mut n0) = ws_sum;
    let mut m = 1.0 - (d2 + n2 + d1 + n1 + d0 + n0);
    m = m.max(0.0);

    let mut is_eps = false;

    if eps.dom() {
        panic!("should not be here.")
    }

    if eps.allsnv() {
        if eps.allcell() {
            // add for all snvs and to allcell

            let ((d2, n2, d1, n1, d0, n0), is_eps) =
                adjust_eps_table7_nondom_allsnv_allcell(d2, n2, d1, n1, d0, n0);

            //d2 = d2 + epsilon_case;
            //d1 = d1 + epsilon_case;
            //d0 = d0 + epsilon_case;
            //n2 = n2 + epsilon_cont;
            //n1 = n1 + epsilon_cont;
            //n0 = n0 + epsilon_cont;
            //is_eps = true;
        } else {
            unimplemented!()
        }
    } else {
        if eps.allcell() {
            // add to d2~n0

            let ((d2, n2, d1, n1, d0, n0), is_eps) =
                adjust_eps_table7_nondom_allcell(d2, n2, d1, n1, d0, n0);

            /*
            if (d2 < CONTINGENCY_TABLE_FILL)
                || (n2 < CONTINGENCY_TABLE_FILL)
                || (d1 < CONTINGENCY_TABLE_FILL)
                || (n1 < CONTINGENCY_TABLE_FILL)
                || (d0 < CONTINGENCY_TABLE_FILL)
                || (n0 < CONTINGENCY_TABLE_FILL)
            {
                d2 = d2 + epsilon_case;
                d1 = d1 + epsilon_case;
                d0 = d0 + epsilon_case;
                n2 = n2 + epsilon_cont;
                n1 = n1 + epsilon_cont;
                n0 = n0 + epsilon_cont;
                is_eps = true;
            }
             */
        } else {
            let ((d2, n2, d1, n1, d0, n0), is_eps) =
                adjust_eps_table7_nondom_add_cat(d2, n2, d1, n1, d0, n0);

            /*
            if (d2 < CONTINGENCY_TABLE_FILL) || (n2 < CONTINGENCY_TABLE_FILL) {
                d2 = d2 + epsilon_case;
                n2 = n2 + epsilon_cont;
                is_eps = true
            }
            if (d1 < CONTINGENCY_TABLE_FILL) || (n1 < CONTINGENCY_TABLE_FILL) {
                d1 = d1 + epsilon_case;
                n1 = n1 + epsilon_cont;
                is_eps = true
            }
            if (d0 < CONTINGENCY_TABLE_FILL) || (n0 < CONTINGENCY_TABLE_FILL) {
                d0 = d0 + epsilon_case;
                n0 = n0 + epsilon_cont;
                is_eps = true
            }
             */
        }
    }

    return (
        ContingencyTable::new_seven((d2, n2, d1, n1, d0, n0, m)),
        is_eps,
    );
}
 */

// TODO: this is not bad but not consistent
// input table4 and return eps table4
pub fn adjust_eps_table4(
    ab: (f64, f64),
    ef: (f64, f64),
    epsilons: (f64, f64),
    eps: Option<Eps>,
) -> (ContingencyTable, bool) {
    // input a,b and return appropriate abcd

    let (epsilon_case, epsilon_cont) = epsilons;
    let (mut a, mut b) = ab;
    let (e, f) = ef;

    let mut c = e - a;
    let mut d = f - b;

    let mut is_eps = false;

    if eps.is_none() {
        return (ContingencyTable::new_four((a, b, c, d)), is_eps);
    }

    if eps.unwrap().allsnv() {
        unimplemented!()
    }

    if eps.unwrap().allcell() {
        if (a < CONTINGENCY_TABLE_FILL)
            || (b < CONTINGENCY_TABLE_FILL)
            || (c < CONTINGENCY_TABLE_FILL)
            || (d < CONTINGENCY_TABLE_FILL)
        {
            a = a + epsilon_case;
            c = c + epsilon_case;
            b = b + epsilon_cont;
            d = d + epsilon_cont;
            is_eps = true;
        }
    } else {
        if (a < CONTINGENCY_TABLE_FILL) || (b < CONTINGENCY_TABLE_FILL) {
            a = a + epsilon_case;
            b = b + epsilon_cont;
            is_eps = true
        }
        if (c < CONTINGENCY_TABLE_FILL) || (d < CONTINGENCY_TABLE_FILL) {
            c = c + epsilon_case;
            d = d + epsilon_cont;
            is_eps = true
        }
    }

    return (ContingencyTable::new_four((a, b, c, d)), is_eps);
}

pub fn calculate_table4_epsilons(
    pred_sm: &[u8],
    ps: &[f64],
    phe: &Phe,
    ef_: (f64, f64),
    epsilons: (f64, f64), // (epsilon_case: f64, epsilon_cont: f64,)
    eps: Option<Eps>,
) -> (ContingencyTable, bool) {
    //let n = phe.n();

    let mut a_ = 0.0f64;
    let mut b_ = 0.0f64;
    //for (ni, p) in ps[..n].iter().enumerate() {
    for (ni, p) in ps.iter().enumerate() {
        let pred = pred_sm[ni] != 0;
        //let pred = bit::bget(pred_sm, ni);
        //let mis = operate::bget(pred_sm, ni);
        //log::debug!("y,mis,p: {},{},{}", y, mis, *p);
        if !pred {
            continue;
        }

        let y = phe.get_unchecked(ni);
        //let y = bit::bget(ys, ni);
        if y {
            a_ += *p;
        } else {
            b_ += *p;
        }
        /*
        if y && pred {
            a_ += *p;
        }
        // faster to use if not else if?
        if (!y) && pred {
            b_ += *p;
        }
        */
    }

    //unimplemented!("eps");

    adjust_eps_table4((a_, b_), ef_, epsilons, eps)
}

/// For T/F,
/// TODO: eps
fn calculate_table2_sm(pred_sm: &[u8], ps: &[f64], phe: &Phe) -> (ContingencyTable, bool) {
    //let n = phe.n();
    // a_=correct
    let mut a_ = 0.0f64;

    //ps[..n]
    ps.iter()
        .zip(pred_sm.iter())
        .zip(phe.iter())
        .for_each(|((p, pred), y)| {
            if y == (*pred != 0) {
                a_ += *p;
            }
        });

    let b_ = 1.0 - a_;

    (ContingencyTable::Two((a_, b_)), false)

    //unimplemented!("eps");
}

pub fn calculate_ef_sum(ps: &[f64], phe: &Phe) -> (f64, f64) {
    //let n = phe.n();
    //let e_ = ps[..n]
    let e_ = ps
        .iter()
        .zip(phe.iter())
        .filter(|(_, ph)| *ph)
        .map(|(p, _)| *p)
        .sum();

    let f_ = 1.0 - e_;
    (e_, f_)
}

/// calculate loss for FreeModelMissing
///  loss = M + 2*(sqrt(D2 * N2)+sqrt(D1 * N1)+sqrt(D0 *N0))
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
pub unsafe fn calculate_table7_sum_simd(
    gsnv: &GenotSnvRef,
    ps_pad: &[f64],
    phe: &Phe,
) -> ContingencyTable {
    // aarch64 is completely different: https://moshg.github.io/rust-std-ja/core/arch/aarch64/index.html
    //#[cfg(target_arch = "aarch64")]
    //use std::arch::aarch64::*;
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;
    //use std::convert::TryInto;

    let n = phe.n();
    let ys = phe.inner();
    let pred_s0m = gsnv.predict_s(0);
    let pred_s1m = gsnv.predict_s(1);

    // b: 0x80808080 = -x: value as i32
    // ~b = x - 1
    // x = ~b + 1
    //log::debug!("0x7f7f7f7f + 1: {:#010x}", 0x7f7f7f7f + 1 as i32);
    //log::debug!("0x7f7f7f7f + 1: {}", 0x7f7f7f7f + 1 as i32);
    //log::debug!("0x7f7f7f80: {}", 0x7f7f7f80 as i32);

    let bit_ext_mask: __m256i = _mm256_set_epi32(
        0x10101010,
        0x01010101,
        0x20202020,
        0x02020202,
        0x40404040,
        0x04040404,
        -0x7f7f7f80, // =0x80808080
        0x08080808,
    );
    let zerod: __m256d = _mm256_setzero_pd();
    let zeros: __m256i = _mm256_setzero_si256();
    let ones: __m256i = _mm256_cmpeq_epi32(zeros, zeros);

    let mut d2_sum_acc = _mm256_setzero_pd();
    let mut n2_sum_acc = _mm256_setzero_pd();
    let mut d1_sum_acc = _mm256_setzero_pd();
    let mut n1_sum_acc = _mm256_setzero_pd();
    let mut d0_sum_acc = _mm256_setzero_pd();
    let mut n0_sum_acc = _mm256_setzero_pd();

    // bi=0-3, shift=24,16,8,0
    //let shift_v = _mm256_set1_epi32(shift);
    let shifts: [__m256i; 4] = [
        _mm256_set1_epi32(24),
        _mm256_set1_epi32(16),
        _mm256_set1_epi32(8),
        _mm256_set1_epi32(0),
    ];

    // when n=33, (n/32+1)=2
    // sample index to look through is [0..64)
    // size of y: n/8+5=9: [0..64), which is included
    for ni in 0..(n / 32 + 1) {
        //log::debug!("ni {}", ni);

        // broadcast 32bit int to 256bit
        // ex. DCBA -> DCBADCBA...DCBA
        // (D=abcdefgh)

        // 1. use _mm256_set_epi8
        // o2. use from_be() and use set1 <- assume to be fast??
        let pred_s0_b32 = u32::from_le_bytes(pred_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let pred_s1_b32 = u32::from_le_bytes(pred_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let mis_s0_b32 = u32::from_be_bytes(mis_s0m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        //let mis_s1_b32 = u32::from_be_bytes(mis_s1m[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let predv_s0_32 = _mm256_set1_epi32(pred_s0_b32 as i32);
        let predv_s1_32 = _mm256_set1_epi32(pred_s1_b32 as i32);

        //log::debug!("mis 0: {:#010b}", mis_s0m[4 * ni]);
        //log::debug!("mis 1: {:#010b}", mis_s0m[4 * ni + 1]);
        //log::debug!("mis 2: {:#010b}", mis_s0m[4 * ni + 2]);
        //log::debug!("mis 3: {:#010b}", mis_s0m[4 * ni + 3]);
        //log::debug!("mis b32: {:#034b}", mis_s0_b32);
        //log::debug!("misv {:?}", misv_s0_32);

        let ys_b32 = u32::from_le_bytes(ys[4 * ni..4 * (ni + 1)].try_into().unwrap());
        let yv_32 = _mm256_set1_epi32(ys_b32 as i32);

        // bitwise not is fastest in xor by xor(v, ones)
        // D2sum = y & pred0 & pred1 = y & ( pred0 & pred1 )
        // N2sum = !y & pred0 & pred1 = !y & ( pred0 & pred1 )
        // D1sum = y & pred0 & !pred1 = y & ( !pred1 & pred0 )
        // N1sum = !y & pred0 & !pred1 = !y & ( !pred1 & pred0 )
        // D0sum = y & !pred0 & !pred1 = !pred0 & ( !pred1 & y )
        // N0sum = !y & !pred0 & !pred1 = !y & ( !pred0 & (!pred1) )
        let flagv_d2_32 = _mm256_and_si256(yv_32, _mm256_and_si256(predv_s0_32, predv_s1_32));
        let flagv_n2_32 = _mm256_andnot_si256(yv_32, _mm256_and_si256(predv_s0_32, predv_s1_32));
        let flagv_d1_32 = _mm256_and_si256(yv_32, _mm256_andnot_si256(predv_s1_32, predv_s0_32));
        let flagv_n1_32 = _mm256_andnot_si256(yv_32, _mm256_andnot_si256(predv_s1_32, predv_s0_32));
        let flagv_d0_32 = _mm256_andnot_si256(predv_s0_32, _mm256_andnot_si256(predv_s1_32, yv_32));
        let flagv_n0_32 = _mm256_andnot_si256(
            yv_32,
            _mm256_andnot_si256(predv_s0_32, _mm256_xor_si256(predv_s1_32, ones)),
        );

        // ex. D=abcdefgh -> extract d,h,c,g,b,f,a,e for each 32 bit
        // abcdefgh(x4)|...
        // -> extracted (at highest position)
        // 000d0000(x4)|...
        // -> mask
        // 11111111|00000000|00000000|11111111|...
        let flagv_d2_32_ext = _mm256_and_si256(flagv_d2_32, bit_ext_mask);
        let flagv_n2_32_ext = _mm256_and_si256(flagv_n2_32, bit_ext_mask);
        let flagv_d1_32_ext = _mm256_and_si256(flagv_d1_32, bit_ext_mask);
        let flagv_n1_32_ext = _mm256_and_si256(flagv_n1_32, bit_ext_mask);
        let flagv_d0_32_ext = _mm256_and_si256(flagv_d0_32, bit_ext_mask);
        let flagv_n0_32_ext = _mm256_and_si256(flagv_n0_32, bit_ext_mask);

        let take_mask_d2_32 = _mm256_cmpeq_epi8(flagv_d2_32_ext, bit_ext_mask);
        let take_mask_n2_32 = _mm256_cmpeq_epi8(flagv_n2_32_ext, bit_ext_mask);
        let take_mask_d1_32 = _mm256_cmpeq_epi8(flagv_d1_32_ext, bit_ext_mask);
        let take_mask_n1_32 = _mm256_cmpeq_epi8(flagv_n1_32_ext, bit_ext_mask);
        let take_mask_d0_32 = _mm256_cmpeq_epi8(flagv_d0_32_ext, bit_ext_mask);
        let take_mask_n0_32 = _mm256_cmpeq_epi8(flagv_n0_32_ext, bit_ext_mask);

        // bi=0-3, shift=24,16,8,0
        //const SHIFTS: &'static [i32] = &[24, 16, 8, 0];

        for bi in 0usize..4 {
            // DCBADCBA...DCBA
            // -> b=1: for B
            // BA00BA00...BA00

            let shift_v = shifts[bi];
            let take_mask_d2 = _mm256_sllv_epi32(take_mask_d2_32, shift_v);
            let take_mask_n2 = _mm256_sllv_epi32(take_mask_n2_32, shift_v);
            let take_mask_d1 = _mm256_sllv_epi32(take_mask_d1_32, shift_v);
            let take_mask_n1 = _mm256_sllv_epi32(take_mask_n1_32, shift_v);
            let take_mask_d0 = _mm256_sllv_epi32(take_mask_d0_32, shift_v);
            let take_mask_n0 = _mm256_sllv_epi32(take_mask_n0_32, shift_v);

            //log::debug!("take_mask a s0 {:?}", take_mask_a_s0);
            //log::debug!("take_mask b s0 {:?}", take_mask_b_s0);
            //log::debug!("take_mask a s1 {:?}", take_mask_a_s1);
            //log::debug!("take_mask b s1 {:?}", take_mask_b_s1);

            let psv_lo_ptr = ps_pad[32 * ni + 8 * bi..32 * ni + 8 * bi + 4].as_ptr();
            let psv_hi_ptr = ps_pad[32 * ni + 8 * bi + 4..32 * ni + 8 * (bi + 1)].as_ptr();

            //log::debug!("ps ind {}", 32 * ni + 8 * bi);
            //log::debug!("ps ind {}", 32 * ni + 8 * bi + 4);
            //log::debug!("ps lo {:?}", &ps[32 * ni + 8 * bi..32 * ni + 8 * bi + 4]);
            //log::debug!(
            //    "ps hi {:?}",
            //    &ps[32 * ni + 8 * bi + 4..32 * ni + 8 * bi + 8]
            //);

            let psv_lo: __m256d = _mm256_load_pd(psv_lo_ptr as *const _);
            let psv_hi: __m256d = _mm256_load_pd(psv_hi_ptr as *const _);

            //log::debug!("ps lo {:?}", psv_lo);
            //log::debug!("ps hi {:?}", psv_hi);

            // first for low
            let ps_masked_d2_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_d2));
            let ps_masked_n2_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_n2));
            let ps_masked_d1_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_d1));
            let ps_masked_n1_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_n1));
            let ps_masked_d0_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_d0));
            let ps_masked_n0_lo =
                _mm256_blendv_pd(zerod, psv_lo, _mm256_castsi256_pd(take_mask_n0));

            d2_sum_acc = _mm256_add_pd(d2_sum_acc, ps_masked_d2_lo);
            n2_sum_acc = _mm256_add_pd(n2_sum_acc, ps_masked_n2_lo);
            d1_sum_acc = _mm256_add_pd(d1_sum_acc, ps_masked_d1_lo);
            n1_sum_acc = _mm256_add_pd(n1_sum_acc, ps_masked_n1_lo);
            d0_sum_acc = _mm256_add_pd(d0_sum_acc, ps_masked_d0_lo);
            n0_sum_acc = _mm256_add_pd(n0_sum_acc, ps_masked_n0_lo);

            //log::debug!("ps a s0 lo {:?}", ps_masked_a_s0_lo);
            //log::debug!("ps a s1 lo {:?}", ps_masked_a_s1_lo);
            //log::debug!("ps b s0 lo {:?}", ps_masked_b_s0_lo);
            //log::debug!("ps b s1 lo {:?}", ps_masked_b_s1_lo);

            // for high
            let take_mask_d2_hi = _mm256_slli_epi64(take_mask_d2, 32);
            let take_mask_n2_hi = _mm256_slli_epi64(take_mask_n2, 32);
            let take_mask_d1_hi = _mm256_slli_epi64(take_mask_d1, 32);
            let take_mask_n1_hi = _mm256_slli_epi64(take_mask_n1, 32);
            let take_mask_d0_hi = _mm256_slli_epi64(take_mask_d0, 32);
            let take_mask_n0_hi = _mm256_slli_epi64(take_mask_n0, 32);

            let ps_masked_d2_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_d2_hi));
            let ps_masked_n2_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_n2_hi));
            let ps_masked_d1_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_d1_hi));
            let ps_masked_n1_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_n1_hi));
            let ps_masked_d0_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_d0_hi));
            let ps_masked_n0_hi =
                _mm256_blendv_pd(zerod, psv_hi, _mm256_castsi256_pd(take_mask_n0_hi));

            //log::debug!("a s0 hi {:?}", ps_masked_a_s0_hi);

            d2_sum_acc = _mm256_add_pd(d2_sum_acc, ps_masked_d2_hi);
            n2_sum_acc = _mm256_add_pd(n2_sum_acc, ps_masked_n2_hi);
            d1_sum_acc = _mm256_add_pd(d1_sum_acc, ps_masked_d1_hi);
            n1_sum_acc = _mm256_add_pd(n1_sum_acc, ps_masked_n1_hi);
            d0_sum_acc = _mm256_add_pd(d0_sum_acc, ps_masked_d0_hi);
            n0_sum_acc = _mm256_add_pd(n0_sum_acc, ps_masked_n0_hi);
        }
    }

    // sum 4 double horizontally to get the whole sum
    d2_sum_acc = _mm256_hadd_pd(d2_sum_acc, d2_sum_acc);
    n2_sum_acc = _mm256_hadd_pd(n2_sum_acc, n2_sum_acc);
    d1_sum_acc = _mm256_hadd_pd(d1_sum_acc, d1_sum_acc);
    n1_sum_acc = _mm256_hadd_pd(n1_sum_acc, n1_sum_acc);
    d0_sum_acc = _mm256_hadd_pd(d0_sum_acc, d0_sum_acc);
    n0_sum_acc = _mm256_hadd_pd(n0_sum_acc, n0_sum_acc);

    // 1. any way to hadd??
    // 2. _mm256_extractf128_pd and _mm256_cvtsd_f64: get 64:0

    let d2_sum: f64 =
        _mm256_cvtsd_f64(d2_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(d2_sum_acc, 1));
    let n2_sum: f64 =
        _mm256_cvtsd_f64(n2_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(n2_sum_acc, 1));
    let d1_sum: f64 =
        _mm256_cvtsd_f64(d1_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(d1_sum_acc, 1));
    let n1_sum: f64 =
        _mm256_cvtsd_f64(n1_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(n1_sum_acc, 1));
    let d0_sum: f64 =
        _mm256_cvtsd_f64(d0_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(d0_sum_acc, 1));
    let n0_sum: f64 =
        _mm256_cvtsd_f64(n0_sum_acc) + _mm_cvtsd_f64(_mm256_extractf128_pd(n0_sum_acc, 1));

    let m_sum: f64 = 1.0 - (d2_sum + n2_sum + d1_sum + n1_sum + d0_sum + n0_sum);
    // TODO: check m_sum>-MACHNE_EPSILON
    let m_sum = m_sum.max(0.0);

    ContingencyTable::new_seven((d2_sum, n2_sum, d1_sum, n1_sum, d0_sum, n0_sum, m_sum))
}

#[cfg(test)]
mod tests {
    use super::*;

    //fn is_eq_f64(v: f64, w: f64, e: f64) -> bool {
    //    (v - w).abs() < e
    //}

    #[test]
    fn test_fixed_table4_sum() {
        // no eps
        let ab = (0.1, 0.2);
        let ef = (0.4, 0.6);
        let epsilons = (0.05, 0.05);
        let (t, b) = adjust_eps_table4(ab, ef, epsilons, Some(Eps::MedLarge2AllCell));
        //assert_eq!(abcd, (0.1, 0.2, 0.3, 0.4));
        assert_eq!(b, false);
        let (a, b, c, d) = t.four();
        assert_float_absolute_eq!(a, 0.1);
        assert_float_absolute_eq!(b, 0.2);
        assert_float_absolute_eq!(c, 0.3);
        assert_float_absolute_eq!(d, 0.4);
    }

    #[test]
    fn test_fixed_table4_sum_eps() {
        // use eps
        let ab = (0.4, 0.0);
        let ef = (0.5, 0.5);
        let epsilons = (0.05, 0.05);
        let (t, b) = adjust_eps_table4(ab, ef, epsilons, Some(Eps::MedLarge2));
        // abcd ori (0.4, 0.0, 0.1, 0.5)
        // abcd (0.45, 0.05, 0.1, 0.5)
        assert_eq!(b, true);
        let (a, b, c, d) = t.four();
        assert_float_absolute_eq!(a, 0.45);
        assert_float_absolute_eq!(b, 0.05);
        assert_float_absolute_eq!(c, 0.1);
        assert_float_absolute_eq!(d, 0.5);
    }

    #[test]
    fn test_fixed_table4_sum_eps_allcell() {
        // use eps
        let ab = (0.4, 0.0);
        let ef = (0.5, 0.5);
        let epsilons = (0.05, 0.05);
        let (t, b) = adjust_eps_table4(ab, ef, epsilons, Some(Eps::MedLarge2AllCell));
        // abcd ori (0.4, 0.0, 0.1, 0.5)
        // abcd (0.45, 0.05, 0.15, 0.55)
        assert_eq!(b, true);
        let (a, b, c, d) = t.four();
        assert_float_absolute_eq!(a, 0.45);
        assert_float_absolute_eq!(b, 0.05);
        assert_float_absolute_eq!(c, 0.15);
        assert_float_absolute_eq!(d, 0.55);
    }

    #[test]
    fn test_fixed_table7_sum() {
        // no eps
        //let ws_sum = (0.02, 0.01, 0.1, 0.2, 0.3, 0.3);
        let table7 = ContingencyTable::new_seven((0.02, 0.01, 0.1, 0.2, 0.3, 0.3, 0.07));
        let epsilons = (0.05, 0.05);
        let (t, b) = adjust_eps_table7_nondom(table7, epsilons, Eps::MedLarge2AllCell);
        assert_eq!(b, false);
        let (d2, n2, d1, n1, d0, n0, m) = t.seven();
        assert_float_absolute_eq!(d2, 0.02);
        assert_float_absolute_eq!(n2, 0.01);
        assert_float_absolute_eq!(d1, 0.1);
        assert_float_absolute_eq!(n1, 0.2);
        assert_float_absolute_eq!(d0, 0.3);
        assert_float_absolute_eq!(n0, 0.3);
        assert_float_absolute_eq!(m, 0.07);
    }

    #[test]
    fn test_fixed_table7_sum_eps() {
        // use eps
        //let ws_sum = (0.0, 0.01, 0.1, 0.2, 0.3, 0.3);
        let table7 = ContingencyTable::new_seven((0.0, 0.01, 0.1, 0.2, 0.3, 0.3, 0.09));
        let epsilons = (0.05, 0.05);
        let (t, b) = adjust_eps_table7_nondom(table7, epsilons, Eps::MedLarge2);
        assert_eq!(b, true);
        let (d2, n2, d1, n1, d0, n0, m) = t.seven();
        assert_float_absolute_eq!(d2, 0.05);
        assert_float_absolute_eq!(n2, 0.06);
        assert_float_absolute_eq!(d1, 0.1);
        assert_float_absolute_eq!(n1, 0.2);
        assert_float_absolute_eq!(d0, 0.3);
        assert_float_absolute_eq!(n0, 0.3);
        assert_float_absolute_eq!(m, 0.09);
    }

    #[test]
    fn test_fixed_table7_sum_eps_allcell() {
        // use eps
        //let ws_sum = (0.0, 0.01, 0.1, 0.2, 0.3, 0.3);
        let table7 = ContingencyTable::new_seven((0.0, 0.01, 0.1, 0.2, 0.3, 0.3, 0.09));
        let epsilons = (0.05, 0.05);
        let (t, b) = adjust_eps_table7_nondom(table7, epsilons, Eps::MedLarge2AllCell);
        assert_eq!(b, true);
        let (d2, n2, d1, n1, d0, n0, m) = t.seven();
        assert_float_absolute_eq!(d2, 0.05);
        assert_float_absolute_eq!(n2, 0.06);
        assert_float_absolute_eq!(d1, 0.15);
        assert_float_absolute_eq!(n1, 0.25);
        assert_float_absolute_eq!(d0, 0.35);
        assert_float_absolute_eq!(n0, 0.35);
        assert_float_absolute_eq!(m, 0.09);
    }

    #[test]
    fn test_calculate_ef_sum() {
        //let ps = [0.1, 0.3, 0.25, 0.05, 0.3, 0.0, 0.0];
        let ps = [0.1, 0.3, 0.25, 0.05, 0.3];
        let phe_v = vec![false, false, true, true, true];
        let phe = Phe::new(&phe_v);
        let (e_, f_) = calculate_ef_sum(&ps, &phe);
        assert_eq!(e_, 0.6);
        assert_eq!(f_, 0.4);
    }
}
