use super::epsilon;
use crate::{BoostType, Eps};
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
    pub fn seven(self) -> (f64, f64, f64, f64, f64, f64, f64) {
        match self {
            ContingencyTable::Seven(x) => x,
            _ => panic!("Not seven."),
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
    pub fn check_sum_abcd(self) {
        let sum = match self {
            ContingencyTable::Four((x0, x1, x2, x3)) => x0 + x1 + x2 + x3,
            _ => panic!("Not four."),
        };
        if (sum - 1.0).abs() > 1e-10 {
            panic!("Sum of abcd is not 1.0, diff is: {}", (sum - 1.0).abs())
        }
    }
}

// pred_s -> GenotSnv seems better?
//      for this, Ada, ConstAda should also give si
// -> ok: now pred is simple enough
pub fn calculate_table_sum(
    pred_s: &[u8],
    //pred_s: &GenotTwinSnvRef,
    //pred_s: &[B8],
    ps: &[f64],
    ys: &Phe,
    eps: Eps,
    boost_type: BoostType,
) -> (ContingencyTable, bool) {
    match boost_type {
        BoostType::Ada => calculate_table2_sum(pred_s, ps, ys),
        BoostType::ConstAda => calculate_table4_sum(pred_s, ps, ys, eps),
        BoostType::FreeModelMissing => calculate_table7_sum(pred_s, ps, ys, eps),
    }
}

pub fn calculate_table2_sum(pred_sm: &[u8], ps: &[f64], ys: &Phe) -> (ContingencyTable, bool) {
    calculate_ab_sum_sm(pred_sm, ps, ys)
}

/// assumed to be called by outside
pub fn calculate_table4_sum(
    pred_sm: &[u8],
    ps: &[f64],
    ys: &Phe,
    eps: Eps,
) -> (ContingencyTable, bool) {
    let epsilons = epsilon::calculate_epsilons(ps, ys, eps);

    // first calculate E:=case_ws_sum and F:=control_ws_sum
    let ef_ = calculate_ef_sum(ps, ys);

    calculate_abcd_sum_sm(pred_sm, ps, ys, ef_, epsilons, eps)
}

pub fn calculate_table7_sum(
    pred_s: &[u8],
    //pred_s: &GenotTwinSnvRef,
    ps: &[f64],
    ys: &Phe,
    eps: Eps,
) -> (ContingencyTable, bool) {
    let epsilons = epsilon::calculate_epsilons(ps, ys, eps);

    calculate_table7_sum_sm(pred_s, ps, ys, epsilons, eps)
}

pub fn fixed_seven_sum(
    ws_sum: (f64, f64, f64, f64, f64, f64),
    epsilons: (f64, f64),
    eps: Eps,
) -> (ContingencyTable, bool) {
    // input a,b and return appropriate abcd

    let (epsilon_case, epsilon_cont) = epsilons;
    let (mut d2, mut n2, mut d1, mut n1, mut d0, mut n0) = ws_sum;
    // This makes is_eps less chosen
    let mut m = 1.0 - (d2 + n2 + d1 + n1 + d0 + n0);
    m = m.max(0.0);

    let mut is_eps = false;

    if eps.allcell() {
        // add to d2~n0
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
    } else {
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
    }

    return (
        ContingencyTable::new_seven((d2, n2, d1, n1, d0, n0, m)),
        is_eps,
    );
}

pub fn fixed_abcd_sum(
    ab: (f64, f64),
    ef: (f64, f64),
    epsilons: (f64, f64),
    eps: Eps,
) -> (ContingencyTable, bool) {
    // input a,b and return appropriate abcd

    let (epsilon_case, epsilon_cont) = epsilons;
    let (mut a, mut b) = ab;
    let (e, f) = ef;

    let mut c = e - a;
    let mut d = f - b;

    let mut is_eps = false;

    if eps.allcell() {
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

pub fn calculate_abcd_sum_sm(
    pred_sm: &[u8],
    ps: &[f64],
    ys: &Phe,
    ef_: (f64, f64),
    epsilons: (f64, f64), // (epsilon_case: f64, epsilon_cont: f64,)
    eps: Eps,
) -> (ContingencyTable, bool) {
    let n = ys.n();

    let mut a_ = 0.0f64;
    let mut b_ = 0.0f64;
    for (ni, p) in ps[..n].iter().enumerate() {
        let pred = pred_sm[ni] != 0;
        //let pred = bit::bget(pred_sm, ni);
        //let mis = operate::bget(pred_sm, ni);
        //println!("y,mis,p: {},{},{}", y, mis, *p);
        if !pred {
            continue;
        }

        let y = ys.get_unchecked(ni);
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

    fixed_abcd_sum((a_, b_), ef_, epsilons, eps)
}

pub fn calculate_table7_sum_sm(
    pred_s: &[u8],
    //pred_s: &GenotTwinSnvRef,
    ps: &[f64],
    ys: &Phe,
    epsilons: (f64, f64), // (epsilon_case: f64, epsilon_cont: f64,)
    eps: Eps,
) -> (ContingencyTable, bool) {
    let n = ys.n();

    let mut d2_ = 0.0f64;
    let mut n2_ = 0.0f64;
    let mut d1_ = 0.0f64;
    let mut n1_ = 0.0f64;
    let mut d0_ = 0.0f64;
    let mut n0_ = 0.0f64;

    ps[..n]
        .iter()
        .zip(pred_s.iter())
        .zip(ys.iter())
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

    fixed_seven_sum((d2_, n2_, d1_, n1_, d0_, n0_), epsilons, eps)
}

/// For T/F,
/// TODO: eps
fn calculate_ab_sum_sm(pred_sm: &[u8], ps: &[f64], ys: &Phe) -> (ContingencyTable, bool) {
    let n = ys.n();
    // a_=correct
    let mut a_ = 0.0f64;

    ps[..n]
        .iter()
        .zip(pred_sm.iter())
        .zip(ys.iter())
        .for_each(|((p, pred), y)| {
            if y == (*pred != 0) {
                a_ += *p;
            }
        });

    let b_ = 1.0 - a_;

    (ContingencyTable::Two((a_, b_)), false)

    //unimplemented!("eps");
}

pub fn calculate_ef_sum(ps: &[f64], ys: &Phe) -> (f64, f64) {
    let n = ys.n();
    let e_ = ps[..n]
        .iter()
        .zip(ys.iter())
        .filter(|(_, ph)| *ph)
        .map(|(p, _)| *p)
        .sum();

    let f_ = 1.0 - e_;
    (e_, f_)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn is_eq_f64(v: f64, w: f64, e: f64) -> bool {
        (v - w).abs() < e
    }

    #[test]
    fn test_fixed_abcd_sum() {
        // no eps
        let ab = (0.1, 0.2);
        let ef = (0.4, 0.6);
        let epsilons = (0.05, 0.05);
        let (t, b) = fixed_abcd_sum(ab, ef, epsilons, Eps::MedLarge2AllCell);
        //assert_eq!(abcd, (0.1, 0.2, 0.3, 0.4));
        assert_eq!(b, false);
        let (a, b, c, d) = t.four();
        assert!(is_eq_f64(a, 0.1, 1e-10));
        assert!(is_eq_f64(b, 0.2, 1e-10));
        assert!(is_eq_f64(c, 0.3, 1e-10));
        assert!(is_eq_f64(d, 0.4, 1e-10));
    }

    #[test]
    fn test_fixed_abcd_sum_eps() {
        // use eps
        let ab = (0.4, 0.0);
        let ef = (0.5, 0.5);
        let epsilons = (0.05, 0.05);
        let (t, b) = fixed_abcd_sum(ab, ef, epsilons, Eps::MedLarge2);
        // abcd ori (0.4, 0.0, 0.1, 0.5)
        // abcd (0.45, 0.05, 0.1, 0.5)
        assert_eq!(b, true);
        let (a, b, c, d) = t.four();
        assert!(is_eq_f64(a, 0.45, 1e-10));
        assert!(is_eq_f64(b, 0.05, 1e-10));
        assert!(is_eq_f64(c, 0.1, 1e-10));
        assert!(is_eq_f64(d, 0.5, 1e-10));
    }

    #[test]
    fn test_fixed_abcd_sum_eps_allcell() {
        // use eps
        let ab = (0.4, 0.0);
        let ef = (0.5, 0.5);
        let epsilons = (0.05, 0.05);
        let (t, b) = fixed_abcd_sum(ab, ef, epsilons, Eps::MedLarge2AllCell);
        // abcd ori (0.4, 0.0, 0.1, 0.5)
        // abcd (0.45, 0.05, 0.15, 0.55)
        assert_eq!(b, true);
        let (a, b, c, d) = t.four();
        assert!(is_eq_f64(a, 0.45, 1e-10));
        assert!(is_eq_f64(b, 0.05, 1e-10));
        assert!(is_eq_f64(c, 0.15, 1e-10));
        assert!(is_eq_f64(d, 0.55, 1e-10));
    }

    #[test]
    fn test_fixed_table7_sum() {
        // no eps
        let ws_sum = (0.02, 0.01, 0.1, 0.2, 0.3, 0.3);
        let epsilons = (0.05, 0.05);
        let (t, b) = fixed_seven_sum(ws_sum, epsilons, Eps::MedLarge2AllCell);
        assert_eq!(b, false);
        let (d2, n2, d1, n1, d0, n0, m) = t.seven();
        assert!(is_eq_f64(d2, 0.02, 1e-10));
        assert!(is_eq_f64(n2, 0.01, 1e-10));
        assert!(is_eq_f64(d1, 0.1, 1e-10));
        assert!(is_eq_f64(n1, 0.2, 1e-10));
        assert!(is_eq_f64(d0, 0.3, 1e-10));
        assert!(is_eq_f64(n0, 0.3, 1e-10));
        assert!(is_eq_f64(m, 0.07, 1e-10));
    }

    #[test]
    fn test_fixed_table7_sum_eps() {
        // use eps
        let ws_sum = (0.0, 0.01, 0.1, 0.2, 0.3, 0.3);
        let epsilons = (0.05, 0.05);
        let (t, b) = fixed_seven_sum(ws_sum, epsilons, Eps::MedLarge2);
        assert_eq!(b, true);
        let (d2, n2, d1, n1, d0, n0, m) = t.seven();
        assert!(is_eq_f64(d2, 0.05, 1e-10));
        assert!(is_eq_f64(n2, 0.06, 1e-10));
        assert!(is_eq_f64(d1, 0.1, 1e-10));
        assert!(is_eq_f64(n1, 0.2, 1e-10));
        assert!(is_eq_f64(d0, 0.3, 1e-10));
        assert!(is_eq_f64(n0, 0.3, 1e-10));
        assert!(is_eq_f64(m, 0.09, 1e-10));
    }

    #[test]
    fn test_fixed_table7_sum_eps_allcell() {
        // use eps
        let ws_sum = (0.0, 0.01, 0.1, 0.2, 0.3, 0.3);
        let epsilons = (0.05, 0.05);
        let (t, b) = fixed_seven_sum(ws_sum, epsilons, Eps::MedLarge2AllCell);
        assert_eq!(b, true);
        let (d2, n2, d1, n1, d0, n0, m) = t.seven();
        assert!(is_eq_f64(d2, 0.05, 1e-10));
        assert!(is_eq_f64(n2, 0.06, 1e-10));
        assert!(is_eq_f64(d1, 0.15, 1e-10));
        assert!(is_eq_f64(n1, 0.25, 1e-10));
        assert!(is_eq_f64(d0, 0.35, 1e-10));
        assert!(is_eq_f64(n0, 0.35, 1e-10));
        assert!(is_eq_f64(m, 0.09, 1e-10));
    }

    #[test]
    fn test_calculate_ef_sum() {
        let ps = [0.1, 0.3, 0.25, 0.05, 0.3, 0.0, 0.0];
        let phe_v = vec![false, false, true, true, true];
        let phe = Phe::new(phe_v);
        let (e_, f_) = calculate_ef_sum(&ps, &phe);
        assert_eq!(e_, 0.6);
        assert_eq!(f_, 0.4);
    }
}
