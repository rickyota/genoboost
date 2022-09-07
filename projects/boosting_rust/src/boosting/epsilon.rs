use super::sample_weight;
use crate::Eps;
use genetics::samples::prelude::*;

fn med_ps_large(ps: &[f64], ys: &Phe, is_large: bool) -> f64 {
    let (med_case, med_cont) = (
        sample_weight::med_ps(ps, ys, true),
        sample_weight::med_ps(ps, ys, false),
    );
    if is_large {
        med_case.max(med_cont)
    } else {
        med_case.min(med_cont)
    }
}

fn calculate_epsilon_medlarge2allcell(ps: &[f64], ys: &Phe) -> f64 {
    let med_large = med_ps_large(ps, ys, true);

    med_large / 2.0
}

fn calculate_epsilon_medsmall(ps: &[f64], ys: &Phe) -> f64 {
    let med_small = med_ps_large(ps, ys, false);

    med_small
}

/// return (epsilon_case, epsilon_cont)
pub fn calculate_epsilons(ps: &[f64], ys: &Phe, eps: Eps) -> (f64, f64) {
    match eps {
        // large ~ eps of case
        Eps::MedLarge2AllCell => {
            let e = calculate_epsilon_medlarge2allcell(ps, ys);
            (e, e)
        }
        // small ~ eps of cont
        Eps::MedSmall => {
            let e = calculate_epsilon_medsmall(ps, ys);
            (e, e)
        }
        _ => unimplemented!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_epsilons_medlarge2allcell() {
        let phe_v = vec![true, true, true, false];
        let phe = Phe::new(phe_v);
        let ps = vec![0.1, 0.2, 0.3, 0.4, 0.0, 0.0];
        let eps = calculate_epsilons(&ps, &phe, Eps::MedLarge2AllCell);
        // med_case=0.2, cont=0.4
        // eps_large=0.4
        // eps=0.2
        assert_eq!(eps, (0.2, 0.2));
    }

    #[test]
    fn test_calculate_epsilons_medsmall() {
        let phe_v = vec![true, true, true, false];
        let phe = Phe::new(phe_v);
        let ps = vec![0.1, 0.2, 0.3, 0.4, 0.0, 0.0];
        let eps = calculate_epsilons(&ps, &phe, Eps::MedLarge2AllCell);
        // med_case=0.2, cont=0.4
        // eps_small=0.2
        assert_eq!(eps, (0.2, 0.2));
    }
}
