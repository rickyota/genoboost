use super::sample_weight;
use crate::Eps;
use genetics::samples::prelude::*;

//TODO: create Epsilons
// if case==cont, just indicate one value to avoid error

fn calculate_epsilon_medlarge2_abs(ps: &[f64], phe: &Phe) -> f64 {
    let med_large = med_ps_large_abs(ps, phe, true);

    med_large / 2.0
}

fn calculate_epsilon_medlarge_abs(ps: &[f64], phe: &Phe) -> f64 {
    let med_large = med_ps_large_abs(ps, phe, true);

    med_large
}


// for Logit
fn med_ps_large_abs(ps: &[f64], phe: &Phe, is_large: bool) -> f64 {
    let (med_case, med_cont) = (
        sample_weight::med_ps(ps, phe, true),
        sample_weight::med_ps(ps, phe, false),
    );
    if is_large {
        med_case.abs().max(med_cont.abs())
    } else {
        med_case.abs().min(med_cont.abs())
    }
}

fn med_ps_large(ps: &[f64], phe: &Phe, is_large: bool) -> f64 {
    let (med_case, med_cont) = (
        sample_weight::med_ps(ps, phe, true),
        sample_weight::med_ps(ps, phe, false),
    );
    if is_large {
        med_case.max(med_cont)
    } else {
        med_case.min(med_cont)
    }
}

fn calculate_epsilon_medlargeallcell(ps: &[f64], phe: &Phe) -> f64 {
    let med_large = med_ps_large(ps, phe, true);

    med_large
}

fn calculate_epsilon_medlarge2allcell(ps: &[f64], phe: &Phe) -> f64 {
    let med_large = med_ps_large(ps, phe, true);

    med_large / 2.0
}

fn calculate_epsilon_medsmall(ps: &[f64], phe: &Phe) -> f64 {
    let med_small = med_ps_large(ps, phe, false);

    med_small
}

/// return (epsilon_case, epsilon_cont)
pub fn calculate_epsilons(ps: &[f64], phe: &Phe, eps: Option<Eps>) -> (f64, f64) {
    // all right?
    if eps.is_none(){
        return (0.0, 0.0);
    } 

    match eps.unwrap() {
        Eps::MedLargeAllCellAllSnvs => {
            let e = calculate_epsilon_medlargeallcell(ps, phe);
            (e, e)
        }
        // large ~ eps of case
        Eps::MedLarge2AllCell | Eps::MedLarge2AllCellAllSnvs => {
            let e = calculate_epsilon_medlarge2allcell(ps, phe);
            (e, e)
        }
        // small ~ eps of cont
        Eps::MedSmall => {
            let e = calculate_epsilon_medsmall(ps, phe);
            (e, e)
        }
        // never used
        Eps::Dom => (0.0, 0.0),
        _ => unimplemented!(),
    }
}


/// BOTH are POSITIVE
pub fn calc_epsilons_logit_wzs(wzs: &[f64], phe: &Phe, eps: Option<Eps>) -> (f64, f64) {

    if eps.is_none(){
        return (0.0, 0.0);
    } 

    match eps.unwrap() {
        Eps::MedLarge2 => {
            let e = calculate_epsilon_medlarge2_abs(wzs, phe);
            (e, e)
        }
        Eps::MedLargeAllCellAllSnvs => {
            let e = calculate_epsilon_medlarge_abs(wzs, phe);
            (e, e)
        }
/*         Eps::None=>{
            (0.0f64,0.0f64)
        } */
        _ => unimplemented!(),
    }
}

pub fn calc_epsilons_logit_wls(wls: &[f64], phe: &Phe, eps: Option<Eps>) -> (f64, f64) {
    if eps.is_none(){
        return (0.0, 0.0);
    } 

    match eps.unwrap() {
        Eps::MedLarge2 => {
            let e = calculate_epsilon_medlarge2allcell(wls, phe);
            (e, e)
        },
        Eps::MedLargeAllCellAllSnvs=>{
            let e = calculate_epsilon_medlargeallcell(wls, phe);
            (e, e)
        }
/*         Eps::None=>{
            (0.0f64,0.0f64)
        } */

        _ => unimplemented!(),
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_epsilons_medlarge2allcell() {
        let phe_v = vec![true, true, true, false];
        let phe = Phe::new(&phe_v);
        let ps = vec![0.1, 0.2, 0.3, 0.4];
        //let ps = vec![0.1, 0.2, 0.3, 0.4, 0.0, 0.0,0.0,0.0];
        let eps = calculate_epsilons(&ps, &phe, Some(Eps::MedLarge2AllCell));
        // med_case=0.2, cont=0.4
        // eps_large=0.4
        // eps=0.2
        assert_eq!(eps, (0.2, 0.2));
    }

    #[test]
    fn test_calculate_epsilons_medsmall() {
        let phe_v = vec![true, true, true, false];
        let phe = Phe::new(&phe_v);
        let ps = vec![0.1, 0.2, 0.3, 0.4];
        let eps = calculate_epsilons(&ps, &phe, Some(Eps::MedLarge2AllCell));
        // med_case=0.2, cont=0.4
        // eps_small=0.2
        assert_eq!(eps, (0.2, 0.2));
    }
}
