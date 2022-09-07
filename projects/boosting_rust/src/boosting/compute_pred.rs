use super::WgtBoost;
use genetics::cov::CovKind;
use genetics::genot::BaseGenot;
use genetics::wgt::WgtKind;
use genetics::Genot;

pub fn compute_pred_const(pred: &mut [u8]) {
    // true as u8
    pred.fill(1);
}

pub fn compute_pred(pred: &mut [u8], wgt: &WgtBoost, predictions: &Genot) {
    match wgt.wgt().kind() {
        WgtKind::Snv(_, si, mi) => {
            let mi = mi.unwrap();
            match si {
                Some(si) => {
                    let v = predictions.vals_snv_s(mi, *si);
                    pred.iter_mut()
                        .zip(v.iter())
                        .for_each(|(p, b)| *p = *b as u8)
                }
                // modelfree
                None => {
                    let v = predictions.vals_snv(mi);
                    pred.copy_from_slice(&v);
                }
            }
        }
        WgtKind::Cov(cov_wgt) => match cov_wgt.kind() {
            CovKind::Const => {
                compute_pred_const(pred);
            }
            CovKind::Cov => {
                panic!("ny");
            }
        },
    }
}

/*
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_coefficients() {}
}
*/
