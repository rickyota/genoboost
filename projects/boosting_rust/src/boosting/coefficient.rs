use super::BoostType;
use crate::ContingencyTable;
use genetics::wgt::Coef;

// assume table after eps
pub fn calculate_coefficients(
    table: ContingencyTable,
    boost_type: BoostType,
    learning_rate: Option<f64>,
) -> Coef {
    let lr = learning_rate.unwrap_or(1.0);
    match boost_type {
        BoostType::Ada => {
            let ab_sum = table.two();
            let (a, b) = ab_sum;
            let alpha = (a / b).ln() / 2.0;

            let alpha = lr * alpha;
            Coef::Binary((0.0, alpha))
            //Coef::Binary((alpha, 0.0))
        }
        BoostType::ConstAda => {
            let abcd_sum = table.four();
            let (a, b, c, d) = abcd_sum;
            let const_ti = ((a * c) / (b * d)).ln() / 4.0;
            let alpha_ti = ((a * d) / (b * c)).ln() / 4.0;

            let const_ti = lr * const_ti;
            let alpha_ti = lr * alpha_ti;
            Coef::Binary((const_ti, alpha_ti))
            //Coef::Binary((alpha_ti, const_ti))
        }
        BoostType::FreeModelMissing => {
            let table7_sum = table.seven();
            let (d2, n2, d1, n1, d0, n0, _) = table7_sum;
            let s0 = (d0 / n0).ln() / 2.0;
            let s1 = (d1 / n1).ln() / 2.0;
            let s2 = (d2 / n2).ln() / 2.0;

            let s0 = lr * s0;
            let s1 = lr * s1;
            let s2 = lr * s2;
            Coef::Score4((s0, s1, s2, 0.0))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_coefficients_freemodelmissing() {
        let t = ContingencyTable::new_seven((0.02, 0.01, 0.1, 0.2, 0.3, 0.3, 0.07));
        let coef = calculate_coefficients(t, BoostType::FreeModelMissing, None);
        assert_eq!(
            coef,
            Coef::Score4((0.0, 0.5f64.ln() / 2.0, 2.0f64.ln() / 2.0, 0.0))
        );
    }

    #[test]
    fn test_calculate_coefficients_freemodelmissing_lr() {
        let t = ContingencyTable::new_seven((0.02, 0.01, 0.1, 0.2, 0.3, 0.3, 0.07));
        let coef = calculate_coefficients(t, BoostType::FreeModelMissing, Some(0.1));
        assert_eq!(
            coef,
            Coef::Score4((0.0, 0.5f64.ln() / 20.0, 2.0f64.ln() / 20.0, 0.0))
        );
    }
}
