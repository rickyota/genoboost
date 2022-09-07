use std::{collections::HashMap, path::Path};

use crate::{text, SnvId, SumStat};

//use io_rust::{text, Snv};

// TODO: variable col
pub fn load_loss(fin_loss: &Path) -> Vec<SumStat> {
    text::check_open_file(fin_loss);

    let m = text::compute_num_line(fin_loss).unwrap();

    let mut sum_stat: Vec<SumStat> = Vec::with_capacity(m);

    // rs, chrom, pos, A1, A2, loss,alpha
    let cols = [7usize, 0, 1, 2, 3, 10, 9];
    //let cols = [8usize, 1, 0, 3, 4, 5, 11, 10];
    let vss: Vec<Vec<String>> = text::load_table_cols(&fin_loss, &cols, true).unwrap();

    for vi in 0..vss[0].len() {
        sum_stat.push(SumStat::construct_sum_stat_string(
            vss[0][vi].clone(),
            &vss[1][vi],
            &vss[2][vi],
            vss[3][vi].clone(),
            vss[4][vi].clone(),
            &vss[5][vi],
            &vss[6][vi],
        ));
    }

    sum_stat
}

// TODO: what if fin_loss does not contain some SNVs?
// load loss_ss and alpha_ss and set to snvs
pub fn convert_sum_stat_from_snv_consume(fin_loss: &Path, snvs: Vec<SnvId>) -> Vec<SumStat> {
    let snvs_loss = load_loss(fin_loss);

    // map: sid_in -> plink index
    let m = snvs.len();
    let mut sida_to_index: HashMap<String, usize> = HashMap::with_capacity(m);
    //let mut sida_to_index: HashMap<String, usize> = HashMap::with_capacity(m);

    // This .to_owned() is necessary
    for (si, s) in snvs_loss.iter().enumerate() {
        sida_to_index.insert(s.sida().to_owned(), si);
    }

    let mut sum_stats: Vec<SumStat> = Vec::with_capacity(snvs.len());

    for snv in snvs.into_iter() {
        let sum_stat: SumStat;

        let sida = snv.sida();
        if let Some(&v) = sida_to_index.get(sida) {
            let snv_loss = &snvs_loss[v];
            let loss = snv_loss.loss();
            let alpha = snv_loss.alpha();
            //let mut snv = &snvs[*v];
            //(*snv).set_loss_ss(loss, alpha);
            //snvs[v].set_loss_ss(loss, alpha);
            sum_stat = SumStat::construct_sum_stat_from_snv_string(snv, loss, alpha);
        } else {
            sum_stat = SumStat::construct_sum_stat_from_snv_string(snv, None, None);
        }

        sum_stats.push(sum_stat);
    }
    sum_stats
}

/*
#[cfg(test)]
mod tests {
    use super::*;
}
 */
