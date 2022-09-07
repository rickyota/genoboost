pub mod io;

use super::{Chrom, SnvId};

#[derive(Clone, Debug, Default)]
pub struct SumStat {
    snv_index: SnvId,
    loss: Option<f64>,
    alpha: Option<f64>,
}

impl SumStat {
    pub fn construct_sum_stat_string(
        rs: String,
        chrom: &str,
        pos: &str,
        a1: String,
        a2: String,
        loss: &str,
        alpha: &str,
    ) -> SumStat {
        let snv = SumStat {
            snv_index: SnvId::construct_snv_index_string(rs, chrom, pos, a1, a2),
            loss: Some(loss.parse::<f64>().unwrap()),
            alpha: Some(alpha.parse::<f64>().unwrap()),
        };
        snv
    }

    pub fn construct_sum_stat_from_snv_string(
        snv_index: SnvId,
        loss: Option<f64>,
        alpha: Option<f64>,
    ) -> SumStat {
        let snv = SumStat {
            snv_index,
            loss,
            alpha,
        };
        snv
    }

    pub fn set_loss(&mut self, loss_ss: Option<f64>, alpha_ss: Option<f64>) {
        self.loss = loss_ss;
        self.alpha = alpha_ss;
    }

    /*
    pub fn snv(&self) -> &Snv {
        &self.snv
    }

    pub fn snv_mut(&mut self) -> &mut Snv {
        &mut self.snv
    }
     */

    pub fn snv_index(&self) -> &SnvId {
        &self.snv_index
    }

    /*
    pub fn snv_index_mut(&mut self) -> &mut SnvIndex {
        self.snv_mut().snv_index_mut()
    }
    */
    pub fn rs(&self) -> &str {
        self.snv_index().rs()
    }

    pub fn chrom(&self) -> &Chrom {
        self.snv_index().chrom()
    }

    pub fn pos(&self) -> usize {
        self.snv_index().pos()
    }

    pub fn a1(&self) -> &str {
        self.snv_index().a1()
    }

    pub fn a2(&self) -> &str {
        self.snv_index().a2()
    }

    pub fn sida(&self) -> &str {
        self.snv_index().sida()
    }

    pub fn loss(&self) -> Option<f64> {
        self.loss
    }
    pub fn alpha(&self) -> Option<f64> {
        self.alpha
    }
}

impl std::fmt::Display for SumStat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.sida())
    }
}

impl PartialEq for SumStat {
    fn eq(&self, other: &Self) -> bool {
        self.snv_index().eq(other.snv_index())
    }
}

impl Eq for SumStat {}

impl PartialEq<SnvId> for SumStat {
    fn eq(&self, other: &SnvId) -> bool {
        self.snv_index().eq(other)
    }
}

impl PartialEq<SumStat> for SnvId {
    fn eq(&self, other: &SumStat) -> bool {
        self.eq(other.snv_index())
    }
}

/*
impl PartialEq<Snv> for SumStat {
    fn eq(&self, other: &Snv) -> bool {
        self.snv_index().eq(other)
    }
}

impl PartialEq<SumStat> for Snv {
    fn eq(&self, other: &SumStat) -> bool {
        self.eq(other.snv_index())
    }
}
 */

impl PartialOrd for SumStat {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SumStat {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.snv_index().cmp(other.snv_index())
    }
}

impl AsRef<SnvId> for SumStat {
    #[inline]
    fn as_ref(&self) -> &SnvId {
        self.snv_index()
    }
}

/*
impl AsMut<SnvIndex> for SumStat {
    #[inline]
    fn as_mut(&mut self) -> &mut SnvIndex {
        self.snv_index_mut()
    }
}
*/

/*
pub fn extract_snvs_consume(snvs_in: Vec<Snv>, use_snvs: &[bool], m: usize) -> Vec<Snv> {
    let mut snvs: Vec<Snv> = Vec::with_capacity(m);

    for (snv, b) in snvs_in.into_iter().zip(use_snvs.iter()) {
        if *b {
            snvs.push(snv);
        }
    }

    snvs
}

pub fn extract_snvs(snvs_in: &[Snv], use_snvs: &[bool], m: usize) -> Vec<Snv> {
    let mut snvs: Vec<Snv> = Vec::with_capacity(m);

    for (snv, b) in snvs_in.iter().zip(use_snvs.iter()) {
        if *b {
            snvs.push(snv.clone());
        }
    }

    snvs
}
*/

/*
#[cfg(test)]
mod tests {
}
*/
