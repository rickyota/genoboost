//! Split SnvIndex from Snv because SnvIndex will be used for several purposes and easy for implement Eq.

pub mod chrom;
pub mod io;
pub mod snv_index;
pub mod snvs;

pub use chrom::Chrom;
pub use snv_index::SnvId;
pub use snvs::Snvs;

pub use io::*;

pub fn extract_snvs_consume(snvs_in: Vec<SnvId>, use_snvs: &[bool], m: usize) -> Vec<SnvId> {
    let mut snvs: Vec<SnvId> = Vec::with_capacity(m);

    for (snv, b) in snvs_in.into_iter().zip(use_snvs.iter()) {
        if *b {
            snvs.push(snv);
        }
    }

    snvs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_snvs_consume() {
        let snvs_in = vec![
            SnvId::new(
                "rs1".to_string(),
                "1",
                "1",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs2".to_string(),
                "1",
                "2",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs3".to_string(),
                "1",
                "3",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs4".to_string(),
                "1",
                "4",
                "A".to_string(),
                "C".to_string(),
            ),
            SnvId::new(
                "rs5".to_string(),
                "1",
                "5",
                "A".to_string(),
                "C".to_string(),
            ),
        ];

        let use_snvs = vec![true, false, true, false, true];
        let m = 3;

        let snvs = extract_snvs_consume(snvs_in, &use_snvs, m);

        assert_eq!(snvs.len(), m);
        assert_eq!(snvs[0].id(), "rs1");
        assert_eq!(snvs[1].id(), "rs3");
        assert_eq!(snvs[2].id(), "rs5");
    }
}
