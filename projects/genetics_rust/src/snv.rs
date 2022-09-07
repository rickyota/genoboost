//! Split SnvIndex from Snv because SnvIndex will be used for several purposes and easy for implement Eq.

pub mod chrom;
pub mod snv_index;

pub use chrom::Chrom;
pub use snv_index::SnvId;

pub fn extract_snvs_consume(snvs_in: Vec<SnvId>, use_snvs: &[bool], m: usize) -> Vec<SnvId> {
    let mut snvs: Vec<SnvId> = Vec::with_capacity(m);

    for (snv, b) in snvs_in.into_iter().zip(use_snvs.iter()) {
        if *b {
            snvs.push(snv);
        }
    }

    snvs
}

/*
#[cfg(test)]
mod tests {
}
*/
