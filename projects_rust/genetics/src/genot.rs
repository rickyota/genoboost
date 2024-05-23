//! Genotype
//!
//!
//! use iterator as possible since this can avoid index boundary check.
//!

pub mod base_genot;
pub mod genot_iterator;
pub mod genot_struct;
pub mod prelude;

//pub use base_genot::{BaseGenot, BaseGenotMut, BaseGenotSnv, BaseGenotSnvMut};
//pub use genot_iterator::GenotIterMut;
//pub use genot_struct::{Genot, GenotMut, GenotSnv, GenotSnvMut, GenotSnvRef};
//pub use genot_struct::{B8, B8_2};
pub use prelude::*;

pub fn maf_to_mode(maf: f64) -> u8 {
    if maf < 1.0 / 3.0f64 {
        0
    } else if maf > 2.0 / 3.0f64 {
        2
    } else {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_maf_to_mode() {
        let mode = maf_to_mode(0.1);
        assert_eq!(mode, 0);

        let mode = maf_to_mode(0.5);
        assert_eq!(mode, 1);

        let mode = maf_to_mode(0.9);
        assert_eq!(mode, 2)
    }
}
