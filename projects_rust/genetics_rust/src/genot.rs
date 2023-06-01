//! Genotype
//!
//!
//! use iterator as possible since this can avoid index boundary check.
//!

pub mod base_genot;
pub mod genot_iterator;
pub mod genot_struct;
pub mod prelude;

pub use base_genot::{BaseGenot, BaseGenotMut, BaseGenotSnv, BaseGenotSnvMut};
pub use genot_iterator::GenotIterMut;
pub use genot_struct::{Genot, GenotMut, GenotSnv, GenotSnvMut, GenotSnvRef};

pub use genot_struct::{B8, B8_2};
