
//! Compressed Matrix crate
//!
//! When using base,  use as `use cmatrix::prelude::*`
//!
//! There are three main structures: [CMatrix], [CVec], and [CBits].
//! CMatrix stores row x column of one or two bits of values.
//! CVec stores a column of CMatrix.
//! CBits stores a bit vector of CVec.
//!
//! Each structures have reference and mutable reference structure, e.g., CMatrixRef and CMatrixMut.
//!
//! There are two traits for each structure, basic trait and mutable trait, e.g., BaseCMatrix and BaseCMatrixMut.
//! CMatrix and CMatrixMut implement both, and CMatrixRef implement BaseCMatrix only.
//!
//! All six traits are independent like below.
//!
//!
//! BaseCMatrix -> BaseCMatrixMut
//!         V                            V
//! BaseCVec     ->  BaseCVecMut
//!         V                            V
//! BaseCBits     ->  BaseCBitsMut
//!
//!
//! This means that CBitsMut implementing BaseCBitsMut can use all functions in all other traits.
//!
//!

//mod alloc;
pub mod base_cbits;
pub mod base_cmatrix;
pub mod base_cvec;
pub mod cbits;
pub mod cbits_mut;
pub mod cbits_ref;
pub mod cmatrix_struct;
pub mod cmatrix_mut;
pub mod cmatrix_ref;
pub mod cvec;
pub mod calc;
pub mod cvec_mut;
pub mod cvec_ref;
pub mod iterator;
pub mod prelude;

// used from inner crate
pub use base_cbits::{BaseCBits, BaseCBitsMut};
pub use base_cmatrix::{BaseCMatrix, BaseCMatrixMut};
pub use base_cvec::{BaseCVec, BaseCVecMut};
pub use cbits::CBits;
pub use cbits_mut::CBitsMut;
pub use cbits_ref::CBitsRef;
pub use cmatrix_struct::CMatrix;
pub use cmatrix_mut::CMatrixMut;
pub use cmatrix_ref::CMatrixRef;
pub use cvec::CVec;
pub use cvec_mut::CVecMut;
pub use cvec_ref::CVecRef;
pub use iterator::{BoolIter, RowIter, RowIterMut, ScalIter};
pub use calc::*;

pub type B8 = u8;
