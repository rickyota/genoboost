// change arg name from v to ar, and for v,i

pub mod alloc;
pub mod cov;
pub mod dataset;
pub mod genot;
pub mod genot_index;
pub mod plink;
pub mod samples;
pub mod snv;
pub mod snvs;
pub mod sum_stat;
pub mod text;
pub mod vec;
pub mod wgt;

// you can call Snv with `use io_rust::Snv` instead of `use io_rust::snv::Snv`
// should call with `use io_rust::wgt::Model` to avoid confusion
// -> no meaning
//pub use wgt::{CovWgt, Model, SnvWgt, Wgt, WgtKind};
// only make pub of Wgt
pub use cov::{CovKind, Var};
pub use dataset::Dataset;
pub use genot::genot_struct;
pub use genot::genot_struct::{Genot, B8, B8_2, THRESHOLD_SNV};
pub use samples::{CovId, Covs, Samples};
pub use snv::{Chrom, SnvId};
pub use snvs::Snvs;
pub use sum_stat::SumStat;
pub use wgt::{CovWgt, Model, ModelType, SnvWgt, Wgt, WgtKind, WgtTrait};
