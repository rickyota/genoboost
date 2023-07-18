//! TODO: when fin does not exist, raise error or return Option<usize>
//! later might be better
// TODO: PathBuf -> T: AsRef<Path>

pub mod text;
pub mod fileio;
pub mod tsv;

pub use text::*;
pub use tsv::*;
pub use fileio::*;


