// TODO: when fin does not exist, raise error or return Option<usize>
// return Option<usize> might be better
// TODO: PathBuf -> T: AsRef<Path>

pub mod fileio;
pub mod text;
pub mod tsv;

pub use fileio::*;
pub use text::*;
//pub use tsv::*;
