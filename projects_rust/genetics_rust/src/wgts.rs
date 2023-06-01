
use super::{wgt,Wgt};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Wgts {
    wgts: Vec<Wgt>,
    //columns: Vec<String>,
    //index_next_write: usize,
}

impl Wgts {
    pub fn new_from_file(fwgt: &Path) -> Self {
        let wgts: Vec<Wgt> = wgt::io::load_wgts_file(fwgt);

        // TODO: check colums match io::wgt_columns()

        Wgts {
            wgts,
        }
    }

    pub fn wgts_mut(&mut self) -> &mut [Wgt] {
        &mut self.wgts
    }

    pub fn wgts(&self) -> &[Wgt] {
        &self.wgts
    }
}