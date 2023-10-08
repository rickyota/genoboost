use super::{wgt, Wgt};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Wgts {
    wgts: Vec<Wgt>,
    //columns: Vec<String>,
    //index_next_write: usize,
}

impl Wgts {
    //pub fn new_from_file(fwgt: &Path,    use_snv_pos: bool,is_nonadd:bool) -> Self {
    pub fn new_from_file(fwgt: &Path, is_nonadd: bool) -> Self {
        let wgts: Vec<Wgt> = wgt::io::load_wgts_file(fwgt, is_nonadd);
        //let wgts: Vec<Wgt> = wgt::io::load_wgts_file(fwgt,use_snv_pos,is_nonadd);

        // TODO: check colums match io::wgt_columns()

        Wgts { wgts }
    }

    pub fn wgts_mut(&mut self) -> &mut [Wgt] {
        &mut self.wgts
    }

    pub fn wgts(&self) -> &[Wgt] {
        &self.wgts
    }
}
