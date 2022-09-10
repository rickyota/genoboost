use std::io::BufWriter;
use std::path::Path;

use crate::wgt_boost;
use crate::wgt_boost::io;
use crate::BoostType;
use crate::WgtBoost;

#[derive(Debug, Clone)]
pub struct WgtBoosts {
    wgts: Vec<WgtBoost>,
    // coef type for whole wgt, which determines column
    //wgt_column: WgtColumn,
    // for column
    boost_type: BoostType,
    columns: Vec<String>,
    index_next_write: usize,
}

impl WgtBoosts {
    pub fn new(boost_type: BoostType) -> Self {
        WgtBoosts {
            wgts: vec![],
            boost_type,
            columns: io::wgt_columns(boost_type),
            index_next_write: 0,
        }
    }

    pub fn new_from_file_pathbuf(fin_wgt: &Path, boost_type: BoostType) -> Self {
        let wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(fin_wgt, boost_type);
        let wgts_len = wgts.len();

        // TODO: check colums match io::wgt_columns()

        WgtBoosts {
            wgts,
            boost_type,
            columns: io::wgt_columns(boost_type),
            index_next_write: wgts_len,
        }
    }

    // TODO: estimate boost_type from column
    pub fn new_from_file(fwgt: &Path, boost_type: BoostType) -> Self {
        let wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts_file(fwgt, boost_type);
        let wgts_len = wgts.len();

        // TODO: check colums match io::wgt_columns()

        WgtBoosts {
            wgts,
            boost_type,
            columns: io::wgt_columns(boost_type),
            index_next_write: wgts_len,
        }
    }

    // TODO: estimate boost_type from column
    pub fn new_from_file_dir(dout: &Path, boost_type: BoostType) -> Self {
        let wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(dout, boost_type);
        let wgts_len = wgts.len();

        // TODO: check colums match io::wgt_columns()

        WgtBoosts {
            wgts,
            boost_type,
            columns: io::wgt_columns(boost_type),
            index_next_write: wgts_len,
        }
    }

    pub fn boost_type(&self) -> BoostType {
        self.boost_type
    }

    pub fn wgts(&self) -> &[WgtBoost] {
        &self.wgts
    }
    pub fn wgts_mut(&mut self) -> &mut [WgtBoost] {
        &mut self.wgts
    }
    pub fn index_next_write(&self) -> usize {
        self.index_next_write
    }
    pub fn columns(&self) -> &[String] {
        &self.columns
    }

    pub fn add_wgt(&mut self, wgt: WgtBoost) {
        assert_eq!(self.wgts().len(), wgt.iteration());
        self.wgts.push(wgt);
    }

    pub fn write_wgt_writer<W: std::io::Write>(&mut self, writer: &mut BufWriter<W>) {
        // write wgt from iter_next_writer to newest
        let len = self.wgts().len();
        for index_w in self.index_next_write()..len {
            // cannot write...
            // -> refleshed after written
            //writer.write("aaa".as_bytes()).unwrap();
            println!("writing index {}", index_w);
            let wgt = &self.wgts()[index_w];
            wgt.write_wgt_writer(writer, self.columns());
            /*             let str = wgt.string_write();
            println!("str wgt {}", str);
            writer.write(str.as_bytes()).unwrap(); */
        }
        self.index_next_write = len;
    }
}
