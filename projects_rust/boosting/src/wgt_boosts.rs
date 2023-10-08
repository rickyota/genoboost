use std::collections::HashSet;
use std::io::BufWriter;
use std::path::Path;

use crate::wgt_boost;
use crate::wgt_boost::io;
use crate::BoostType;
use crate::WgtBoost;
use genetics::textfile;

#[derive(Debug, Clone)]
pub struct WgtBoosts {
    wgts: Vec<WgtBoost>,
    // coef type for whole wgt, which determines column
    //wgt_column: WgtColumn,
    // for column
    // None when loading in score
    boost_type: Option<BoostType>,
    columns: Vec<String>,
    index_next_write: usize,
}

impl WgtBoosts {
    pub fn new(boost_type: BoostType) -> Self {
        WgtBoosts {
            wgts: vec![],
            boost_type: Some(boost_type),
            columns: io::wgt_columns(boost_type),
            index_next_write: 0,
        }
    }

    //pub fn new_from_file_pathbuf(fin_wgt: &Path, boost_type: BoostType) -> Self {
    //    let wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(fin_wgt, boost_type);
    //    let wgts_len = wgts.len();

    //    // TODO: check colums match io::wgt_columns()

    //    WgtBoosts {
    //        wgts,
    //        boost_type,
    //        columns: io::wgt_columns(boost_type),
    //        index_next_write: wgts_len,
    //    }
    //}

    // TODO: estimate boost_type from column
    /// call on validate
    /// boost_type is not determined (since only info is columns)...
    //pub fn new_from_file(fwgt: &Path, boost_type: BoostType) -> Self {
    pub fn new_from_file(fwgt: &Path) -> Self {
        let wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts_file(fwgt);
        //let wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts_file(fwgt, boost_type);
        let wgts_len = wgts.len();

        // TODO: check colums match io::wgt_columns()

        WgtBoosts {
            wgts,
            boost_type: None,
            //columns: io::wgt_columns(boost_type),
            columns: textfile::load_table_header(fwgt),
            index_next_write: wgts_len,
        }
    }

    // TODO: estimate boost_type from column
    /// call on training
    /// boost_type is required
    pub fn new_from_file_dir(dout: &Path, boost_type: BoostType) -> Self {
        let wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(dout);
        //let wgts: Vec<WgtBoost> = wgt_boost::io::load_wgts(dout, boost_type);
        let wgts_len = wgts.len();

        // TODO: check colums match io::wgt_columns()

        let fwgt = wgt_boost::io::get_fname_wgt(dout);
        if io::wgt_columns(boost_type) != textfile::load_table_header(&fwgt) {
            panic!(
                "Existing wgtboost does not match boost_type: {:?}",
                boost_type
            );
        }

        WgtBoosts {
            wgts,
            boost_type: Some(boost_type),
            //columns: io::wgt_columns(boost_type),
            columns: textfile::load_table_header(&fwgt),
            index_next_write: wgts_len,
        }
    }

    pub fn boost_type(&self) -> BoostType {
        self.boost_type.unwrap()
        //self.boost_type
    }

    pub fn wgts(&self) -> &[WgtBoost] {
        &self.wgts
    }
    pub fn wgts_mut(&mut self) -> &mut [WgtBoost] {
        &mut self.wgts
    }

    pub fn wgts_n(&self) -> usize {
        self.wgts().len()
    }

    pub fn index_next_write(&self) -> usize {
        self.index_next_write
    }
    pub fn columns(&self) -> &[String] {
        &self.columns
    }

    //pub fn use_missing(&self) -> bool {
    //    // TODO: ok?
    //    self.columns.contains(&"scorem".to_string())
    //}

    pub fn fill_missing(&self) -> bool {
        // TODO: ok?
        !self.columns.contains(&"scorem".to_string())
    }

    // TODO: test
    pub fn count_unique_snv(&self) -> usize {
        // https://www.reddit.com/r/rust/comments/b4cxrj/how_to_count_number_of_unique_items_in_an_array/
        self.wgts()
            .iter()
            .filter(|x| x.wgt().kind().is_snv())
            .map(|x| x.wgt().kind().snv_index().rs().to_string())
            .collect::<HashSet<String>>()
            .len()
    }

    pub fn last_wgt(&self) -> Option<&WgtBoost> {
        self.wgts().last()
        //if self.wgts_n()==0{
        //    return None;
        //}
        //Some(self.wgts())
    }

    /// check if last n items are all snv
    /// None if items are less than n
    pub fn is_last_snv(&self, n: usize) -> Option<bool> {
        if self.wgts_n() < n {
            return None;
        }

        Some(
            self.wgts()
                .iter()
                .rev()
                .take(n)
                .filter(|x| x.wgt().kind().is_snv())
                .count()
                == n,
        )
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
            // -> refreshed after written
            //writer.write("aaa".as_bytes()).unwrap();
            log::debug!("writing index {}", index_w);
            let wgt = &self.wgts()[index_w];
            wgt.write_wgt_writer(writer, self.columns());
            /*             let str = wgt.string_write();
            log::debug!("str wgt {}", str);
            writer.write(str.as_bytes()).unwrap(); */
        }
        self.index_next_write = len;
    }
}
