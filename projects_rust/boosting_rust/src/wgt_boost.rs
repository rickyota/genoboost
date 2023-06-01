// remove "get" from "get_model()"
// where to implement functions?
// 1. wgt.func()
// 2. wgt.get_wgt_kind().func()
// 3. match get_wgt_kind()...  ~.func()

pub mod io;

use std::collections::HashMap;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::ContingencyTable;
use genetics::wgt::{Coef, Model, Wgt, WgtKind, WgtTrait};
use genetics::CovKind;

#[derive(Debug, Clone)]
pub struct WgtBoost {
    wgt: Wgt,
    iteration: usize,
    loss: Option<f64>,
    contingency_table: Option<ContingencyTable>,
    is_eps: Option<bool>,
    is_eff_eps: Option<bool>,
}

impl WgtTrait for WgtBoost {
    fn wgt(&self) -> &Wgt {
        self.wgt()
    }
    fn wgt_mut(&mut self) -> &mut Wgt {
        self.wgt_mut()
    }

    fn set_snv_index(&mut self, mi: Option<usize>) {
        let kind = self.wgt_mut().kind_mut();
        // TODO: move to WgtKind::
        if let WgtKind::Snv(_, _, ref mut index) = kind {
            *index = mi;
        }
    }
}

impl WgtBoost {
    pub fn new(
        wgt: Wgt,
        iteration: usize,
        loss: Option<f64>,
        contingency_table: Option<ContingencyTable>,
        is_eps: Option<bool>,
        is_eff_eps: Option<bool>,
    ) -> WgtBoost {
        let wgt_boost = WgtBoost {
            wgt,
            iteration,
            loss,
            contingency_table,
            is_eps,
            is_eff_eps,
        };
        //wgt_boost.check_sum_abcd();
        wgt_boost
    }

    pub fn construct_wgt(
        wgt: Wgt,
        iteration: usize,
        loss: f64,
        contingency_table: ContingencyTable,
        is_eps: bool,
        is_eff_eps: Option<bool>,
    ) -> WgtBoost {
        let wgt_boost = WgtBoost {
            wgt,
            iteration,
            loss: Some(loss),
            contingency_table: Some(contingency_table),
            is_eps: Some(is_eps),
            is_eff_eps,
        };
        //wgt_boost.check_sum_abcd();
        wgt_boost
    }
    pub fn construct_wgt_loss(wgt: Wgt, iteration: usize, loss: Option<f64>) -> WgtBoost {
        let wgt_boost = WgtBoost {
            wgt,
            iteration,
            loss,
            contingency_table: None,
            //contingency_table: Some(ContingencyTable::new_four(contingency_table.unwrap())), //abcd_sum,
            is_eps: None,
            is_eff_eps: None,
        };
        //wgt_boost.check_sum_abcd();
        wgt_boost
    }

    pub fn construct_wgt_iteration(wgt: Wgt, iteration: usize) -> WgtBoost {
        let wgt_boost = WgtBoost {
            wgt,
            iteration,
            loss: None,
            contingency_table: None,
            is_eps: None,
            is_eff_eps: None,
        };
        wgt_boost
    }

    /*
    /// check sum ~1.0
    fn check_sum_abcd(&self) {
        //return ();
        // temporary not use
        let (a, b, c, d) = self.abcd_sum();
        if ((a + b + c + d) - 1.0).abs() > 1e-10 {
            panic!(
                "Sum of abcd is not 1.0, diff is: {}",
                ((a + b + c + d) - 1.0).abs()
            )
        }
    }
     */

    pub fn wgt(&self) -> &Wgt {
        &self.wgt
    }

    pub fn wgt_mut(&mut self) -> &mut Wgt {
        &mut self.wgt
    }

    // FIXME: Opion<f64>
    pub fn loss(&self) -> f64 {
        self.loss.unwrap()
    }

    pub fn loss_option(&self) -> Option<f64> {
        self.loss
    }

    pub fn set_coef(&mut self, coef: Coef) {
        self.wgt_mut().model_mut().set_coef(coef);
    }

    pub fn set_contingency_table(&mut self, contingency_table: ContingencyTable, is_eps: bool) {
        self.contingency_table = Some(contingency_table);
        self.is_eps = Some(is_eps);
    }

    pub fn set_is_eps(&mut self, is_eps: bool) {
        self.is_eps = Some(is_eps);
    }

    pub fn set_is_eff_eps(&mut self, is_eff_eps: bool) {
        self.is_eff_eps = Some(is_eff_eps);
    }

    pub fn contingency_table(&self) -> Option<ContingencyTable> {
        self.contingency_table
    }

    pub fn table4_sum(&self) -> (f64, f64, f64, f64) {
        self.contingency_table().unwrap().four()
        /*
        let table = self.contingency_table.unwrap();
        match table {
            ContingencyTable::Four((x0, x1, x2, x3)) => (x0, x1, x2, x3),
            _ => panic!("Not Four."),
        }
         */
    }

    pub fn table2_sum(&self) -> (f64, f64) {
        self.contingency_table().unwrap().two()
        /*
        let table = self.contingency_table.unwrap();
        match table {
            ContingencyTable::Two((x0, x1)) => (x0, x1),
            _ => panic!("Not Two."),
        }
         */
    }

    pub fn table7_sum(&self) -> (f64, f64, f64, f64, f64, f64, f64) {
        self.contingency_table().unwrap().seven()
        /*
        let table = self.contingency_table.unwrap();
        print!("table {:?}", table);
        match table {
            ContingencyTable::Seven((x0, x1, x2, x3, x4, x5, x6)) => (x0, x1, x2, x3, x4, x5, x6),
            _ => panic!("Not Seven."),
        }
         */
    }

    pub fn iteration(&self) -> usize {
        self.iteration
    }

    fn iteration_string(&self) -> String {
        self.iteration.to_string()
    }

    fn is_eps_string(&self) -> String {
        //self.is_eps.unwrap().to_string()
        match self.is_eps {
            Some(x) => x.to_string(),
            None => "none".to_string(),
        }
    }
    fn is_eff_eps_string(&self) -> String {
        //self.is_eps.unwrap().to_string()
        match self.is_eff_eps {
            Some(x) => x.to_string(),
            None => "none".to_string(),
        }
    }

    fn content_hash(&self) -> HashMap<String, String> {
        //let content_hash: HashMap<String, String> = match self.wgt().kind() {
        match self.wgt().kind() {
            WgtKind::Snv(snv, _, _) => {
                //let snv = snv_wgt.snv_index();
                let model = self.wgt().model();

                let mut hash = HashMap::from([
                    ("iteration".to_owned(), self.iteration_string()),
                    ("kind".to_owned(), "SNV".into()),
                    ("var".to_owned(), snv.rs().into()),
                    ("model".to_owned(), model.model_name().into()),
                    ("threshold".to_owned(), model.threshold_string().into()),
                    ("eps".to_owned(), self.is_eps_string()),
                    ("eff_eps".to_owned(), self.is_eff_eps_string()),
                    //("alpha".to_owned(), model.alpha_string().into()),
                    //("const".to_owned(), model.const_string().into()),
                    ("chrom".to_owned(), snv.chrom().to_string().into()),
                    ("pos".to_owned(), snv.pos().to_string().into()),
                    ("a1".to_owned(), snv.a1().into()),
                    ("a2".to_owned(), snv.a2().into()),
                    //("nan".to_owned(), "NaN".into()),
                ]);

                let hash_snv = model.to_string_hash();
                log::debug!("hash_snv {:?}", hash_snv);
                hash.extend(hash_snv);

                hash
            }

            WgtKind::Cov(cov_id) => {
                let kind_name = match cov_id.kind() {
                    CovKind::Const => "CONST",
                    //VarKind::Var => "VAR",
                    CovKind::Cov => "COV",
                };
                let model = self.wgt().model();

                // make "NaN" as const string
                let mut hash = HashMap::from([
                    ("iteration".to_owned(), self.iteration_string()),
                    ("kind".to_owned(), kind_name.into()),
                    ("var".to_owned(), cov_id.name().into()),
                    ("model".to_owned(), model.model_name().into()),
                    ("threshold".to_owned(), model.threshold_string().into()),
                    //("alpha".to_owned(), model.alpha_string().into()),
                    //("const".to_owned(), model.const_string().into()),
                    ("eps".to_owned(), "NaN".into()),
                    ("eff_eps".to_owned(), "NaN".into()),
                    ("chrom".to_owned(), "NaN".into()),
                    ("pos".to_owned(), "NaN".into()),
                    ("a1".to_owned(), "N".into()),
                    ("a2".to_owned(), "N".into()),
                    //("nan".to_owned(), "NaN".into()),
                ]);

                let hash_snv = model.to_string_hash();
                hash.extend(hash_snv);
                hash

                /*
                //let cov_id = cov_wgt.var();
                match self.wgt().model().coef().model_type() {
                    ModelType::Linear => {
                        let kind_name = match cov_id.kind() {
                            CovKind::Const => "CONST",
                            //VarKind::Var => "VAR",
                            CovKind::Cov => "COV",
                        };
                        let model = self.wgt().model();

                        // make "NaN" as const string
                        let mut hash = HashMap::from([
                            ("iteration".to_owned(), self.iteration_string()),
                            ("kind".to_owned(), kind_name.into()),
                            ("var".to_owned(), cov_id.name().into()),
                            ("model".to_owned(), model.model_name().into()),
                            ("threshold".to_owned(), model.threshold_string().into()),
                            //("alpha".to_owned(), model.alpha_string().into()),
                            //("const".to_owned(), model.const_string().into()),
                            ("eps".to_owned(), "NaN".into()),
                            // eff_eps?
                            ("eps_eps".to_owned(), "NaN".into()),
                            ("chrom".to_owned(), "NaN".into()),
                            ("pos".to_owned(), "NaN".into()),
                            ("a1".to_owned(), "N".into()),
                            ("a2".to_owned(), "N".into()),
                            //("nan".to_owned(), "NaN".into()),
                        ]);

                        let hash_snv = model.to_string_hash();
                        hash.extend(hash_snv);
                        hash
                    }
                    ModelType::Binary => {
                        let kind_name = match cov_id.kind() {
                            CovKind::Const => "CONST",
                            //VarKind::Var => "VAR",
                            CovKind::Cov => "COV",
                        };
                        let model = self.wgt().model();

                        // make "NaN" as const string
                        let mut hash = HashMap::from([
                            ("iteration".to_owned(), self.iteration_string()),
                            ("kind".to_owned(), kind_name.into()),
                            ("var".to_owned(), cov_id.name().into()),
                            ("model".to_owned(), model.model_name().into()),
                            ("threshold".to_owned(), model.threshold_string().into()),
                            //("alpha".to_owned(), model.alpha_string().into()),
                            //("const".to_owned(), model.const_string().into()),
                            ("eps".to_owned(), "NaN".into()),
                            ("eps_eps".to_owned(), "NaN".into()),
                            ("chrom".to_owned(), "NaN".into()),
                            ("pos".to_owned(), "NaN".into()),
                            ("a1".to_owned(), "N".into()),
                            ("a2".to_owned(), "N".into()),
                            //("nan".to_owned(), "NaN".into()),
                        ]);

                        let hash_snv = model.to_string_hash();
                        hash.extend(hash_snv);
                        hash
                    }
                    _ => unimplemented!(),
                    */
            }
        }
    }

    // TODO: move to boost_wgts and make columns the same
    // distinguish const by var.get_kind()
    pub fn string_write(&self, columns: &[String]) -> String {
        //let boost_type = BoostType::FreeModelMissing;

        let content_hash = self.content_hash();

        //let boost_type = BoostType::Ada;
        // only indicates columns if col and col in content are different
        //let content_to_col = io::content_to_col(self.wgt().kind(), boost_type);
        //match self.kind() {

        //let columns = io::wgt_columns(boost_type);
        log::debug!("columns {:?}", &columns);
        let strings: Vec<String> = columns
            .iter()
            .map(|content_col| {
                content_hash
                    .get(content_col)
                    .unwrap_or(&"NaN".to_owned())
                    .into()
            })
            .collect::<Vec<String>>();
        /*
        let strings: Vec<String> = columns
            .iter()
            .map(|col| {
                content_to_col.get(col).unwrap_or(col)
                //.into()
            })
            .map(|content_col| content_hash.get(content_col).unwrap().into())
            .collect::<Vec<String>>();
             */
        strings.join("\t") + "\n"
    }

    /*     // TODO: move to boost_wgts and make columns the same
       // distinguish const by var.get_kind()
       pub fn string_write(&self) -> String {
           let boost_type = BoostType::FreeModelMissing;
           //let boost_type = BoostType::Ada;
           // only indicates columns if col and col in content are different
           //let content_to_col = io::content_to_col(self.wgt().kind(), boost_type);
           //match self.kind() {
           let content_hash: HashMap<String, String> = match self.wgt().kind() {
               WgtKind::Snv(snv, _, _) => {
                   //let snv = snv_wgt.snv_index();
                   let model = self.wgt().model();

                   let mut hash = HashMap::from([
                       ("iteration".to_owned(), self.iteration_string()),
                       ("kind".to_owned(), "SNV".into()),
                       ("var".to_owned(), snv.rs().into()),
                       ("model".to_owned(), model.model_name().into()),
                       ("threshold".to_owned(), model.threshold_string().into()),
                       //("alpha".to_owned(), model.alpha_string().into()),
                       //("const".to_owned(), model.const_string().into()),
                       ("chrom".to_owned(), snv.chrom().to_string().into()),
                       ("pos".to_owned(), snv.pos().to_string().into()),
                       ("a1".to_owned(), snv.a1().into()),
                       ("a2".to_owned(), snv.a2().into()),
                       //("nan".to_owned(), "NaN".into()),
                   ]);

                   let hash_snv = model.to_string_hash();
                   log::debug!("hash_snv {:?}", hash_snv);
                   hash.extend(hash_snv);

                   hash
               }

               WgtKind::Cov(cov_id) => {
                   //let cov_id = cov_wgt.var();
                   match self.wgt().model().coef().model_type() {
                       ModelType::Linear => {
                           let kind_name = match cov_id.kind() {
                               CovKind::Const => "CONST",
                               //VarKind::Var => "VAR",
                               CovKind::Cov => "COV",
                           };
                           let model = self.wgt().model();

                           // make "NaN" as const string
                           let mut hash = HashMap::from([
                               ("iteration".to_owned(), self.iteration_string()),
                               ("kind".to_owned(), kind_name.into()),
                               ("var".to_owned(), cov_id.name().into()),
                               ("model".to_owned(), model.model_name().into()),
                               ("threshold".to_owned(), model.threshold_string().into()),
                               //("alpha".to_owned(), model.alpha_string().into()),
                               //("const".to_owned(), model.const_string().into()),
                               ("chrom".to_owned(), "NaN".into()),
                               ("pos".to_owned(), "NaN".into()),
                               ("a1".to_owned(), "N".into()),
                               ("a2".to_owned(), "N".into()),
                               //("nan".to_owned(), "NaN".into()),
                           ]);

                           let hash_snv = model.to_string_hash();
                           hash.extend(hash_snv);
                           hash
                       }
                       ModelType::Binary => {
                           let kind_name = match cov_id.kind() {
                               CovKind::Const => "CONST",
                               //VarKind::Var => "VAR",
                               CovKind::Cov => "COV",
                           };
                           let model = self.wgt().model();

                           // make "NaN" as const string
                           let mut hash = HashMap::from([
                               ("iteration".to_owned(), self.iteration_string()),
                               ("kind".to_owned(), kind_name.into()),
                               ("var".to_owned(), cov_id.name().into()),
                               ("model".to_owned(), model.model_name().into()),
                               ("threshold".to_owned(), model.threshold_string().into()),
                               //("alpha".to_owned(), model.alpha_string().into()),
                               //("const".to_owned(), model.const_string().into()),
                               ("chrom".to_owned(), "NaN".into()),
                               ("pos".to_owned(), "NaN".into()),
                               ("a1".to_owned(), "N".into()),
                               ("a2".to_owned(), "N".into()),
                               //("nan".to_owned(), "NaN".into()),
                           ]);

                           let hash_snv = model.to_string_hash();
                           hash.extend(hash_snv);
                           hash
                       }
                       _ => unimplemented!(),
                   }
               }
           };
           let columns = io::wgt_columns(boost_type);
           log::debug!("columns {:?}", &columns);
           let strings: Vec<String> = columns
               .iter()
               .map(|content_col| {
                   content_hash
                       .get(content_col)
                       .unwrap_or(&"NaN".to_owned())
                       .into()
               })
               .collect::<Vec<String>>();
           /*
           let strings: Vec<String> = columns
               .iter()
               .map(|col| {
                   content_to_col.get(col).unwrap_or(col)
                   //.into()
               })
               .map(|content_col| content_hash.get(content_col).unwrap().into())
               .collect::<Vec<String>>();
                */
           strings.join("\t") + "\n"
       }
    */

    /// should write CONST not const
    pub fn write_wgt(&self, dout: &Path, columns: &[String]) {
        let fwgt = io::get_fname_wgt(dout);
        let file = OpenOptions::new().append(true).open(fwgt).unwrap();
        //let mut writer = BufWriter::new(File::create(fout).unwrap());
        let mut writer = BufWriter::new(file);
        let string = self.string_write(columns);
        //let str = "a\tb\n".to_owned();
        writer.write(string.as_bytes()).unwrap();
    }

    pub fn write_wgt_writer<W: std::io::Write>(
        &self,
        writer: &mut BufWriter<W>,
        columns: &[String],
    ) {
        /*
        let fwgt = io::get_fname_wgt(fout);
        let file = OpenOptions::new().append(true).open(fwgt).unwrap();
        //let mut writer = BufWriter::new(File::create(fout).unwrap());
        let mut writer = BufWriter::new(file);
         */
        let string = self.string_write(columns);
        log::debug!("str wgt {}", &string);
        //let str = "a\tb\n".to_owned();
        writer.write(string.as_bytes()).unwrap();
        // capture error
        writer.flush().unwrap();
    }
}
