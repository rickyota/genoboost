use std::path::{Path, PathBuf};

//use boosting_rust::predict;
use boosting_rust::boosting::loss::LossStruct;
use genetics::{alloc, Dataset};

fn setup_vars(
    fin: &Path,
    fin_snv: Option<&Path>,
    fin_sample: Option<&Path>,
) -> (Dataset, Vec<f64>, LossStruct) {
    let dataset: Dataset = Dataset::new(fin, fin_snv, fin_sample, None, false);
    let n = dataset.samples().samples_n();
    let m = dataset.snvs().snvs_n();

    let mut ps: Vec<f64> = alloc::with_capacity_align_f64(n + 32);
    ps.resize(n, 1.0f64 / (n as f64));
    ps.resize(32, 0.0);

    let losss = LossStruct::new(boosting_rust::BoostType::ConstAda, m);

    (dataset, ps, losss)
}

pub fn setup_test_n100k() -> (Dataset, Vec<f64>, LossStruct) {
    let fin = PathBuf::from("../../test/data/1kg_n100000/1kg_n100000");
    //let fin = String::from("../../test/data/1kg_n100000/genot");
    // too large
    //let fin = String::from("../../test/boosting/data/1kg_n1000000/1kg_n1000000");
    //let fin = String::from("../../test/boosting/data/toy1/toy1");
    //let fin_snv = String::from("../../data/boosting/1000g/1kg_phase1_all/use_snvs_n100000.snvs");
    //let fin_snv = Some(fin_snv.as_str());
    // all snvs are too much
    let fin_snv = None;
    let fin_sample = None;

    setup_vars(&fin, fin_snv, fin_sample)
}

pub fn setup_test_n10k() -> (Dataset, Vec<f64>, LossStruct) {
    let fin = PathBuf::from("../../test/data/1kg_n10000/1kg_n10000");
    //let fin = String::from("../../test/boosting/data/1kg_n10000/genot");
    // too large
    //let fin = String::from("../../test/boosting/data/1kg_n1000000/1kg_n1000000");
    let fin_snv = None;
    let fin_sample = None;

    setup_vars(&fin, fin_snv, fin_sample)
}
