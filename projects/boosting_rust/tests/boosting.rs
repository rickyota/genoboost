// add ./common/
// https://doc.rust-lang.org/rust-by-example/testing/integration_testing.html

//use boosting_rust::{self, BoostMethod, BoostType};

//use std::env;
//use std::path::Path;

// commented out since do not want to create unittest dir from tests
/*
#[test]
fn test_run_boosting() {
    let path = env::current_dir().unwrap();
    println!("current dir: {}", path.display());

    let fout = String::from("../../test/result/boosting/toy1/unittest/toy1");
    //let fout = String::from("../../tests/boosting/result/toy1/tests_toy1");
    use io_rust::text;
    use std::fs;
    use std::path::Path;
    let fwgt = fout.to_owned() + ".wgt";
    if text::exist_file(&fwgt) {
        // delete file
        fs::remove_file(fwgt).unwrap_or_else(|why| {
            println!("! {:?}", why.kind());
        });
    }

    // tmp: better to put the file in this crate
    let fin = String::from("../../test/result/boosting/toy1/genot");
    //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

    //use std::path::Path;
    let fin_bim: String = path.to_str().unwrap().to_owned() + "/" + &fin + ".bim"; // + ".bim";
    println!("fin_bim {:?}", fin_bim);
    let path = Path::new(&fin_bim);
    println!("cano {:?}", path.canonicalize().unwrap());

    let t = 10;
    let boost_type = BoostType::Logit;

    boosting_rust::run_boosting(
        &fout,
        &fin,
        t,
        &boost_type,
        None,
        None,
        None,
        None,
        false,
        BoostMethod::Classic,
        None,
    );

    println!("Done!!");
}
*/

/*
#[test]
fn test_run_boosting_withregcov() {
    let path = env::current_dir().unwrap();
    println!("current dir: {}", path.display());

    let fout = String::from("../../test/result/boosting/toy1/unittest/toy1_regcov");
    //let fout = String::from("../../tests/boosting/result/toy1/tests_regcov_toy1");
    use io_rust::text;
    use std::fs;
    let fwgt = fout.to_owned() + ".wgt";
    if text::exist_file(&fwgt) {
        // delete file
        fs::remove_file(fwgt).unwrap_or_else(|why| {
            println!("! {:?}", why.kind());
        });
    }

    // tmp: better to put the file in this crate
    let fin = String::from("../../test/data/toy1/genot");
    //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

    let fin_wgt_cov = Some(String::from(
        "../../test/data/toy1/genot.regcov.wgtcov",
    ));

    let fin_cov = Some(String::from("../../test/data/toy1/genot.cov"));

    println!("fin_wgt_cov {:?}", fin_wgt_cov);

    let t = 10;
    let boost_type = BoostType::Logit;

    boosting_rust::run_boosting(
        &fout,
        &fin,
        t,
        &boost_type,
        None,
        None,
        fin_cov.as_deref(),
        fin_wgt_cov.as_deref(),
        false,
        BoostMethod::Classic,
        None,
    );

    println!("Done!!");
}
*/

/*
#[test]
fn test_run_boosting_toy2_withregcov() {
    let path = env::current_dir().unwrap();
    println!("current dir: {}", path.display());

    let fout = String::from("../../test/result/boosting/toy2/unittest/genot_regcov");
    //let fout = String::from("../../tests/boosting/result/toy2/tests_regcov_toy2");
    use io_rust::text;
    use std::fs;
    let fwgt = fout.to_owned() + ".wgt";
    if text::exist_file(&fwgt) {
        // delete file
        fs::remove_file(fwgt).unwrap_or_else(|why| {
            println!("! {:?}", why.kind());
        });
    }

    // tmp: better to put the file in this crate
    let fin = String::from("../../test/data/toy2/genot.chrom%");
    //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

    /*
    //use std::path::Path;
    let fin_bim: String = path.to_str().unwrap().to_owned() + "/" + &fin + ".bim"; // + ".bim";
    println!("fin_bim {:?}", fin_bim);
    let path = Path::new(&fin_bim);
    println!("cano {:?}", path.canonicalize().unwrap());
    */

    let fin_wgt_cov = Some(String::from(
        "../../test/data/toy2/genot.regcov.wgtcov",
    ));

    let fin_cov = Some(String::from("../../test/data/toy2/genot.cov"));
    let fin_snv = Some(String::from("../../test/data/toy2/genot.snvs"));
    let fin_sample = Some(String::from(
        "../../test/data/toy2/genot.cv0.samples",
    ));

    println!("fin_wgt_cov {:?}", fin_wgt_cov);

    let t = 10;
    let boost_type = BoostType::Logit;

    boosting_rust::run_boosting(
        &fout,
        &fin,
        t,
        &boost_type,
        fin_snv.as_deref(),
        fin_sample.as_deref(),
        fin_cov.as_deref(),
        fin_wgt_cov.as_deref(),
        false,
        BoostMethod::Classic,
        None,
    );

    println!("Done!!");
}
*/
