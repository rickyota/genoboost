// add ./common/
// https://doc.rust-lang.org/rust-by-example/testing/integration_testing.html

//use boosting_rust;

//use std::env;
//use std::path::Path;

/*
#[test]
fn test_run_boosting_score() {
    let path = env::current_dir().unwrap();
    println!("current dir: {}", path.display());

    let fout = String::from("../../test/result/boosting/toy1/unittest/toy1_regcov");
    let fin_wgt = String::from("../../test/data/toy1/genot_for_score");
    //let fin_wgt = String::from("../../test/boosting/result/toy1/unittest/toy1_for_score");

    // tmp: better to put the file in this crate
    let fin = String::from("../../test/data/toy1/genot");
    let fin_cov = Some(String::from("../../test/data/toy1/genot.cov"));

    let iters = [1, 2, 10];

    boosting_rust::run_boosting_score(&fout, &fin, &iters, &fin_wgt, fin_cov.as_deref(), None);

    println!("Done!!");
}

#[test]
fn test_run_boosting_score_withregcov() {
    let path = env::current_dir().unwrap();
    println!("current dir: {}", path.display());

    let fout = String::from("../../test/result/boosting/toy1/unittest/toy1");
    let fin_wgt = String::from("../../test/data/toy1/genot_regcov_for_score");
    //let fin_wgt = String::from("../../test/boosting/result/toy1/unittest/toy1_for_score");

    // tmp: better to put the file in this crate
    let fin = String::from("../../test/data/toy1/genot");
    //let fin = String::from("./data/boosting/1000g/1kg_phase1_all/1kg_phase1_all");

    let iters = [1, 2, 10];

    boosting_rust::run_boosting_score(&fout, &fin, &iters, &fin_wgt, None, None);

    println!("Done!!");
}

#[test]
fn test_run_boosting_score_toy2_withregcov() {
    let path = env::current_dir().unwrap();
    println!("current dir: {}", path.display());

    let fout = String::from("../../test/result/boosting/toy2/unittest/toy2");
    //let fin_wgt = String::from("../../test/boosting/result/toy2/unittest/toy2_regcov_for_score");
    let fin_wgt = String::from("../../test/data/toy2/genot_regcov_for_score");

    // tmp: better to put the file in this crate
    let fin = String::from("../../test/data/toy2/genot.chrom%");
    let fin_cov = Some(String::from("../../test/data/toy2/genot.cov"));
    let fin_sample = Some(String::from(
        "../../test/data/toy2/genot.cv0.samples",
    ));

    let iters = [1, 2, 10];

    boosting_rust::run_boosting_score(
        &fout,
        &fin,
        &iters,
        &fin_wgt,
        fin_cov.as_deref(),
        fin_sample.as_deref(),
    );

    println!("Done!!");
}

*/
