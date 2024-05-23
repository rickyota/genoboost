//use std::path::Path;
//use std::path::PathBuf;

// move to Dout_file
//pub fn fscore_from_fwgt(dscore: &Path, fwgt: &Path) -> PathBuf {
//    let fname_score = fwgt.file_stem().unwrap().to_str().unwrap().to_owned() + ".score";
//    //let fname_score = fwgt.file_name().unwrap().to_str().unwrap().to_owned() + ".score";
//    let fout_score = dscore.join(fname_score);
//    fout_score
//}

//pub fn fscore_concat(dscore: &Path, para: &str, fwgt: &Path) -> PathBuf {
//    // fwgt to get method name
//    // ex. 'clump_p-0.1_n-100' -> 'clump_p-0.1_n.score'
//
//    let method = fwgt
//        .file_stem()
//        .unwrap()
//        .to_str()
//        .unwrap()
//        .split(&("_".to_string() + para + "-"))
//        .collect::<Vec<&str>>()[0]
//        .to_string();
//
//    // [method]_n.score
//    //let method = fwgt
//    //    .file_name()
//    //    .unwrap()
//    //    .to_str()
//    //    .unwrap()
//    //    .split("_")
//    //    .collect::<Vec<&str>>()[0]
//    //    .to_string();
//    let fname_score = method + "_" + para + ".score";
//    let fout_score = dscore.join(fname_score);
//    fout_score
//}

// move to dout_file?
//pub fn para_from_fwgt(fwgt: &Path, para: &str) -> String {
//    let para = fwgt
//        .file_stem() // exclude .wgt here
//        .unwrap()
//        .to_str()
//        .unwrap()
//        .split(&("_".to_string() + para + "-"))
//        .collect::<Vec<&str>>()[1]
//        .to_string();
//    para
//}
//
//pub fn is_fwgt_concat(fwgt: &Path, concat_para: &str) -> bool {
//    let is_concat = fwgt
//        .file_stem()
//        .unwrap()
//        .to_str()
//        .unwrap()
//        .contains(&("_".to_string() + concat_para + "-"));
//
//    if is_concat {
//        // check if file stem name end with _(para)-*
//        if fwgt
//            .file_stem()
//            .unwrap()
//            .to_str()
//            .unwrap()
//            .split(&("_".to_string() + concat_para + "-"))
//            .collect::<Vec<&str>>()
//            .last()
//            .unwrap()
//            .contains("_")
//        {
//            panic!("File stem name should end with _(para)-*. ");
//        }
//    };
//    is_concat
//}
