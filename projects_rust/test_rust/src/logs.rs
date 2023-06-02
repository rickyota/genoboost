//use log::{error, info, warn};
// THIS IS NOT NECESSARY!!!!
// since in Cargo.toml?
//use log;
// use env_logger;

pub fn test() {
    //let verbose = true;
    let verbose = false;
    println!("in log test");

    // only way to change level is set env
    if verbose {
        std::env::set_var("RUST_LOG", "info");
    } else {
        std::env::set_var("RUST_LOG", "warn");
    }

    // in main
    env_logger::init();

    log::error!("error");
    log::warn!("warn");
    log::info!("info");
}
