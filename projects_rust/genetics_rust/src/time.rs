// move to other crate

// `let (_, elapsed_time) = elapsed!(function_we_want_to_eval());`
macro_rules! elapsed {
    ($a:expr) => {{
        let start = std::time::Instant::now();
        let return_value = $a;
        let end = std::time::Instant::now();
        (return_value, (end - start))
    }};
}
