/*
// return & needs lifetime
pub fn generate_str() -> Vec<&str> {
    vec!["abc", "def"];
}
 */

pub fn test() {
    // vec<&str> can be alive?
    let a = vec!["abc", "def"];
}
