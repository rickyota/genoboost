pub fn test() {
    let mut v: Vec<u8> = Vec::with_capacity(3);
    v.resize(3, 5);
    /*
    unsafe {
        v.set_len(3);
    }
    // no work without set_len
    v.fill(5);
     */
    for x in &v {
        println!("x {}", x);
    }
}
