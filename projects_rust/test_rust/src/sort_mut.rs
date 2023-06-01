fn sort_v(v: &mut [i32]) {
    v.sort()
}

pub fn test() {
    let mut v = vec![2, 1, 3, 5, 4];
    // ok
    //v.sort();
    //println!("{:?}", v);

    sort_v(&mut v);
    println!("{:?}", v);
}
