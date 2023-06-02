pub fn test() {
    let v = vec![0usize, 1, 2];

    let strings = v
        .iter()
        .map(|(&x)| x.to_string() + &x.to_string())
        .collect::<Vec<String>>();

    println!("{}", strings[0]);


    let strings = v
        .iter()
        .map(|(x)| x.to_string() + &x.to_string())
        .collect::<Vec<String>>();

    println!("{}", strings[0]);
}
