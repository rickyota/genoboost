
pub fn test() {

	let acc_monitor=vec![0.05f64, 0.08,0.02,-0.01,f64::NAN,f64::NAN];


    let index_max = acc_monitor
        .iter()
        .enumerate()
        .filter(|(_,a)| !a.is_nan() )
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(index, _)| index)
        .unwrap();

    println!("index_max {}",index_max);

    let index_max = acc_monitor
        .iter()
        .enumerate()
        .filter(|(_,a)| !a.is_nan() )
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(index, _)| index)
        .unwrap();


    println!("index_max {}",index_max);

}
