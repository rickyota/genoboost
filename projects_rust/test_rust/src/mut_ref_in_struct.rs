pub struct Predict {
    predict_v: Vec<u8>,
    n: usize,
}

pub struct PredictSnv<'a> {
    predict_snv_v: &'a mut [u8],
    n: usize,
}

pub fn test() {
    let n = 5;
    //let mut x: Vec<u8> = Vec::with_capacity(n);
    //for _ in 0..n {
    //    x.push(0);
    //}
    let x = (0..n).map(|i| i as u8).collect::<Vec<u8>>();

    let mut pred = Predict { predict_v: x, n };

    /*
    let pred_snv = PredictSnv {
        predict_snv_v: &pred.predict_v[0..n - 1],
        n,
    };
    let y = pred_snv.predict_snv_v;
    println!("{:?}", y);
    */

    // mut is unneccesary in mut pred_snv since mut is inside
    let mut pred_snv = PredictSnv {
        predict_snv_v: &mut pred.predict_v[0..n - 1],
        n,
    };

    let mut z = pred_snv.predict_snv_v;
    println!("{:?}", z);
    (*z)[0] = 2;
    println!("{:?}", z);

    // raise error
    //println!("{:?}", y);
}
