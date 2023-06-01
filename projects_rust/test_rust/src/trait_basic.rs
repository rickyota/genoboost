// Dirty implementation for trait

// Do not implement this Trait and use AsRef instead
pub trait PredictTrait {
    fn n(&self) -> usize;
    fn val(&self) -> &[u8];
}

// Corresponds to Vec, String
pub struct PredictVec {
    predict_v: Vec<u8>,
    n: usize,
}

//impl PredictVec{
//    pub fn to_predict(&self)->&Predict{
//        Predict
//    }
//}

impl PredictTrait for PredictVec {
    fn n(&self) -> usize {
        self.n
    }
    fn val(&self) -> &[u8] {
        &self.predict_v
    }
}

// Corresponds to [], str
pub struct Predict<'a> {
    predict_snv_v: &'a [u8],
    n: usize,
}

// AsRef<str> for String is able since str is naive type
// This Predict might be impossible
// IMPOSSIBLE
/*
impl AsRef<Predict> for PredictVec {
    fn as_ref(&self) -> &Predict {
        &&Predict {
            predict_snv_v: self.val(),
            n: self.n(),
        }
    }
}
*/

// https://users.rust-lang.org/t/trait-impl-lifetime-nightmare/54735
impl<'a> PredictTrait for Predict<'a> {
    fn n(&self) -> usize {
        self.n
    }
    // Don't I need lifetime here???
    fn val(&self) -> &[u8] {
        self.predict_snv_v
    }
}

// How to separately use these two func??
// to use get_first_item_1, should impl PredictTrait for Predit and PredictSNv
// to use item_2, should impl AsRef<PredictSnv> for Predict
// ex. Twitter and Article example in rust book is item_1
// ex. fread(Path) and String is item_2
// -> these two ex. illustrate how to use!!
// -> item_2: there is one core struct and other derived struct can be used.
// -> item_1: independent traits sharing features.
// item_1: broadly used since irrelevant two types with just one Trait can be used
// Q1. For Predict and PredictSnv, which to use?
// ->  item_2 might be better!! refer to struct_asref.rs
// Q2. For SnvIndex and Snv, which to use?
//  -> item_2 is better since SnvIndex is core
// SHOULD use item_2 if possible since this is more intuitive
fn get_first_item_1<T>(pred: &T) -> u8
where
    T: PredictTrait,
{
    pred.val()[0]
}

fn get_fist_item_2<'a, T>(pred: &T) -> u8
where
    T: AsRef<Predict<'a>>,
{
    pred.as_ref().val()[0]
}

pub fn test() {
    let n = 5;
    //let mut x: Vec<u8> = Vec::with_capacity(n);
    //for ni in 0..n {
    //    x.push(ni as u8);
    //}
    let x = (0..n).map(|i| i as u8).collect::<Vec<u8>>();

    let pred = PredictVec { predict_v: x, n };

    let pred_snv = Predict {
        predict_snv_v: &pred.predict_v[0..n - 1],
        n,
    };
    let y = pred_snv.predict_snv_v;
    println!("y {:?}", y);

    println!("pred_snv.n {:?}", pred_snv.n());
    println!("pred_snv.val {:?}", pred_snv.val());

    let pred_snv_2 = Predict {
        predict_snv_v: &pred.predict_v[2..n - 1],
        n,
    };
    let z = pred_snv_2.predict_snv_v;
    println!("z {:?}", z);

    println!("pred_snv_2.n {:?}", pred_snv_2.n());
    println!("pred_snv_2.val {:?}", pred_snv_2.val());

    println!("y {:?}", y);

    println!("get_first {:?}", get_first_item_1(&pred));
    println!("get_first {:?}", get_first_item_1(&pred_snv));
    println!("get_first {:?}", get_first_item_1(&pred_snv_2));
}
