//! implement Predict{Vec<u8>} and Predict{&[u8]} simultaneously
//! There are two ways
//! Method 2 is more scalable
//! 1. Predict<T>
//! 2. use trait

/*
Vec in struct should not be &[u8]

struct Predict<T: >>{
    inner: T
}

// This is dataset, so once created, never be mut
struct Dataset{
    predict: Option<Predict<Vec<u8>>>,
    y: Option<Y>,
}

// when using, &Predict is fine, should not use Predict<&[u8]>
// ok
take(&predict)
// ng
take(predict.to_borrow())

*/

/*
// Do not implement this Trait and use AsRef instead
pub trait PredictTrait {
    fn val(&self) -> &[u8];
    fn n(&self) -> usize;
}
 */

// 1. One way is make T itself to asref
pub struct Predict<T: AsRef<[u8]>> {
    inner: T,
}

impl Predict<Vec<u8>> {
    // Self=Vec<u8>
    fn new(x: Vec<u8>) -> Self {
        Self { inner: x }
    }
    fn to_borrow(&self) -> Predict<&[u8]> {
        Predict {
            inner: self.inner.as_slice(),
        }
    }
    fn to_borrow_mut(&mut self) -> Predict<&mut [u8]> {
        Predict {
            inner: &mut self.inner,
        }
    }
}

// use as_ref inside the func
fn take_both<T: AsRef<[u8]>>(x: &Predict<T>) -> usize {
    x.inner.as_ref().len()
}

// only for non-mut
fn take_ref(x: &Predict<&[u8]>) -> usize {
    x.inner.len()
}

// only for mut
fn take_ref_mut(x: &mut Predict<&mut [u8]>) -> usize {
    x.inner[3] = 10;
    x.inner.len()
}

/*
// impossible?
impl AsRef<Predict<&[u8]>> for Predict<Vec<u8>> {
    fn as_ref(&self) -> &Predict<&[u8]> {
        &self.to_borrow()
    }
}
 */

// 2. create two struct and use trait
// Second is return to_B() not to_B_ref()
struct A {
    inner: Vec<u8>,
}

impl A {
    fn to_b(&self) -> B {
        B { inner: &self.inner }
    }
}

trait HasVec {
    fn get_vec_ref(&self) -> &[u8];
}

impl HasVec for A {
    fn get_vec_ref(&self) -> &[u8] {
        &self.inner
    }
}

struct B<'a> {
    inner: &'a [u8],
}

impl<'a> HasVec for B<'a> {
    fn get_vec_ref(&self) -> &[u8] {
        self.inner
    }
}

fn take_both_A_B<T: HasVec>(x: &T) -> usize {
    x.get_vec_ref().len()
}

pub fn test1() {
    // first way
    let n = 5;
    //let mut x: Vec<u8> = Vec::with_capacity(n);
    //for ni in 0..n {
    //    x.push(ni as u8);
    //}
    let x = (0..n).map(|i| i as u8).collect::<Vec<u8>>();

    let a = Predict::new(x);
    //let pred = Predict { inner: x };

    let b = a.to_borrow();

    println!("{}", take_both(&a));
    println!("{}", take_both(&b));
    //println!("{}", take_ref(&a));
    println!("{}", take_ref(&b));

    // Is this heavy?
    // -> not heavy but not ideal too since [u8] is better
    //let c = &a;

    // what if for mutable? -> ok
    //let mut x: Vec<u8> = Vec::with_capacity(n);
    //for ni in 0..n {
    //    x.push(ni as u8);
    //}
    let x = (0..n).map(|i| i as u8).collect::<Vec<u8>>();
    let mut a = Predict::new(x);
    //let pred = Predict { inner: x };

    println!("{}", take_both(&mut a));

    let mut b = a.to_borrow_mut();

    println!("{}", take_both(&b));
    println!("{}", take_ref_mut(&mut b));
}

pub fn test2() {
    let a = A { inner: vec![1; 10] };
    take_both_A_B(&a);

    take_both_A_B(&a.to_b());
}

pub fn test() {
    test1();
    test2();
}
