//! what if inner is option

pub struct Predict<T: AsRef<[u8]>> {
    inner: Option<T>,
}

impl Predict<Vec<u8>> {
    // Self=Vec<u8>
    fn new(x: Vec<u8>) -> Self {
        Self { inner: Some(x) }
    }
    //fn to_borrow(&self) -> Predict<&[u8]> {
    fn to_borrow(&self) -> Result<Predict<&[u8]>, String> {
        match &self.inner {
            Some(x) => Ok(Predict {
                inner: Some(x.as_slice()),
            }),
            None => Err("error".to_owned()),
        }
        /*
        Predict {
            // error
            //inner: self.inner.as_ref().ok_or("error".to_owned()), // self.inner.as_slice(),
        }
         */
    }
}

/*
// impossible?
impl AsRef<Predict<&[u8]>> for Predict<Vec<u8>> {
    fn as_ref(&self) -> &Predict<&[u8]> {
        &self.to_borrow()
    }
}
 */

// use as_ref inside the func
fn take_both<T: AsRef<[u8]>>(x: &Predict<T>) -> usize {
    match &x.inner {
        Some(y) => y.as_ref().len(),
        None => 1,
    }
}

fn take_ref(x: &Predict<&[u8]>) -> usize {
    match x.inner {
        Some(y) => y.as_ref().len(),
        None => 1,
    }
}

pub fn test() {
    let n = 5;
    //let mut x: Vec<u8> = Vec::with_capacity(n);
    //for ni in 0..n {
    //    x.push(ni as u8);
    //}
    let x = (0..n).map(|i| i as u8).collect::<Vec<u8>>();

    let a = Predict::new(x);
    //let pred = Predict { inner: x };

    let b = a.to_borrow().unwrap();

    println!("{}", take_both(&a));
    println!("{}", take_both(&b));
    //println!("{}", take_ref(&a));
    println!("{}", take_ref(&b));

    // Is this heavy?
    // -> not heavy but not ideal too since [u8] is better
    //let c = &a;
}
