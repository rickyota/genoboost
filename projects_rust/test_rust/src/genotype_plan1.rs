//! create struct for B8
//!
//! directly use MB8 in ys is troublesome
//!  -> genotype_plan2.rs

type B8 = u8;

// for three dim matrix
// Vec, m, s, n
// if ys, m=s=1
pub struct MB8 {
    inner: Vec<B8>,
    m: usize,
    s: usize,
    n: usize,
}
//pub struct VB8(Vec<B8>, usize, usize, usize);
//pub struct VB8Ref<'a>(&'a [B8], usize);
// all &[B8] should be rewritten to VB8Ref
// mainly for loading and loss func
pub struct MB8Ref<'a> {
    inner: &'a mut [B8],
    m: usize,
    s: usize,
    n: usize,
}

// for genot
// GenotTwinRef is unnecessary
pub struct GenotTwin(MB8);
//pub struct GenotTwin {
//    genot_v: Vec<B8>,
//    m: usize,
//    n: usize,
//}

/*
// not necessary
pub struct VB8 {
    pred: Vec<B8>,
    n: usize,
}
pub struct GB8 {
    pred: Vec<B8>,
    s: usize,
    n: usize,
}
 */

impl MB8 {
    fn inner(&self) -> &Vec<B8> {
        &self.inner
        //&self.0
    }
    fn inner_mut(&mut self) -> &mut Vec<B8> {
        &mut self.inner
        //&self.0
    }
    fn access(&self, mi: usize, si: usize, ni: usize) -> bool {
        // assert ni < self.1
        self.inner[ni] == 1
        // bget()
    }
    fn set(&mut self, x: bool, mi: usize, si: usize, ni: usize) {
        let pred = self.inner_mut();
        pred[ni] = x as u8;
    }
    fn to_ref_mut(&mut self) -> MB8Ref {
        MB8Ref {
            inner: &mut self.inner,
            m: self.m,
            s: self.s,
            n: self.n,
        }
    }
}

impl<'a> MB8Ref<'a> {
    fn inner_mut(&mut self) -> &mut [B8] {
        self.inner
        //&self.0
    }
    fn access(&self, mi: usize, si: usize, ni: usize) -> bool {
        self.inner[ni] == 1
    }
    fn sub_mut(&mut self, mi_start: usize, mi_end: usize) -> MB8Ref {
        MB8Ref {
            inner: &mut self.inner[0..1],
            m: self.m,
            s: self.s,
            n: self.n,
        }
    }
    fn set(&mut self, x: bool, mi: usize, si: usize, ni: usize) {
        self.inner[ni] = x as u8;
        let pred = self.inner_mut();
        pred[ni] = x as u8;
    }
}

pub struct SampleTwin {
    y: MB8,
    cov: Vec<u8>,
}

impl SampleTwin {
    fn set_y(&mut self, x: bool, ni: usize) {
        self.y.set(x, 1, 1, ni);
    }
}

impl GenotTwin {
    fn extract(&self, mi: usize, si: usize, ni: usize) -> bool {
        self.0.access(mi, si, ni)
    }
}

pub fn test() {
    let y = MB8 {
        inner: vec![1, 2, 3],
        m: 1,
        s: 1,
        n: 20,
    };
    let mut sample = SampleTwin { y, cov: vec![] };
    sample.set_y(false, 2);

    let preds = MB8 {
        inner: vec![1, 2, 3, 4, 5, 6],
        m: 3,
        s: 2,
        n: 4,
    };
    let mut genot = GenotTwin(preds);
    println!("access genot {}", genot.extract(0, 0, 0));

    let x = genot.extract(2, 1, 1);

    let mut genot_ref = genot.0.to_ref_mut();
    let mut genot_sub_ref = genot_ref.sub_mut(0, 1);
    println!("access {}", genot_sub_ref.access(0, 0, 0));
    genot_sub_ref.set(false, 0, 0, 0);
    println!("access afr set {}", genot_sub_ref.access(0, 0, 0));
    println!("access genot afr set {}", genot.extract(0, 0, 0));

    //let g = GenotTwin {
    //    genot_v: vec![1, 0],
    //    m: 10,
    //    n: 3,
    //};
    //let gr = g.to_ref();
    //let loss = loss(&gr);
    //println!("loss {}", loss);
}
