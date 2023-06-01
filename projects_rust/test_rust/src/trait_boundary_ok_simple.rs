//! create SampleTwin but y is YB8 not Vec<B8>
//! complicated since we need VB8, VB8Ref, GenotTwinM...
//! Can I make them integrated?
//! -> just prepare VB8 and GenotTwinM

type B8 = u8;

// for y
pub struct VB8(Vec<B8>, usize);
//pub struct VB8Ref<'a>(&'a [B8], usize);

impl VB8 {
    fn inner(&self) -> &Vec<B8> {
        &self.0
    }
    fn access(&self, ni: usize) -> u8 {
        // assert ni < self.1
        self.inner()[0]
        // bget()
    }
}

pub struct SampleTrait {
    y: VB8,
    cov: Vec<u8>,
}

// for genot
pub struct GenotTwin {
    genot_v: Vec<B8>,
    m: usize,
    n: usize,
}

pub struct GenotTwinM<'a>(&'a [B8], usize);

impl<'a> GenotTwinM<'a> {
    pub fn access(&self, ni: usize) -> u8 {
        self.0[ni]
    }
}

/*
impl<'a> GenotTwinM<'a> {
    pub fn to_vb8ref(&self) -> VB8Ref {
        VB8Ref(self.0, self.1)
    }
}
 */

impl GenotTwin {
    pub fn to_ref(&self) -> GenotTwinM {
        GenotTwinM(&self.genot_v, self.n)
    }
}

pub fn loss(gr: &GenotTwinM) -> f64 {
    gr.0[0] as f64
}

pub fn test() {
    let g = GenotTwin {
        genot_v: vec![1, 0],
        m: 10,
        n: 3,
    };
    let gr = g.to_ref();
    let loss = loss(&gr);
    println!("loss {}", loss);
}
