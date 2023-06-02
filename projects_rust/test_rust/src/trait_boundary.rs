//! getting better, but IsYTrait::access() needs m, which in in Genot, so access() should be in Genot

type B8 = u8;

use std::marker::PhantomData;

pub struct VB8(Vec<B8>);
pub struct Vu8(Vec<u8>);

// for exome
// Vec<B8> for dense matrix and Vec<usize> for sparse
pub struct VE8(Vec<B8>, Vec<usize>);

pub trait IsYTrait {
    type Y;
    fn inner(&self) -> &Self::Y;
    // for index access not sample_i
    //fn access_i(&self, i: usize) -> u8;
    fn access(&self, ni: usize, n: usize) -> u8;
}
impl IsYTrait for VB8 {
    type Y = Vec<B8>;
    fn inner(&self) -> &Self::Y {
        &self.0
    }
    fn access(&self, ni: usize, n: usize) -> u8 {
        self.inner()[ni + n]
    }
}
impl IsYTrait for Vu8 {
    type Y = Vec<u8>;
    fn inner(&self) -> &Self::Y {
        &self.0
    }
    fn access(&self, ni: usize, n: usize) -> u8 {
        self.inner()[ni + n]
    }
}

impl IsYTrait for VE8 {
    // for generalized returned values
    type Y = Vec<u8>;
    fn inner(&self) -> &Self::Y {
        &self.0
    }
    fn access(&self, ni: usize, n: usize) -> u8 {
        self.inner()[ni + n]
    }
}

impl VE8 {
    // split for sparse and dense
    fn a(&self) -> usize {
        1
    }
}

pub struct Genot<Bi>
where
    Bi: IsYTrait,
{
    raw: Bi, //VB8 or Vu8
    n: usize,
}
pub trait GT<Bi: IsYTrait> {
    fn raw_gt(&self) -> &Bi;
    fn raw_gt_in(&self) -> &Bi::Y;
    fn access(&self, ni: usize) -> u8;
}
impl<Bi: IsYTrait> GT<Bi> for Genot<Bi> {
    fn raw_gt(&self) -> &Bi {
        &self.raw
    }
    //fn raw_gt_in(&self) -> &<Bi as IsYTrait>::Y {
    fn raw_gt_in(&self) -> &Bi::Y {
        &self.raw.inner()
    }
    fn access(&self, ni: usize) -> u8 {
        // here needs n in Genot
        self.raw_gt().access(ni, self.n)
    }
}

pub struct Dataset<G, V, S, Bi> {
    g: G,
    v: V,
    s: S,
    bi: PhantomData<Bi>,
}

impl<G, V, S, Bi> Dataset<G, V, S, Bi> {
    fn g(&self) -> &G {
        &self.g
    }
}

pub trait DT
where
    Self::X: IsYTrait,
{
    type X; //=Bi
    fn slice(&self) -> &Self::X;
    // let user do
    //fn slice_in(&self) -> &<<Self as DT>::X as IsYTrait>::Y;
    //fn slice_in(&self) -> &Self::V::Y;
    fn access(&self, ni: usize) -> u8;
}

impl<G: GT<Bi>, V, S, Bi: IsYTrait> Dataset<G, V, S, Bi> {
    fn raw(&self) -> &Bi {
        &self.g().raw_gt()
    }
}

// should separate imple for Bi
impl<V, S, Bi> DT for Dataset<Genot<Bi>, V, S, Bi>
where
    Bi: IsYTrait,
{
    type X = Bi;
    fn slice(&self) -> &Self::X {
        self.raw()
    }
    /*
    // if you want this, you can just return slice() and let user write .slice().inner
    fn slice_in(&self) -> &<<Self as DT>::X as IsYTrait>::Y {
        self.slice().inner()
        //= self.raw().inner()
    }
     */
    fn access(&self, ni: usize) -> u8 {
        // ex. indicate snv_index and sample_id and get value
        self.g().access(ni)
    }
}

pub fn test() {
    let raw = VB8(vec![0; 10]);
    let g = Genot { raw, n: 3 };
    let d: Dataset<_, _, _, VB8> = Dataset {
        g,
        v: 0,
        s: 0,
        bi: PhantomData,
    };
    d.access(3);
}
