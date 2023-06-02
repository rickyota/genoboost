//! getting better, but IsYTrait::access() needs m, which in in Genot, so access() should be in Genot

type B8 = u8;

use std::marker::PhantomData;

pub struct VB8(Vec<B8>);
pub struct Vu8(Vec<u8>);

pub trait IsYTrait {
    type Y;
    fn inner(&self) -> &Self::Y;
    fn access(&self, ni: usize) -> u8;
}
impl IsYTrait for VB8 {
    type Y = Vec<B8>;
    fn inner(&self) -> &Self::Y {
        &self.0
    }
    fn access(&self, ni: usize) -> u8 {
        self.inner()[0]
    }
}
impl IsYTrait for Vu8 {
    type Y = Vec<u8>;
    fn inner(&self) -> &Self::Y {
        &self.0
    }
    fn access(&self, ni: usize) -> u8 {
        self.inner()[0]
    }
}

pub struct Genot<Bi>
where
    Bi: IsYTrait,
{
    raw: Bi, //VB8 or Vu8
}
pub trait GT<Bi: IsYTrait> {
    fn raw_gt(&self) -> &Bi;
    fn raw_gt_in(&self) -> &Bi::Y;
}
impl<Bi: IsYTrait> GT<Bi> for Genot<Bi> {
    fn raw_gt(&self) -> &Bi {
        &self.raw
    }
    //fn raw_gt_in(&self) -> &<Bi as IsYTrait>::Y {
    fn raw_gt_in(&self) -> &Bi::Y {
        &self.raw.inner()
    }
}

pub struct Dataset<G, V, S, Bi> {
    g: G,
    v: V,
    s: S,
    bi: PhantomData<Bi>,
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
        &self.g.raw_gt()
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
        // TODO: self.slice().inner().get(si, mi, ni)
        self.slice().access(ni) as u8
    }
}

pub fn test() {
    let raw = VB8(vec![0]);
    let g = Genot { raw };
    let d: Dataset<_, _, _, VB8> = Dataset {
        g,
        v: 0,
        s: 0,
        bi: PhantomData,
    };
    d.access(3);
}
