// trait

type B8 = u8;

use std::{marker::PhantomData, vec};

pub struct YB8(Vec<B8>);
pub struct Yu8(Vec<u8>);

pub trait IsYTrait {}
impl IsYTrait for YB8 {}
impl IsYTrait for Yu8 {}

pub struct Genot<Bi>
where
    Bi: IsYTrait,
{
    raw: Bi,
}
pub trait GT<Bi: IsYTrait> {
    fn raw_gt(&self) -> &Bi;
}
impl<Bi: IsYTrait> GT<Bi> for Genot<Bi> {
    fn raw_gt(&self) -> &Bi {
        &self.raw
    }
}

pub struct Dataset<G, V, S, Bi> {
    g: G,
    v: V,
    s: S,
    bi: PhantomData<Bi>,
}
pub trait DT {
    type V;
    fn access(&self) -> u8;
    fn slice(&self) -> &Self::V;
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
    type V = Bi;
    fn access(&self) -> u8 {
        3
    }
    fn slice(&self) -> &Self::V {
        self.raw()
    }
}

pub fn test() {
    let raw = YB8(vec![0]);
    let g = Genot { raw };
    let d: Dataset<_, _, _, YB8> = Dataset {
        g,
        v: 0,
        s: 0,
        bi: PhantomData,
    };
    d.access();
}
