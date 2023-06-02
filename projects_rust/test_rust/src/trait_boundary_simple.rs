// trait

pub struct Genot {}
pub trait GT {}
impl GT for Genot {}

pub struct Dataset<G, V, S> {
    g: G,
    v: V,
    s: S,
}
pub trait DT {
    fn access(&self) -> u8;
}
impl<G, V, S> DT for Dataset<G, V, S>
where
    G: GT,
{
    fn access(&self) -> u8 {
        3
    }
}

pub fn test() {
    let g = Genot {};
    let d = Dataset { g, v: 0, s: 0 };
    d.access();
}
