// How to implement genot trait

struct Mat(Vec<u8>);
struct MatRef<'a>(&'a [u8]);

trait BaseMat {
    fn get(&self, ni: usize) -> u8;
}

impl BaseMat for Mat {
    fn get(&self, ni: usize) -> u8 {
        self.0[ni]
    }
}

impl<'a> BaseMat for MatRef<'a> {
    fn get(&self, ni: usize) -> u8 {
        self.0[ni]
    }
}

struct Genot(Mat);
struct GenotRef<'a>(MatRef<'a>);

trait BaseGenot {
    type Inner: BaseMat;
    fn genot_inner(&self) -> &Self::Inner;

    fn get(&self, ni: usize) -> u8 {
        self.genot_inner().get(ni)
    }
}

/* trait BaseGenot {
    fn genot_inner<T>(&self) -> T
    where
        T: BaseMat;

    fn get<T>(&self, ni: usize) -> u8
    where
        T: BaseMat,
    {
        self.genot_inner::<T>().get(ni)
    }
} */

/* trait BaseGenot {
    fn genot_inner<T>(&self) -> T
    where
        T: BaseMat;

    fn get<T>(&self, ni: usize) -> u8
    where
        T: BaseMat,
    {
        self.genot_inner::<T>().get(ni)
    }
} */

impl BaseGenot for Genot {
    type Inner = Mat;
    fn genot_inner(&self) -> &Self::Inner {
        &self.0
    }

    /*     fn genot_inner(&self) -> Mat {
        self.0
    } */
}

pub fn test() {
    let v = vec![1, 2, 3];
    let g = Genot(Mat(v));
    println!("get {}", g.get(2));
}
