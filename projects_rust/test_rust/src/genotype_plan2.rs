//! create struct for B8
//!
//! Use this!!
//! simpler ver. is plan3
//!
//!
//! preparing TB8, VB8, MB8 is good since these should be separatedly used
//!
//!
//! seems great but not yet perfect...
//! creating Y and YRef are troublesome
//! -> but inevitable?
//!
//! This ref is great!! matrix
//! https://docs.rs/rulinalg/latest/src/rulinalg/matrix/mod.rs.html#44
//! https://athemathmo.github.io/rulinalg/doc/src/rulinalg/matrix/mod.rs.html#44
//! Matrix<T>, MatrixSlice<T>, MatrixSliceMut<T> are necessary
//! PhantomData is for 'a
//!
//! using *T might be better
//! -> plan4
//!
//!

type B8 = u8;

// regard as Matrix (2-dim)
// for three dim compressed matrix
// ax0 is compressed
// Vec, m, s, n
// if ys, m=s=1
/// Compressed Matrix
pub struct CMatrix {
    inner: Vec<B8>,
    // n
    ax0: usize,
    // s
    ax1: usize,
    // m
    ax2: usize,
}

//pub struct VB8(Vec<B8>, usize, usize, usize);
//pub struct VB8Ref<'a>(&'a [B8], usize);
// all &[B8] should be rewritten to VB8Ref
// mainly for loading and loss func
// for calculating loss, mut is not necessary
// -> create both mut and not, when loading it is mut and when returning it is ref
pub struct CMatrixMut<'a> {
    //inner: &'a [B8],
    inner: &'a mut [B8],
    ax0: usize,
    ax1: usize,
    ax2: usize,
}

pub struct CMatrixRef<'a> {
    inner: &'a [B8],
    ax0: usize,
    ax1: usize,
    ax2: usize,
}

/*
 // same as CVec and digit=1
 // not bad creating CBits too
pub struct CBits(CMatrix);
pub struct CBitsRef<'a>(CMatrixRef<'a>);
pub struct CBitsMut<'a>(CMatrixMut<'a>);
 */

pub struct CVec(CMatrix);
pub struct CVecRef<'a>(CMatrixRef<'a>);
pub struct CVecMut<'a>(CMatrixMut<'a>);

trait BaseCMatrix {
    fn inner(&self) -> &[B8];
    fn ax0(&self) -> usize;
}
trait BaseCVec {}
trait BaseCBits {}

// for mut
// T: BaseMatrix
trait BaseRef {}
// T: BaseVB
trait BaseMut {}

//pub struct Y(CBits);
// T=VB8Ref, VB8RefMut
pub struct YRef<T>(T);
/*
pub struct YRef<'a>(VB8Ref<'a>);
pub struct YRefMut<'a>(VB8RefMut<'a>);
 */

// for genot
pub struct GenotTwin(CMatrix);
/*
pub struct GenotTwinRef<'a>(TB8Ref<'a>);
/// one snv
pub struct GenotTwinRefM<'a>(MB8Ref<'a>);
/// one snv and one model
pub struct GenotTwinRefMS<'a>(VB8Ref<'a>);
 */

// T=TB8Ref, MB8Ref, VB8Ref, TB8RefMut, MB8RefMut, VB8RefMut
pub struct GenotTwinRef<T>(T);

impl CMatrix {
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
    fn to_ref_mut(&mut self) -> CMatrixMut {
        CMatrixMut {
            inner: &mut self.inner,
            ax0: self.ax0,
            ax1: self.ax1,
            ax2: self.ax2,
        }
    }
}

impl<'a> CMatrixMut<'a> {
    fn inner_mut(&mut self) -> &mut [B8] {
        self.inner
        //self.inner
        //&self.0
    }
    fn access(&self, mi: usize, si: usize, ni: usize) -> bool {
        self.inner[ni] == 1
    }
    fn sub_mut(&mut self, mi_start: usize, mi_end: usize) -> CMatrixMut {
        CMatrixMut {
            inner: &mut self.inner[mi_start..mi_end],
            ax0: self.ax0,
            ax1: self.ax1,
            ax2: self.ax2,
        }
    }
    fn set(&mut self, x: bool, mi: usize, si: usize, ni: usize) {
        self.inner[ni] = x as u8;
        let pred = self.inner_mut();
        pred[ni] = x as u8;
    }
}

impl<'a> GenotTwinRef<CMatrixMut<'a>> {
    fn a(&self) -> bool {
        true
    }
    fn inner(&self) -> &[B8] {
        self.0.inner
    }
    fn inner_mut(&mut self) -> &mut [B8] {
        self.0.inner
    }
    fn to_GMB8(&'a mut self) -> GenotTwinRef<CVecMut<'a>> {
        let m = self.0.ax0;
        let s = self.0.ax1;
        let n = self.0.ax2;
        //let inner = self.0.inner;
        let inner = self.inner_mut();
        let t = CMatrixMut {
            inner,
            ax0: m,
            ax1: s,
            ax2: n,
        };
        let m = CVecMut(t);
        GenotTwinRef(m)
    }

    // OK!!
    // 'a of returned value should be removed!!
    // What is the lifetime of returned value???
    fn sub_mut(&mut self, mi_start: usize, mi_end: usize) -> GenotTwinRef<CMatrixMut> {
        let m = self.0.ax0;
        let s = self.0.ax1;
        let n = self.0.ax2;
        //let inner = self.0.inner;
        let inner = self.inner_mut();
        let t = CMatrixMut {
            inner: &mut inner[mi_start..mi_end],
            ax0: m,
            ax1: s,
            ax2: n,
        };
        GenotTwinRef(t)
    }
    fn access(&'a mut self, mi: usize) -> u8 {
        self.inner()[mi]
    }
    fn set_val(&'a mut self, val: u8, mi: usize) {
        self.inner_mut()[mi] = val;
    }
    fn set_val2<'b>(&'b mut self, val: u8, mi: usize) {
        self.inner_mut()[mi] = val;
    }
}

pub struct SampleTwin {
    y: CMatrix,
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

fn set_val<'a>(g: &'a mut GenotTwinRef<CMatrixMut<'a>>, val: u8, mi_start: usize, mi_end: usize) {
    g.sub_mut(mi_start, mi_end).set_val(val, 0);
}

fn f(v: &mut [u8]) {
    v[0] = 0;
}

/*
pub struct ChunksMut<'a, T: 'a> {
    v: &'a mut [T],
    chunk_size: usize,
}

// https://doc.rust-lang.org/src/core/slice/iter.rs.html#1583-1586
impl<'a, T> Iterator for ChunksMut<'a, T> {
    type Item = &'a mut [T];

    fn next(&mut self) -> Option<&'a mut [T]> {
        if self.v.is_empty() {
            None
        } else {
            let sz = std::cmp::min(self.v.len(), self.chunk_size);
            let tmp = std::mem::replace(&mut self.v, &mut []);
            //let (head, tail) = tmp.split_at_mut(sz);
            let tmp = &mut self.v;
            let (head, tail) = tmp.split_at_mut(sz);
            self.v = tail;
            Some(head)
        }
    }
}
 */

pub fn test() {
    /*
    let mut v = vec![1, 2, 3, 4, 5, 6, 7, 8];
    let a = ChunksMut {
        v: &mut v,
        chunk_size: 2,
    };

    for x in a {
        x[0] = 100;
    }
    println!("x {:?}", v);
     */

    /*
    let preds = TB8 {
        inner: vec![1, 2, 3, 4, 5, 6],
        ax0: 3,
        ax1: 2,
        ax2: 4,
    };
    let mut genot = GenotTwin(preds);
    println!("access genot {}", genot.extract(0, 0, 0));

    let tb8 = TB8RefMut {
        inner: &mut genot.0.inner,
        m: genot.0.ax2,
        s: genot.0.ax1,
        n: genot.0.ax2,
    };
    let mut genot_ref = GenotTwinRef(tb8);
    //let genot_ref_mut = &mut genot_ref;

    //println!("access genot {}", genot_ref.access(0));
    //println!("access genot {}", (&mut genot_ref).access(0));
    //println!("genot_ref {}", &genot_ref);
    //genot_ref.sub_mut(0, 2).inner_mut()[0] = 8;
    //println!("access2 genot {}", genot_ref.access(0));

    let mut v = vec![10; 100];
    f(&mut v);
    f(&mut v);

    genot_ref.sub_mut4(0, 1);
    genot_ref.sub_mut4(0, 1);

    //let genot_ref_mut = &mut genot_ref;
    {
        //genot_ref.set_val2(0, 1);
        genot_ref.sub_mut3(0, 1);
        //genot_ref.sub_mut(0, 1);
    }
    {
        //genot_ref.set_val2(0, 1);
        genot_ref.sub_mut3(0, 1);
    }

    //genot_ref.sub_mut(0, 1).set_val(8, 0);
    //genot_ref.sub_mut(0, 1).set_val(8, 0);

    for _ in [0, 1].iter() {
        genot_ref.sub_mut3(0, 1).set_val(8, 0);
    }

    // how to do this...?
    [0, 2, 4].iter().for_each(|mi| {
        //let mut genot_ref = &mut genot_ref;
        genot_ref.sub_mut3(0, 1).set_val(8, 0);
    })

    /*
    for _ in [0, 1].iter() {
        genot_ref.sub_mut(0, 1).set_val(8, 0);
    }

    // how to do this...?
    [0, 2, 4].iter().for_each(|mi| {
        //let mut genot_ref = &mut genot_ref;
        genot_ref.sub_mut(0, 1).set_val(8, 0);
    })
    //.for_each(|&mi| set_val(genot_ref_mut, 8, mi, mi + 2))
     */
     */
}
