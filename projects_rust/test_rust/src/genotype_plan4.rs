//! create struct for B8
//!
//! Not bad but plan2,3 is more easy.
//!
//!
//!
//!
//!
//! TODO:
//! iterator
//! paralleliterator
//!
//!
//! refer to Matrix in rulinalg
//!
//!
//! Use pointer might be easier? but should use unsafe
//!
//! preparing TB8, VB8, MB8 is good since these should be separatedly used
//!
//!
//! This ref is great!! matrix
//! https://docs.rs/rulinalg/latest/src/rulinalg/matrix/mod.rs.html#44
//! https://github.com/AtheMathmo/rulinalg
//! Matrix<T>, MatrixSlice<T>, MatrixSliceMut<T> are necessary
//! PhantomData is for 'a
//!
//! Most func is in base
//! https://docs.rs/rulinalg/0.4.2/src/rulinalg/matrix/base/mod.rs.html
//!
//!
//!
//!
//! using *T might be better
//!
//! Ndarray
//! https://docs.rs/ndarray/0.15.4/src/ndarray/lib.rs.html#1268-1282
//! IntoParallelIterator
//! https://docs.rs/ndarray/0.15.4/src/ndarray/parallel/par.rs.html#22-25
//!
//!
//! Both uses pointer
//!
//!
//!
//! Iterator
//! https://aloso.github.io/2021/03/09/creating-an-iterator
//!
//!

use rayon::iter::plumbing::*;
use rayon::prelude::*;
use std::marker::PhantomData;

type B8 = u8;

// for three dim matrix
// Vec, m, s, n
// if ys, m=s=1
pub struct TB8 {
    inner: Vec<B8>,
    // this order for VB8
    // n
    ax0: usize,
    // s
    ax1: usize,
    // m
    ax2: usize,
}

impl TB8 {
    fn as_ptr(&self) -> *const B8 {
        self.inner.as_ptr()
    }
    fn to_tb8mut(&mut self) -> TB8Mut {
        TB8Mut {
            ptr: self.as_ptr() as *mut B8,
            //ptr: self.as_ptr(),
            ax0: self.ax0,
            ax1: self.ax1,
            ax2: self.ax2,
            marker: PhantomData::<&mut B8>,
        }
    }
}

// for both ref and mut
trait BaseMatrix {}
// for mut
// T: BaseMatrix
trait BaseMatrixMut {}

// whole ref or
// part of mi_begin ~ mi_end
// ptr at mi_begin, ax0=mi_end-mi_begin
#[derive(Copy, Clone)]
pub struct TB8Mut<'a> {
    // ptr is at mi_start
    ptr: *mut B8,
    //ptr: *const B8, //for ref
    ax0: usize,
    ax1: usize,
    ax2: usize,
    marker: PhantomData<&'a mut B8>,
}

impl<'a> TB8Mut<'a> {
    fn ax0(&self) -> usize {
        self.ax0
    }
    fn set_ax0(&mut self) {
        self.ax0 = 5;
    }
}

// ax1=ax2=0
pub struct VB8(TB8);
pub struct VB8Mut<'a>(TB8Mut<'a>);

pub struct Genot(TB8);
pub struct GenotRef<'a>(TB8Mut<'a>);

impl<'a> GenotRef<'a> {
    fn n(&self) -> usize {
        self.0.ax0()
    }

    // ok
    fn set_n(&mut self) {
        self.0.set_ax0()
    }

    fn set_n2(&mut self) {
        self.inner_mut().set_ax0()
    }

    fn inner_mut(&mut self) -> &mut TB8Mut<'a> {
        &mut self.0
    }

    fn tb8mut(&mut self) -> &mut TB8Mut<'a> {
        &mut self.0
    }
    // same as tb8mut()
    // unncessary: where 'a: 'b,
    // x since no way 'b is longer than 'a
    // -> lifetime::tb8mut_52() does compile
    fn tb8mut2<'b>(&'b mut self) -> &'b mut TB8Mut<'a> {
        &mut self.0
    }
    fn tb8mut3<'b>(&'b mut self) -> &'b mut TB8Mut<'a>
    where
        'a: 'b,
    {
        &mut self.0
    }
}

fn f(tb8mut: &mut TB8Mut) {
    tb8mut.ax0 = 5
}

pub fn test() {
    let tb8 = TB8 {
        inner: vec![1, 2, 3, 4, 5, 6],
        ax0: 3,
        ax1: 2,
        ax2: 4,
    };

    println!("TB8 {:?}", tb8.inner);

    let mut g = Genot(tb8);
    let mut gref = GenotRef(g.0.to_tb8mut());
    println!("n==3 {:?}", gref.n());

    gref.set_n2();
    gref.set_n2();
    //gref.set_n();
    //gref.set_n();
    println!("n==5 {:?}", gref.n());

    f(gref.tb8mut3());
    f(gref.tb8mut3());
}
