//! lifetime
//!
//! I think now I got the point
//!
//! ```ignore
//! struct T(U<'a>);
//!
//!  impl<'a> T<'a>{
//!  fn<'b,'c>(&'b self) -> &'c S<'a>{}
//! }
//! ```
//! 'a is lifetime of U.
//! 'b is lifetime of (input) T.
//! 'c is lifetime of S.
//! The point is 'c should not be longer than 'a or 'b.
//!
//!
//! Q. How come tb8mut_52() can compile??
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
            ptr: self.as_ptr(),
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
trait BaseMatrixMut {}

// whole ref or
// part of mi_begin ~ mi_end
// ptr at mi_begin, ax0=mi_end-mi_begin
#[derive(Copy, Clone)]
pub struct TB8Mut<'a> {
    // ptr is at mi_start
    ptr: *const B8,
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

    // This is only pattern that I found to call twice
    // what the lifetime would be if not elision
    // -> same as tb8mut_53()
    fn tb8mut(&mut self) -> &mut TB8Mut<'a> {
        &mut self.0
    }
    // how to make the returned lifetime less than 'a?
    // so that be able to call tb8mut() twice in a row
    // impossible? since GenotRef<'a> has lifetime 'a
    // self and ret are the same
    // ng
    //fn tb8mut_2<'b: 'a>(&'b mut self) -> &'b mut TB8Mut<'b> where 'b: 'a{
    fn tb8mut_2<'b>(&'b mut self) -> &'b mut TB8Mut<'b>
    where
        // here said left > right
        // https://blog-mk2.d-yama7.com/2020/12/20201230_rust_lifetime/
        //
        // xx 'b < 'a
        // xx 'b is shorter than 'a
        // 'b > 'a
        // 'b is longer than 'a
        'b: 'a,
        //'a: 'b, // error
    {
        &mut self.0
    }
    // ng
    fn tb8mut_4(&'a mut self) -> &'a mut TB8Mut<'a> {
        &mut self.0
    }
    // ng
    // why not?
    fn tb8mut_51<'b>(&'a mut self) -> &'b mut TB8Mut<'a>
    where
        'a: 'b,
    {
        &mut self.0
    }

    // ng
    // HOW COME this can compile??
    // Is it possible to 'b > 'a ?
    fn tb8mut_52<'b>(&'b mut self) -> &'b mut TB8Mut<'a>
    where
        'b: 'a,
    {
        &mut self.0
    }
    // !!!!!!!!!!!!!!!!!!
    // This!!!
    fn tb8mut_53<'b>(&'b mut self) -> &'b mut TB8Mut<'a>
    where
        // 'a is longer than 'b
        'a: 'b,
    {
        &mut self.0
    }
    // This should be non-omitted ver. of tb8mut()
    // returned 'b can be the same as input 'b since input 'b is deleted after this function is called
    // !!!!!!!!!!!!!!!
    // same as 53?
    fn tb8mut_54<'b>(&'b mut self) -> &'b mut TB8Mut<'a> {
        &mut self.0
    }

    // ok
    fn tb8mut_61<'b, 'c>(&'c mut self) -> &'b mut TB8Mut<'a>
    where
        // 'a > 'c > 'b
        'a: 'c,
        'c: 'b,
    {
        &mut self.0
    }

    /*
    // error
    fn tb8mut_8<'b, 'c, 'd>(&'c mut self) -> &'b mut TB8Mut<'d>
    where
        // 'b < 'c < 'd < 'a
        'a: 'd,
        'd: 'c,
        'c: 'b,
    {
        &mut self.0
    }
     */
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

    f(gref.tb8mut());
    f(gref.tb8mut());

    f(gref.tb8mut_61());
    f(gref.tb8mut_61());

    f(gref.tb8mut_53());
    f(gref.tb8mut_53());

    /*
    // this does not change since copy
    let mut tb8mut = gref.0;
    f(&mut tb8mut);
    println!("n==5 {:?}", tb8mut.ax0());
    println!("n==5 {:?}", gref.n());
     */
}
