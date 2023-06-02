//! Iterator
//!
//! iterator for plan4
//! but using ptr is a littele nervous
//! Is there anything that couldn't be achieved by plan2??
//! I think no...
//!
//! best
//! https://aloso.github.io/2021/03/09/creating-an-iterator
//!
//! https://stackoverflow.com/questions/30218886/how-to-implement-iterator-and-intoiterator-for-a-simple-struct
//! matrix
//! https://athemathmo.github.io/rulinalg/doc/src/rulinalg/matrix/iter.rs.html#11-36
//!
//!
//! Ndarray
//! IntoParallelIterator
//! https://docs.rs/ndarray/0.15.4/src/ndarray/parallel/par.rs.html#22-25
//!

use rayon::iter::plumbing::*;
use rayon::iter::ParallelBridge;
use rayon::prelude::ParallelIterator;
use rayon::prelude::*;
use std::marker::PhantomData;
use std::marker::{Send, Sync};
use std::sync::mpsc::channel;

type B8 = u8;

struct TB8 {
    inner: Vec<u8>,
    m: usize,
}

impl TB8 {
    fn as_ptr(&mut self) -> *mut B8 {
        // or as_mut_ptr()?
        //self.inner.as_ptr() as *mut B8
        self.inner.as_mut_ptr()
    }
    fn to_tb8mut(&mut self) -> TB8Mut {
        TB8Mut {
            ptr: self.as_ptr(),
            m: self.m,
            marker: PhantomData::<&mut B8>,
        }
    }
}

impl TB8 {
    fn iter<'a>(&'a mut self) -> TB8MutIter<'a> {
        TB8MutIter {
            ptr: self.as_ptr(),
            m: self.m,
            mi_start: 0,
            marker: PhantomData::<&'a mut B8>,
        }
    }
}

#[derive(Copy, Clone)]
pub struct TB8Mut<'a> {
    // ptr is at mi_start
    ptr: *mut B8,
    m: usize,
    marker: PhantomData<&'a mut B8>,
}

pub struct TB8MutIter<'a> {
    // ptr is at mi_start
    ptr: *mut B8,
    m: usize,
    mi_start: usize,
    marker: PhantomData<&'a mut B8>,
}

// similar is Iterator for ColsMut
// https://docs.rs/rulinalg/latest/src/rulinalg/matrix/iter.rs.html#113-168
// It is impossible to use different lifetime for returned value
// but it should be ok since all returned values are independent
// How to run in rayon??
impl<'a> Iterator for TB8MutIter<'a> {
    type Item = TB8Mut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi_start < self.m {
            let mi = self.mi_start;
            self.mi_start += 1;
            let m = self.m;
            self.m = m - 1;
            let ptr = self.ptr;
            unsafe {
                let val = self.ptr.read();
                self.ptr = self.ptr.offset(1);

                return Some(TB8Mut {
                    ptr,
                    m,
                    marker: PhantomData::<&'a mut B8>,
                });
            }
        }
        None
    }
    // implement size_hint() is better
}

/*
// perfect for item=u8
impl<'a> Iterator for TB8MutIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi_start < self.m {
            let mi = self.mi_start;
            self.mi_start += 1;
            let ptr = self.ptr;
            unsafe {
                let val = self.ptr.read();
                self.ptr = self.ptr.offset(1);
                return Some(val);
            }
        }
        None
    }
}
 */

pub fn test() {
    let a = vec![1, 2, 3, 4, 5, 6];
    let mut g = TB8 { inner: a, m: 6 };

    println!("iter()");
    for x in g.iter() {
        println!("x {}", x.m);
    }

    /*
    println!("par_bridge()");
    for x in g.iter().par_bridge() {
        println!("x {}", x.m);
    }
     */

    // impl<'a> IntoIterator for &'a Genot
    //for x in g.into_iter() {
    //for x in &g {
    //    println!("x {}", x);
    //}

    let a = vec![1, 2, 3, 4, 5, 6];
    let g = TB8 { inner: a, m: 6 };

    // impl IntoIterator for Genot
    //for x in g.into_iter() {
    //    println!("x {}", x);
    //}

    //println!("g {}", g);
}
