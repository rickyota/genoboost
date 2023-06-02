//! Simple ver. for plan2
//! To investigate mut second time
//!
//! Now, all problems are solved!!
//! Use this structure!!!
//!
//!
//! Not and ideal idea
//! 1. This is not ideal even if this can be achieved since &mut slice can be created once at a time
//! and never be used by ParallelIterator.
//! -> DONE: iter_mut() is unsafe, and I did implement it.
//! -> DONE: using par_bridge()
//! TODO LATER: paralleliterator
//!
//! 2. Not easy to implement IntoIterator
//! -> able but not necessary now
//!
//!

use rayon::prelude::*;

type B8 = u8;

// for three dim matrix
// Vec, m, s, n
// if ys, m=s=1
pub struct TB8 {
    inner: Vec<B8>,
    m: usize,
}

pub struct TB8RefMut<'a> {
    //inner: &'a [B8],
    inner: &'a mut [B8],
    m: usize,
}

pub struct GB8(TB8);

pub struct GB8Ref<'a>(TB8RefMut<'a>);

impl TB8 {
    fn iter_mut(&mut self) -> TB8IterMut {
        TB8IterMut {
            tb8refmut: self,
            mi: 0,
        }
    }
    fn val_iter_mut(&mut self) -> TB8ValIterMut {
        TB8ValIterMut {
            tb8refmut: self,
            mi: 0,
        }
    }
}

impl<'a> TB8RefMut<'a> {
    fn inner_mut(&mut self) -> &mut [B8] {
        self.inner
    }
    fn access(&self, mi: usize) -> u8 {
        self.inner[mi]
    }
    fn sub_mut3(&mut self, mi_start: usize, mi_end: usize) -> TB8RefMut {
        let tmp: &mut [u8] = std::mem::replace(&mut &mut self.inner[mi_start..mi_end], &mut []);
        TB8RefMut {
            inner: tmp,
            //inner: &mut self.inner[mi_start..mi_end],
            m: self.m,
        }
    }
    fn sub_mut(&mut self, mi_start: usize, mi_end: usize) -> TB8RefMut {
        TB8RefMut {
            inner: &mut self.inner[mi_start..mi_end],
            m: self.m,
        }
    }
    fn sub_mut2(&'a mut self, mi_start: usize, mi_end: usize) -> TB8RefMut {
        TB8RefMut {
            inner: &mut self.inner[mi_start..mi_end],
            m: self.m,
        }
    }
    fn set_val(&mut self, mi: usize, val: u8) {
        self.inner_mut()[mi] = val;
    }
}

// What I want for iter is
// 1. mut for &'a TB8 itself
// 2. ref for &TB8Ref<'a>
// into_iter() is unnecessary

// IterMut
// .col_iter_mut()
pub struct TB8IterMut<'a> {
    tb8refmut: &'a mut TB8,
    mi: usize,
}

// .iter_mut()
pub struct TB8ValIterMut<'a> {
    tb8refmut: &'a mut TB8,
    mi: usize,
}

// for &'a TB8
// Item=TB8Ref<'a>
// https://stackoverflow.com/questions/63437935/in-rust-how-do-i-create-a-mutable-iterator
// slice::split_at_mut
// https://stackoverflow.com/questions/60072498/how-to-implement-iterator-yielding-mutable-references
// std::mem::transmute
impl<'a> Iterator for TB8IterMut<'a> {
    type Item = TB8RefMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi + 2 <= self.tb8refmut.m {
            let mi = self.mi;
            self.mi += 2;

            unsafe {
                // recommended here
                // https://doc.rust-lang.org/stable/std/mem/fn.transmute.html
                // use like split_at_stdlib()
                //let s = &mut self.tb8refmut.inner[mi..mi + 2];
                // s should be slice
                let s = &mut self.tb8refmut.inner[mi..mi + 1];
                let ptr = s.as_mut_ptr();
                let p = std::slice::from_raw_parts_mut(ptr, 2);

                //let s = &mut self.tb8refmut.inner[mi..mi + 2];
                /*
                // how to avoid transmute?
                // -> above
                let p = std::mem::transmute::<&mut [u8], &mut [u8]>(
                    &mut self.tb8refmut.inner[mi..mi + 2],
                );
                 */
                /*
                // error for lifetime
                let (_, right) = self.tb8refmut.inner.as_slice().split_at_mut(mi);
                let (p, _) = right.split_at_mut(2);
                 */
                let r = TB8RefMut { inner: p, m: 2 };
                return Some(r);
            }
        }
        None
    }
}

/*
// ok
// for &'a TB8
// Item=&'a
// https://stackoverflow.com/questions/63437935/in-rust-how-do-i-create-a-mutable-iterator
// slice::split_at_mut
// https://stackoverflow.com/questions/60072498/how-to-implement-iterator-yielding-mutable-references
// std::mem::transmute
impl<'a> Iterator for TB8IterMut<'a> {
    // start from easy ver.
    type Item = &'a mut [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi + 2 < self.tb8refmut.m {
            let mi = self.mi;
            self.mi += 2;

            unsafe {
                return Some(std::mem::transmute(&mut self.tb8refmut.inner[mi..mi + 2]));
            }
        }
        None
    }
}
 */

/*
// ok
// for &'a TB8
// Item=u8
// https://stackoverflow.com/questions/63437935/in-rust-how-do-i-create-a-mutable-iterator
// slice::split_at_mut
// https://stackoverflow.com/questions/60072498/how-to-implement-iterator-yielding-mutable-references
// std::mem::transmute
impl<'a> Iterator for TB8IterMut<'a> {
    // start from easy ver.
    type Item = &'a mut u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi < self.tb8refmut.m {
            let mi = self.mi;
            self.mi += 1;

            unsafe {
                return Some(std::mem::transmute(&mut self.tb8refmut.inner[mi]));
            }
        }
        None
    }
}
 */

/*
pub struct TB8RefMutIntoIterU8<'a> {
    tb8refmut: TB8RefMut<'a>,
    mi: usize,
}

// ok for Item=u8
impl<'a> IntoIterator for TB8RefMut<'a> {
    // start from easy ver.
    type Item = u8;
    //type Item = TB8RefMut<'a>;
    type IntoIter = TB8RefMutIntoIterU8<'a>;
    fn into_iter(self) -> Self::IntoIter {
        TB8RefMutIntoIterU8 {
            tb8refmut: self,
            mi: 0,
        }
    }
}

impl<'a> Iterator for TB8RefMutIntoIterU8<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi < self.tb8refmut.m {
            let res = self.tb8refmut.access(self.mi);
            self.mi += 1;
            Some(res)
        } else {
            None
        }
    }
}
 */

/*
// for parallel
impl<'a> Iterator for TB8RefMutIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        let res = self.tb8refmut.access(self.mi);
        self.mi += 1;
        Some(res)
    }
}
 */

/*
// implement iter for struct
// https://stackoverflow.com/questions/30218886/how-to-implement-iterator-and-intoiterator-for-a-simple-struct
pub struct TB8RefMutIter<'a> {
    tb8refmut: TB8RefMut<'a>,
    mi_start: usize,
    mi_end: usize,
}

// iterator for Matrix
// https://athemathmo.github.io/rulinalg/doc/src/rulinalg/matrix/iter.rs.html#11-36


impl<'a> IntoIterator for TB8RefMut<'a> {
    // start from easy ver.
    //type Item = u8;
    type Item = TB8RefMut<'a>;
    type IntoIter = TB8RefMutIter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        TB8RefMutIter {
            tb8refmut: self,
            mi_start: 0,
            mi_end: 2,
        }
    }
}

// How to do...?
// could be done similary to ChunksMut but seems difficult...
impl<'a> Iterator for TB8RefMutIter<'a> {
    type Item = TB8RefMut<'a>;
    //type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi_end + 2 < self.tb8refmut.m {
            let res = self.tb8refmut.sub_mut3(self.mi_start, self.mi_end);
            self.mi_start = self.mi_end;
            self.mi_end = self.mi_start + 2;
            Some(res)
        } else {
            None
        }
    }
}
 */

impl<'a> GB8Ref<'a> {
    fn inner_mut(&mut self) -> &mut [B8] {
        self.0.inner_mut()
    }
    fn sub_mut(&mut self, mi_start: usize, mi_end: usize) -> GB8Ref {
        let m = self.0.m;
        GB8Ref(TB8RefMut {
            inner: &mut self.inner_mut()[mi_start..mi_end],
            //inner: &mut self.0.inner[..mi_end],
            m,
        })
    }
    fn set_val(&mut self, val: u8, mi: usize) {
        self.inner_mut()[mi] = val;
    }
}

pub fn test1() {
    let mut tb8 = TB8 {
        inner: vec![1, 2, 3, 4, 5, 6],
        m: 6,
    };

    // ok
    println!("mut bridge bfr {:?}", tb8.inner);
    tb8.iter_mut().par_bridge().for_each(|x| x.inner[0] = 100);
    println!("mut bridge afr {:?}", tb8.inner);

    // ok
    println!("mut bfr {:?}", tb8.inner);
    for x in tb8.iter_mut() {
        println!("inner {:?}", x.inner);
        if x.inner[0] == 1 {
            x.inner[0] = 100;
        }
    }
    println!("mut afr {:?}", tb8.inner);

    let mut tb8ref = TB8RefMut {
        inner: &mut tb8.inner,
        m: tb8.m,
    };
    //let mut genot_ref = GenotTwinRef(tb8);
    //let genot_ref_mut = &mut genot_ref;

    //println!("access genot {}", genot_ref.access(0));
    //println!("access genot {}", (&mut genot_ref).access(0));
    //println!("genot_ref {}", &genot_ref);
    //genot_ref.sub_mut(0, 2).inner_mut()[0] = 8;
    //println!("access2 genot {}", genot_ref.access(0));

    //let genot_ref_mut = &mut genot_ref;
    //genot_ref.set_val2(0, 1);
    tb8ref.sub_mut(0, 1);
    //genot_ref.set_val2(0, 1);
    tb8ref.sub_mut(0, 1);

    /*
    // into_iter()
    // ok for Item=u8
    println!("into_iter()");
    for x in tb8ref {
        println!("x {}", x);
    }
     */

    /*
    println!("iter()");
    for x in tb8ref.iter() {
        println!("x {}", x);
    }
     */

    /*
    for x in tb8ref {
        println!("y {}", x.access(0));
    }
     */
}

pub fn test2() {
    let mut genot = GB8(TB8 {
        inner: vec![1, 2, 3, 4, 5, 6],
        m: 6,
    });

    let mut genot_ref = GB8Ref(TB8RefMut {
        inner: &mut genot.0.inner,
        m: genot.0.m,
    });
    //let mut genot_ref = GenotTwinRef(tb8);
    //let genot_ref_mut = &mut genot_ref;

    //println!("access genot {}", genot_ref.access(0));
    //println!("access genot {}", (&mut genot_ref).access(0));
    //println!("genot_ref {}", &genot_ref);
    //genot_ref.sub_mut(0, 2).inner_mut()[0] = 8;
    //println!("access2 genot {}", genot_ref.access(0));

    //let genot_ref_mut = &mut genot_ref;
    //genot_ref.set_val2(0, 1);
    genot_ref.sub_mut(0, 1);
    //genot_ref.set_val2(0, 1);
    genot_ref.sub_mut(0, 1);

    [0, 2, 4].iter().for_each(|&mi| {
        //let mut genot_ref = &mut genot_ref;
        genot_ref.sub_mut(mi, mi + 2).set_val(8, 0);
    });

    // try par_mut
    // of cource, not able
    // should implement iterator by myself
    /*
    [0, 2, 4].par_iter().for_each(|&mi| {
        //let mut genot_ref = &mut genot_ref;
        genot_ref.sub_mut(mi, mi + 2).set_val(8, 0);
    });
     */

    // par_bridge? for my iterator?
    // -> NO, implement IntoParallelIterator
    // https://stackoverflow.com/questions/70637421/how-to-turn-a-string-into-a-parallel-iterator-using-rayon-and-rust
}

pub fn test() {
    test1();
    test2();
}
