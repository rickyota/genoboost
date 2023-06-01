//! Iterator
//!
//!
//!
//! best
//! https://aloso.github.io/2021/03/09/creating-an-iterator
//!
//! https://stackoverflow.com/questions/30218886/how-to-implement-iterator-and-intoiterator-for-a-simple-struct
//! matrix
//! https://athemathmo.github.io/rulinalg/doc/src/rulinalg/matrix/iter.rs.html#11-36
//!

struct Genot {
    inner: Vec<u8>,
    m: usize,
}

struct GenotIter<'a> {
    genot: &'a Genot,
    mi: usize,
}

struct GenotIntoIter {
    genot: Genot,
    mi: usize,
}

impl Genot {
    fn iter(&self) -> GenotIter {
        GenotIter { genot: self, mi: 0 }
    }
}

// IntoIterator: to get Iterator
//
// for x in &g
// for x in g.into_iter()
impl<'a> IntoIterator for &'a Genot {
    type Item = u8;
    type IntoIter = GenotIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
        //GenotIntoIter { genot: self, mi: 0 }
    }
}

impl<'a> Iterator for GenotIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi < self.genot.m {
            let mi = self.mi;
            self.mi += 1;
            return Some(self.genot.inner[mi]);
        }
        None
    }
}

/*
// without ref
// for x in g
// for x in g.into_iter() // both in upper?
impl IntoIterator for Genot {
    type Item = u8;
    type IntoIter = GenotIntoIter;

    fn into_iter(self) -> Self::IntoIter {
        GenotIntoIter { genot: self, mi: 0 }
        //self.iter()
    }
}

impl Iterator for GenotIntoIter {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi < self.genot.m {
            let mi = self.mi;
            self.mi += 1;
            return Some(self.genot.inner[mi]);
        }
        None
    }
}
 */

/*
struct GenotIntoIter {
    genot: Genot,
    mi: usize,
}

impl Iterator for GenotIntoIter {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        if self.mi < self.genot.m {
            let mi = self.mi;
            self.mi += 1;
            return Some(self.genot.inner[mi]);
        }
        None
    }
}
 */

/*
struct GenotSliceIter<'a> {
    slice: &'a [u8],
    mi_start: usize,
    mi_end: usize,
    n: usize,
}

impl<'a> Iterator for GenotSliceIter<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.mi_start + 2 < self.n {
            let s = &self.slice[self.mi_start..self.mi_end];
            self.mi_start = self.mi_end;
            self.mi_end = self.mi_end + 2;
            return Some(s);
        }
        None
    }
}
*/

pub fn test() {
    let a = vec![1, 2, 3, 4, 5, 6];
    let g = Genot { inner: a, m: 6 };

    for x in g.iter() {
        println!("x {}", x);
    }

    // impl<'a> IntoIterator for &'a Genot
    //for x in g.into_iter() {
    for x in &g {
        println!("x {}", x);
    }

    let a = vec![1, 2, 3, 4, 5, 6];
    let g = Genot { inner: a, m: 6 };

    // impl IntoIterator for Genot
    for x in g.into_iter() {
        println!("x {}", x);
    }

    //println!("g {}", g);
}
