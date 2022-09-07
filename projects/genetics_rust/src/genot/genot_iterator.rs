use super::genot_struct::{GenotSnvMut, GenotSnvRef};
use cmatrix::prelude::*;

pub struct GenotCountIter<'a>(ScalIter<'a>);

impl<'a> GenotCountIter<'a> {
    pub fn new(iter: ScalIter) -> GenotCountIter {
        GenotCountIter(iter)
    }
}

impl<'a> Iterator for GenotCountIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        return self.0.next();
    }
}

pub struct GenotIter<'a>(RowIter<'a>);

impl<'a> GenotIter<'a> {
    pub fn new(iter: RowIter) -> GenotIter {
        GenotIter(iter)
    }
}

impl<'a> Iterator for GenotIter<'a> {
    type Item = GenotSnvRef<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(x) = self.0.next() {
            // GenotRef::new(x)
            return Some(GenotSnvRef::from_cvecref(x));
        }
        None
    }
}

pub struct GenotIterMut<'a>(RowIterMut<'a>);

impl<'a> GenotIterMut<'a> {
    pub fn new(iter: RowIterMut) -> GenotIterMut {
        GenotIterMut(iter)
    }
}

impl<'a> Iterator for GenotIterMut<'a> {
    type Item = GenotSnvMut<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(x) = self.0.next() {
            // GenotRef::new(x)
            return Some(GenotSnvMut::new(x));
        }
        None
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
    use crate::genot::genot_struct::Genot;
    use crate::genot::{BaseGenot, BaseGenotMut, BaseGenotSnvMut};
    use rayon::prelude::*;

    #[test]
    fn test_iter_mut() {
        let vec = vec![1, 2, 3, 0, 1, 2];
        let mut m = Genot::new(2, 3, vec);

        m.iter_snv_mut().for_each(|mut row| row.set(3, 0));

        assert_eq!(m.vals_snv(0), vec![3, 2, 3]);
        assert_eq!(m.vals_snv(1), vec![3, 1, 2]);
    }

    #[test]
    fn test_iter_mut_par() {
        let vec = vec![1, 2, 3, 0, 1, 2];
        let mut m = Genot::new(2, 3, vec);

        m.iter_snv_mut()
            .par_bridge()
            .for_each(|mut row| row.set(3, 0));

        assert_eq!(m.vals_snv(0), vec![3, 2, 3]);
        assert_eq!(m.vals_snv(1), vec![3, 1, 2]);
    }
}
