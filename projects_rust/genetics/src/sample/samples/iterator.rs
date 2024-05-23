//use cmatrix::prelude::*;
//
//pub struct PheIter<'a>(BoolIter<'a>);
//
//impl<'a> PheIter<'a> {
//    pub fn new(iter: BoolIter) -> PheIter {
//        PheIter(iter)
//    }
//}
//
//impl<'a> Iterator for PheIter<'a> {
//    type Item = bool;
//    fn next(&mut self) -> Option<Self::Item> {
//        return self.0.next();
//    }
//}
//
