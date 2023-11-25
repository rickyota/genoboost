use std::{convert::TryFrom, str::FromStr};

use super::Chrom;

type Alleles = (String, String);

// TODO: remove allele_flip and generate when needed?
/// reverse: A1 <-> A2
/// flip : in A1, A2, A <-> T
/// TODO: register ref/alt and minor/major alleles
/// all programs uses minor/major but input plink might have different one
#[derive(Clone, Hash, Debug, Default)]
pub struct SnvId {
    // TODO: rs->id
    rs: String,
    chrom: Chrom,
    pos: usize,
    alleles: Alleles,
    // only for one letter
    alleles_flip: Alleles,
    //alleles_ref_alt: Alleles,
    sida: String,
    sid: String,
    //vid: String,
}

// how to implement?
// 1. hash["A"]="T" -> cannot create const hash
// 2. match "A" => "T"
/// assume str is A,C,G,T
fn flip_allele(a: &str) -> String {
    match a {
        "A" => "T".to_owned(),
        "C" => "G".to_owned(),
        "G" => "C".to_owned(),
        "T" => "A".to_owned(),
        _ => panic!("allele is not one of A,C,G,T"),
    }
}

impl SnvId {
    // a1, a2 cannot include other than "ACGT".
    pub fn construct_snv_index(
        rs: String,
        chrom: &str,
        pos: &str,
        a1: String,
        a2: String,
        // this is ok but need changes all over the codes; instead use in vid()
        //use_snv_pos: bool,
    ) -> SnvId {
        let mut snv = SnvId {
            rs,
            chrom: Chrom::from_str(&chrom).unwrap(),
            pos: pos.parse::<usize>().unwrap(),
            alleles: (a1, a2),
            alleles_flip: ("".to_owned(), "".to_owned()),
            sida: "".to_string(),
            sid: "".to_string(),
            //vid: "".to_string(),
        };
        //snv.check_alleles();
        snv.set_alleles_flip();
        snv.set_sida();
        snv.set_sid();
        //snv.set_vid(use_snv_pos);
        snv
    }

    // for use_snvs
    pub fn construct_snv_index_rs(rs: String) -> SnvId {
        let snv = SnvId {
            rs,
            // dummy for chrom1
            chrom: Chrom::try_from(1).unwrap(),
            pos: 0,
            alleles: ("".to_owned(), "".to_owned()),
            alleles_flip: ("".to_owned(), "".to_owned()),
            sida: "".to_string(),
            sid: "".to_string(),
        };
        //snv.check_alleles();
        //snv.set_alleles_revcomp();
        //snv.set_sida();
        snv
    }

    // TODO: better way
    // TODO: add N etc.
    // TODO: 1kg contained '<DEL>'
    /// should consist of A,C,G,T
    #[allow(dead_code)]
    fn check_alleles(&self) {
        fn check_allele(a: &str) {
            if !(a
                .chars()
                .all(|v| (v == 'A') || (v == 'C') || (v == 'G') || (v == 'T')))
            {
                panic!("Alleles should be one of A, C, G, or T: {}.", a);
            }
        }
        check_allele(self.a1());
        check_allele(self.a2());
    }

    fn set_alleles_flip(&mut self) {
        if self.is_one_letter() {
            self.alleles_flip = (flip_allele(self.a1()), flip_allele(self.a2()))
            //self.alleles_rev = (complement_allele(self.a2()), complement_allele(self.a1()))
        }
        // else stay ("","")
    }

    fn set_sida(&mut self) {
        self.sida = self.chrom.to_string()
            + ":"
            + &self.pos.to_string()
            + ":"
            + &self.a1()
            + ":"
            + &self.a2();
    }

    fn set_sid(&mut self) {
        self.sid = self.chrom.to_string() + ":" + &self.pos.to_string();
    }

    pub fn reverse_alleles(&mut self) {
        // update alleles, alleles_flip, sida
        self.alleles = (self.alleles.1.clone(), self.alleles.0.clone());
        self.set_alleles_flip();
        self.set_sida();
    }

    pub fn rs(&self) -> &str {
        &self.rs
    }

    pub fn chrom(&self) -> &Chrom {
        &self.chrom
    }

    pub fn pos(&self) -> usize {
        self.pos
    }

    pub fn sida(&self) -> &str {
        &self.sida
    }

    fn _sid(&self) -> &str {
        &self.sid
    }

    pub fn vid(&self, use_snv_pos: bool) -> &str {
        if use_snv_pos {
            &self.sid
        } else {
            &self.rs
        }
    }

    fn alleles(&self) -> (&str, &str) {
        (&self.alleles.0, &self.alleles.1)
    }

    fn alleles_rev(&self) -> (&str, &str) {
        (&self.alleles.1, &self.alleles.0)
    }

    fn alleles_flip(&self) -> (&str, &str) {
        (&self.alleles_flip.0, &self.alleles_flip.1)
    }

    fn alleles_rev_flip(&self) -> (&str, &str) {
        (&self.alleles_flip.1, &self.alleles_flip.0)
    }

    pub fn a1(&self) -> &str {
        &self.alleles.0
    }
    pub fn a2(&self) -> &str {
        &self.alleles.1
    }
    // flip
    fn _a1f(&self) -> &str {
        &self.alleles_flip.0
    }
    // flip
    fn _a2f(&self) -> &str {
        &self.alleles_flip.1
    }

    //pub fn to_sid(&self) -> String {
    //    self.chrom.to_string() + ":" + &self.pos.to_string()
    //}

    //pub fn to_sida_rev(&self) -> String {
    //    self.chrom.to_string() + ":" + &self.pos.to_string() + ":" + &self.a2() + ":" + &self.a1()
    //}
    //// Only for one letter
    //pub fn to_sida_flip(&self) -> Option<String> {
    //    if self.is_one_letter() {
    //        Some(
    //            self.chrom.to_string()
    //                + ":"
    //                + &self.pos.to_string()
    //                + ":"
    //                + &self.a1f()
    //                + ":"
    //                + &self.a2f(),
    //        )
    //    } else {
    //        None
    //    }
    //}
    //// Only for one letter
    //pub fn to_sida_rev_flip(&self) -> Option<String> {
    //    if self.is_one_letter() {
    //        Some(
    //            self.chrom.to_string()
    //                + ":"
    //                + &self.pos.to_string()
    //                + ":"
    //                + &self.a2f()
    //                + ":"
    //                + &self.a1f(),
    //        )
    //    } else {
    //        None
    //    }
    //}

    fn is_one_letter(&self) -> bool {
        (self.a1().len() == 1) && (self.a2().len() == 1)
    }

    // list all candidates
    // TODO: renmae
    fn flip_or_rev(
        &self,
    ) -> (
        (&str, &str),
        (&str, &str),
        Option<((&str, &str), (&str, &str))>,
    ) {
        let a = self.alleles();
        let a_rev = (a.1, a.0);

        if !self.is_one_letter() {
            return (a, a_rev, None);
        }
        let a_flip = self.alleles_flip();
        let a_rev_flip = (a_flip.1, a_flip.0);

        (a, a_rev, Some((a_flip, a_rev_flip)))
    }

    // to reverse genotype
    pub fn is_rev(&self, snv: &SnvId, use_snv_pos: bool) -> bool {
        //if self.sid() != snv.sid() {
        if self.vid(use_snv_pos) != snv.vid(use_snv_pos) {
            return false;
        }
        if self.is_one_letter() {
            return (self.alleles() == snv.alleles_rev())
                | (self.alleles() == snv.alleles_rev_flip());
            // otherwise not rev or alleles do not match
        } else {
            return self.alleles() == snv.alleles_rev();
        }
    }
}

/*
impl Default for SnvIndex {
    fn default() -> Self {
        Self {
            rs: "".to_owned(),
            chrom: Chrom::Auto(1),
            pos: 0,
            alleles: ("".to_owned(), "".to_owned()),
            alleles_revcomp: ("".to_owned(), "".to_owned()),
            sida: "".to_owned(),
        }
    }
}
*/

// This auto-implement to_string()
impl std::fmt::Display for SnvId {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.sida())
    }
}

impl PartialEq for SnvId {
    fn eq(&self, other: &Self) -> bool {
        if self.chrom() != other.chrom() || self.pos() != other.pos() {
            return false;
        }
        // same pos
        let a_other = other.alleles();
        let (a, a_flip, alleles_rev) = self.flip_or_rev();

        if (a == a_other) || (a_flip == a_other) {
            return true;
        }

        // check rev for one letter
        if let Some((a_rev, a_rev_flip)) = alleles_rev {
            if (a_rev == a_other) || (a_rev_flip == a_other) {
                return true;
            }
        }
        return false;
    }
}

impl Eq for SnvId {}

// partial order is the same as order
impl PartialOrd for SnvId {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SnvId {
    // TODO: introducing self.eq is ok?
    // result of sorting would not be unique?
    // but not using eq is also strange
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // if eq, return Order.Equal first
        if self.eq(other) {
            return std::cmp::Ordering::Equal;
        }

        let ord = self.chrom().cmp(other.chrom());
        if ord.is_ne() {
            return ord;
        }
        let ord = self.pos().cmp(&other.pos());
        if ord.is_ne() {
            return ord;
        }
        let ord = self.a1().cmp(other.a1());
        if ord.is_ne() {
            return ord;
        }
        self.a2().cmp(other.a2())
    }
}

impl AsRef<SnvId> for SnvId {
    #[inline]
    fn as_ref(&self) -> &SnvId {
        self
    }
}

impl AsMut<SnvId> for SnvId {
    #[inline]
    fn as_mut(&mut self) -> &mut SnvId {
        self
    }
}

// later implement FromStr, TryFrom
// from 1:123:A:C or 1_123_A_C
// but ambiguous on A1 and A2

#[cfg(test)]
mod tests {
    use super::*;

    /// test of TryFrom
    #[test]
    fn test_construct_snv_index_string() {
        let snv_index = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        assert_eq!(snv_index.rs(), "rs1");
        assert_eq!(snv_index.chrom(), &Chrom::Auto(1));
        assert_eq!(snv_index.pos(), 123);
        assert_eq!(snv_index.a1(), "A");
        assert_eq!(snv_index.a2(), "C");
    }

    #[test]
    #[should_panic]
    fn test_construct_snv_index_string_panic() {
        let _ = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "N".to_owned(),
            "N".to_owned(),
        );
    }

    #[test]
    fn test_rev_flip() {
        let snv_index = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        assert_eq!(snv_index.alleles_rev(), ("C", "A"));
        assert_eq!(snv_index.alleles_flip(), ("T", "G"));
        assert_eq!(snv_index.alleles_rev_flip(), ("G", "T"));
    }

    #[test]
    fn test_rev_flip_long() {
        let snv_index = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "AAT".to_owned(),
            "C".to_owned(),
        );
        assert_eq!(snv_index.alleles_rev(), ("C", "AAT"));
        assert_eq!(snv_index.alleles_flip(), ("", ""));
        assert_eq!(snv_index.alleles_rev_flip(), ("", ""));
    }

    #[test]
    fn test_is_rev() {
        let snv_index1 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_index2 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        assert!(!snv_index1.is_rev(&snv_index2, false));

        let snv_index2 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "T".to_owned(),
            "G".to_owned(),
        );
        assert!(!snv_index1.is_rev(&snv_index2, false));

        let snv_index2 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "C".to_owned(),
            "A".to_owned(),
        );
        assert!(snv_index1.is_rev(&snv_index2, false));

        let snv_index2 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "G".to_owned(),
            "T".to_owned(),
        );
        assert!(snv_index1.is_rev(&snv_index2, false));

        // alleles do not match
        let snv_index2 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "G".to_owned(),
        );
        assert!(!snv_index1.is_rev(&snv_index2, false));

        // snv does not match
        let snv_index2 = SnvId::construct_snv_index(
            "rs2".to_owned(),
            "1",
            "124",
            "C".to_owned(),
            "A".to_owned(),
        );
        assert!(!snv_index1.is_rev(&snv_index2, false));
    }

    /// test of Display
    #[test]
    fn test_snv_index_to_string() {
        let snv_index = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        assert_eq!(snv_index.to_string(), "1:123:A:C".to_owned());
    }

    #[test]
    fn test_check_alleles() {
        let snv_index = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        snv_index.check_alleles()
    }

    #[test]
    fn test_eq() {
        let snv_index1 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_index2 = SnvId::construct_snv_index(
            //"rs1".to_owned(),
            "1:123:G:T".to_owned(),
            "1",
            "123",
            "G".to_owned(),
            "T".to_owned(),
        );

        assert_eq!(snv_index1, snv_index2);

        let snv_index3 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "G".to_owned(),
        );

        assert_ne!(snv_index1, snv_index3);
    }

    #[test]
    fn test_ord() {
        let snv_index1 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_index2 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "G".to_owned(),
        );

        assert!(snv_index1 < snv_index2);

        let snv_index3 = SnvId::construct_snv_index(
            "rs1".to_owned(),
            "1",
            "123",
            "G".to_owned(),
            "T".to_owned(),
        );
        // should be ==
        assert!(!(snv_index1 < snv_index3));
        assert!((snv_index1 <= snv_index3));
        assert!(!(snv_index1 > snv_index3));
    }
}
