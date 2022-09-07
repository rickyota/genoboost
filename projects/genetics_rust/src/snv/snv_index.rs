use std::str::FromStr;

use super::Chrom;

type Alleles = (String, String);

#[derive(Clone, Hash, Debug, Default)]
pub struct SnvId {
    rs: String,
    chrom: Chrom,
    pos: usize,
    alleles: Alleles,
    // only for one letter
    alleles_revcomp: Alleles,
    sida: String,
}

// how to implement?
// 1. hash["A"]="T" -> cannot create const hash
// 2. match "A" => "T"
/// assume str is A,C,G,T
fn complement_allele(a: &str) -> String {
    match a {
        "A" => "T".to_owned(),
        "C" => "G".to_owned(),
        "G" => "C".to_owned(),
        "T" => "A".to_owned(),
        _ => panic!("allele is not one of A,C,G,T"),
    }
}

impl SnvId {
    // TODO: Should change to
    pub fn construct_snv_index_string(
        rs: String,
        chrom: &str,
        pos: &str,
        a1: String,
        a2: String,
    ) -> SnvId {
        let mut snv = SnvId {
            rs,
            chrom: Chrom::from_str(&chrom).unwrap(),
            pos: pos.parse::<usize>().unwrap(),
            alleles: (a1, a2),
            alleles_revcomp: ("".to_owned(), "".to_owned()),
            sida: "".to_string(),
        };
        snv.check_alleles();
        snv.set_alleles_revcomp();
        snv.set_sida();
        snv
    }

    // TODO: better way
    // TODO: add N etc.
    /// should consist of A,C,G,T
    fn check_alleles(&self) {
        fn check_allele(a: &str) {
            a.chars()
                .all(|v| (v == 'A') || (v == 'C') || (v == 'G') || (v == 'T'));
        }
        check_allele(self.a1());
        check_allele(self.a2());
    }

    fn set_alleles_revcomp(&mut self) {
        if self.is_one_letter() {
            self.alleles_revcomp = (complement_allele(self.a2()), complement_allele(self.a1()))
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

    pub fn alleles(&self) -> (&str, &str) {
        (&self.alleles.0, &self.alleles.1)
    }
    pub fn alleles_revcomp(&self) -> (&str, &str) {
        (&self.alleles_revcomp.0, &self.alleles_revcomp.1)
    }

    pub fn a1(&self) -> &str {
        &self.alleles.0
    }
    pub fn a2(&self) -> &str {
        &self.alleles.1
    }

    pub fn is_one_letter(&self) -> bool {
        (self.a1().len() == 1) && (self.a2().len() == 1)
    }

    // list all candidates
    pub fn flip_or_complement(
        &self,
    ) -> Option<((&str, &str), (&str, &str), (&str, &str), (&str, &str))> {
        if !self.is_one_letter() {
            return None;
        }
        let a = self.alleles();
        let a_flip = (a.1, a.0);
        let a_revcomp = self.alleles_revcomp();
        let a_revcomp_flip = (a_revcomp.1, a_revcomp.0);

        Some((a, a_flip, a_revcomp, a_revcomp_flip))
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
        match self.flip_or_complement() {
            None => self.alleles() == a_other,
            Some((a, a_flip, a_revcomp, a_revcomp_flip)) => {
                (a == a_other)
                    || (a_flip == a_other)
                    || (a_revcomp == a_other)
                    || (a_revcomp_flip == a_other)
            }
        }
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
        let snv_index = SnvId::construct_snv_index_string(
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

        let alleles = snv_index.flip_or_complement().unwrap();
        assert_eq!(alleles.0, ("A", "C"));
        assert_eq!(alleles.1, ("C", "A"));
        assert_eq!(alleles.2, ("G", "T"));
        assert_eq!(alleles.3, ("T", "G"));
    }

    #[test]
    #[should_panic]
    fn test_construct_snv_index_string_panic() {
        let _ = SnvId::construct_snv_index_string(
            "rs1".to_owned(),
            "1",
            "123",
            "N".to_owned(),
            "N".to_owned(),
        );
    }

    /// test of Display
    #[test]
    fn test_snv_index_to_string() {
        let snv_index = SnvId::construct_snv_index_string(
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
        let snv_index = SnvId::construct_snv_index_string(
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
        let snv_index1 = SnvId::construct_snv_index_string(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_index2 = SnvId::construct_snv_index_string(
            "rs1".to_owned(),
            "1",
            "123",
            "G".to_owned(),
            "T".to_owned(),
        );

        assert_eq!(snv_index1, snv_index2);

        let snv_index3 = SnvId::construct_snv_index_string(
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
        let snv_index1 = SnvId::construct_snv_index_string(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "C".to_owned(),
        );
        let snv_index2 = SnvId::construct_snv_index_string(
            "rs1".to_owned(),
            "1",
            "123",
            "A".to_owned(),
            "G".to_owned(),
        );

        assert!(snv_index1 < snv_index2);

        let snv_index3 = SnvId::construct_snv_index_string(
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
