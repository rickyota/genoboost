use std::convert::TryFrom;
use std::str::FromStr;

#[derive(Eq, PartialEq, Ord, PartialOrd, Clone, Hash, Debug)]
pub enum Chrom {
    Auto(usize), // 1-22, or more if sex is converted into usize
    Sex(String), // X,Y,XY
    Mt,          // MT
}

impl Chrom {
    // 30 is for when sex is iconverted into usize
    // auto: 1-22
    // call Chrom::auto_max()
    pub fn auto_max() -> usize {
        30
    }

    // TODO: create iterator
    pub fn variants() -> Vec<Chrom> {
        let vec: Vec<usize> = (1..=Chrom::auto_max()).collect();
        let mut vec: Vec<Chrom> = vec.iter().map(|&v| Chrom::try_from(v).unwrap()).collect();

        let vec_ = ["X", "Y", "XY"];
        let vec_: Vec<Chrom> = vec_.iter().map(|&v| Chrom::from_str(v).unwrap()).collect();
        vec.extend(vec_);

        let vec_: Vec<Chrom> = [Chrom::Mt].to_vec();
        vec.extend(vec_);

        vec
    }
}

impl Default for Chrom {
    fn default() -> Self {
        Chrom::Auto(1)
    }
}

// This auto-implement to_string()
impl std::fmt::Display for Chrom {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let chrom_string = match self {
            Chrom::Auto(chrom) => chrom.to_string(),
            Chrom::Sex(chrom) => chrom.to_owned(),
            Chrom::Mt => "MT".to_owned(),
        };
        write!(f, "{}", chrom_string.as_str())
    }
}

impl FromStr for Chrom {
    // FromStr shoudl be `Err`
    type Err = String;
    // should be consistent to Chrom::variants()
    fn from_str(str: &str) -> Result<Self, Self::Err> {
        // if able to parse into usize
        if let Ok(chrom) = str.parse::<usize>() {
            if chrom >= 1 && chrom <= Chrom::auto_max() {
                return Ok(Chrom::Auto(chrom));
            } else {
                return Err(format!(
                    "Autosomal chromosome should be bet. 1 and {}
                 not {}",
                    Chrom::auto_max(),
                    str
                ));
            }
        }

        // string
        if str == "MT" {
            return Ok(Chrom::Mt);
        } else if (str == "X") || (str == "Y") || (str == "XY") {
            return Ok(Chrom::Sex(str.to_owned()));
        } else {
            return Err(format!(
                "Chrom string should be one of X,Y,XY,MT not {}",
                str
            ));
        }
    }
}

// only admit 1 - Chrom::auto_max()
impl TryFrom<usize> for Chrom {
    // TryFrom should be `Error`
    type Error = String;

    //type Err = String;
    fn try_from(val: usize) -> Result<Self, Self::Error> {
        if val >= 1 && val <= Chrom::auto_max() {
            Ok(Chrom::Auto(val))
        } else {
            return Err(format!(
                "Autosomal chromosome should be bet. 1 and {}
                 not {}",
                Chrom::auto_max(),
                val
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::convert::{TryFrom, TryInto};

    #[test]
    fn test_chrom_variants() {
        let vec = Chrom::variants();

        let vec_exp = [
            Chrom::Auto(1),
            Chrom::Auto(2),
            Chrom::Auto(3),
            Chrom::Auto(4),
            Chrom::Auto(5),
            Chrom::Auto(6),
            Chrom::Auto(7),
            Chrom::Auto(8),
            Chrom::Auto(9),
            Chrom::Auto(10),
            Chrom::Auto(11),
            Chrom::Auto(12),
            Chrom::Auto(13),
            Chrom::Auto(14),
            Chrom::Auto(15),
            Chrom::Auto(16),
            Chrom::Auto(17),
            Chrom::Auto(18),
            Chrom::Auto(19),
            Chrom::Auto(20),
            Chrom::Auto(21),
            Chrom::Auto(22),
            Chrom::Auto(23),
            Chrom::Auto(24),
            Chrom::Auto(25),
            Chrom::Auto(26),
            Chrom::Auto(27),
            Chrom::Auto(28),
            Chrom::Auto(29),
            Chrom::Auto(30),
            Chrom::Sex("X".to_owned()),
            Chrom::Sex("Y".to_owned()),
            Chrom::Sex("XY".to_owned()),
            Chrom::Mt,
        ]
        .to_vec();

        assert_eq!(vec, vec_exp);
    }

    /// test of TryFrom
    #[test]
    fn test_chrom_from_str() {
        let chrom = "10";
        assert_eq!(Chrom::from_str(chrom).unwrap(), Chrom::Auto(10));

        let chrom = "X";
        assert_eq!(Chrom::from_str(chrom).unwrap(), Chrom::Sex("X".to_owned()));

        let chrom = "MT";
        assert_eq!(Chrom::from_str(chrom).unwrap(), Chrom::Mt);
    }

    /*
    /// test of From
    #[test]
    fn test_chrom_from_usize() {
        let chrom: usize = 10;
        assert_eq!(Chrom::from(chrom), Chrom::Auto(10));
    }
    */

    /// test of TryFrom
    #[test]
    fn test_chrom_try_from_usize() {
        let chrom: usize = 10;
        assert_eq!(Chrom::try_from(chrom).unwrap(), Chrom::Auto(10));
    }

    /// test of TryFrom
    #[test]
    #[should_panic]
    fn test_chrom_try_from_usize_panic() {
        let chrom: usize = 31;
        Chrom::try_from(chrom).unwrap();
    }

    #[test]
    fn test_chrom_usize_try_into() {
        // Chrom::from() is more stable
        let chrom: usize = 10;
        let chrom_: Chrom = chrom.try_into().unwrap();
        assert_eq!(chrom_, Chrom::Auto(10));

        let chrom: usize = 10;
        assert_eq!(Chrom::Auto(10), chrom.try_into().unwrap());

        // looks impossible https://stackoverflow.com/questions/41207885/using-generic-trait-methods-like-into-when-type-inference-is-impossible
        //let chrom: usize = 10;
        //assert_eq!(chrom.try_into().unwrap(), Chrom::Auto(10)); // requires type annotation
        //assert_eq!(chrom.try_into::<Chrom>().unwrap(), Chrom::Auto(10)); // requires type annotation
    }

    /// test of Display
    #[test]
    fn test_chrom_to_string() {
        let chrom = Chrom::Auto(10);
        assert_eq!(chrom.to_string(), "10".to_owned());

        let chrom = Chrom::Sex("X".to_owned());
        assert_eq!(chrom.to_string(), "X".to_owned());

        let chrom = Chrom::Mt;
        assert_eq!(chrom.to_string(), "MT".to_owned());
    }
}
