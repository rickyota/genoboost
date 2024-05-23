#[derive(PartialEq, Debug, Clone, Copy)]
pub enum LdCriteria {
    Radius(usize),
    R2(f64),
}

impl LdCriteria {
    pub fn new_radius(radius: usize) -> Self {
        Self::Radius(radius)
    }

    pub fn new_r2(r2: f64) -> Self {
        if r2 < 0.0 || r2 > 1.0 {
            panic!("r2 must be between 0 and 1.");
        }
        Self::R2(r2)
    }

    pub fn new(r2: Option<f64>, radius: Option<usize>) -> Self {
        match (r2, radius) {
            (Some(r2), None) => Self::new_r2(r2),
            (None, Some(radius)) => Self::new_radius(radius),
            (None, None) => panic!("Cannot specify neither r2 nor radius."),
            (Some(_), Some(_)) => panic!("Cannot specify both r2 and radius."),
        }
    }
}
