struct Atom {
    symbol: Option<String>,
    atomic_number: Option<u32>,
    mass: Option<f64>,
    coordinates: (f64, f64, f64),
}

impl Atom {
    fn new(symbol: Option<String>, atomic_number: Option<u32>, mass: Option<f64>, coordinates: (f64, f64, f64)) -> Self {
        Atom {
            symbol,
            atomic_number,
            mass,
            coordinates,
        }
    }

    pub fn distance(&self, other: &Atom) -> f64 {
        let dx = self.coordinates.0 - other.coordinates.0;
        let dy = self.coordinates.1 - other.coordinates.1;
        let dz = self.coordinates.2 - other.coordinates.2;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    pub fn angle(&self, other1: &Atom, other2: &Atom) -> f64 {
        let d1 = (
            other1.coordinates.0 - self.coordinates.0,
            other1.coordinates.1 - self.coordinates.1,
            other1.coordinates.2 - self.coordinates.2,
        );
        let d2 = (
            other2.coordinates.0 - self.coordinates.0,
            other2.coordinates.1 - self.coordinates.1,
            other2.coordinates.2 - self.coordinates.2,
        );
        let dot_product = d1.0 * d2.0 + d1.1 * d2.1 + d1.2 * d2.2;
        let dist1 = Atom::distance(self, other1);
        let dist2 = Atom::distance(self, other2);
        (dot_product / (dist1 * dist2)).acos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 0.0, 0.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (1.0, 1.0, 1.0));
        assert_eq!(Atom::distance(&atom1, &atom2), (3.0_f64).sqrt());
    }

    #[test]
    fn test_angle() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 0.0, 0.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (1.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 1.0, 0.0));
        assert_eq!(Atom::angle(&atom1, &atom2, &atom3), std::f64::consts::FRAC_PI_2);
    }
}