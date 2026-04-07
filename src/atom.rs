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

    fn from_coords(coordinates: (f64, f64, f64)) -> Self {
        Atom {
            symbol: None,
            atomic_number: None,
            mass: None,
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
        // Calculate the angle between three atoms. Atom1-Atom2-Atom3, where Atom2 is the vertex. Return the angle in radians in f64.
        let d1 = (
            self.coordinates.0 - other1.coordinates.0,
            self.coordinates.1 - other1.coordinates.1,
            self.coordinates.2 - other1.coordinates.2,
        );
        let d2 = (
            other2.coordinates.0 - other1.coordinates.0,
            other2.coordinates.1 - other1.coordinates.1,
            other2.coordinates.2 - other1.coordinates.2,
        );
        let dot_product = d1.0 * d2.0 + d1.1 * d2.1 + d1.2 * d2.2;
        let dist1 = Atom::distance(self, other1);
        let dist2 = Atom::distance(other1, other2);
        ((dot_product / (dist1 * dist2)).clamp(-1.0, 1.0)).acos()
    }

    pub fn out_of_plane_angle(&self, other1: &Atom, other2: &Atom, other3: &Atom) -> f64 {
        let dkj = (
            other2.coordinates.0 - other1.coordinates.0,
            other2.coordinates.1 - other1.coordinates.1,
            other2.coordinates.2 - other1.coordinates.2,
        );
        let dkl = (
            other3.coordinates.0 - other2.coordinates.0,
            other3.coordinates.1 - other2.coordinates.1,
            other3.coordinates.2 - other2.coordinates.2,
        );
        let dij = (
            self.coordinates.0 - other1.coordinates.0,
            self.coordinates.1 - other1.coordinates.1,
            self.coordinates.2 - other1.coordinates.2,
        );

        let cross_product = (
            dkj.1 * dkl.2 - dkj.2 * dkl.1,
            dkj.2 * dkl.0 - dkj.0 * dkl.2,
            dkj.0 * dkl.1 - dkj.1 * dkl.0,
        );

        let norm_cross = (cross_product.0 * cross_product.0 + cross_product.1 * cross_product.1 + cross_product.2 * cross_product.2).sqrt();
        if norm_cross < 1e-10 {
            eprintln!("Warning: Atoms are collinear, out-of-plane angle is undefined. Returning 0.0.");
            return 0.0; // Avoid division by zero for collinear atoms
        }

        let dot_product = cross_product.0 * dij.0 + cross_product.1 * dij.1 + cross_product.2 * dij.2;
        let sin_angle = (dot_product / (norm_cross * (dij.0 * dij.0 + dij.1 * dij.1 + dij.2 * dij.2).sqrt())).clamp(-1.0, 1.0);
        sin_angle.asin()
    }

    pub fn dihedral_angle(&self, other1: &Atom, other2: &Atom, other3: &Atom) -> f64 {
        let dij = (
            other1.coordinates.0 - self.coordinates.0,
            other1.coordinates.1 - self.coordinates.1,
            other1.coordinates.2 - self.coordinates.2,
        );
        let djk = (
            other2.coordinates.0 - other1.coordinates.0,
            other2.coordinates.1 - other1.coordinates.1,
            other2.coordinates.2 - other1.coordinates.2,
        );
        let dkl = (
            other3.coordinates.0 - other2.coordinates.0,
            other3.coordinates.1 - other2.coordinates.1,
            other3.coordinates.2 - other2.coordinates.2,
        );
        let cross1 = (
            dij.1 * djk.2 - dij.2 * djk.1,
            dij.2 * djk.0 - dij.0 * djk.2,
            dij.0 * djk.1 - dij.1 * djk.0,
        );
        let cross2 = (
            djk.1 * dkl.2 - djk.2 * dkl.1,
            djk.2 * dkl.0 - djk.0 * dkl.2,
            djk.0 * dkl.1 - djk.1 * dkl.0,
        );
        let norm_cross1 = (cross1.0 * cross1.0 + cross1.1 * cross1.1 + cross1.2 * cross1.2).sqrt();
        let norm_cross2 = (cross2.0 * cross2.0 + cross2.1 * cross2.1 + cross2.2 * cross2.2).sqrt();
        if norm_cross1 < 1e-10 || norm_cross2 < 1e-10 {
            eprintln!("Warning: Atoms are collinear, dihedral angle is undefined. Returning 0.0.");
            return 0.0; // Avoid division by zero for collinear atoms
        }
        let cos_angle = (cross1.0 * cross2.0 + cross1.1 * cross2.1 + cross1.2 * cross2.2) / (norm_cross1 * norm_cross2);
        let cross_product = (
            cross1.1 * cross2.2 - cross1.2 * cross2.1,
            cross1.2 * cross2.0 - cross1.0 * cross2.2,
            cross1.0 * cross2.1 - cross1.1 * cross2.0,
        );
        let djk_norm = (djk.0 * djk.0 + djk.1 * djk.1 + djk.2 * djk.2).sqrt();
        let sin_angle = (cross_product.0 * djk.0 + cross_product.1 * djk.1 + cross_product.2 * djk.2) / (norm_cross1 * norm_cross2 * djk_norm);
        sin_angle.atan2(cos_angle)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance_zero() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (1239652.7359, -659301792.563987654, 0.00045));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (1239652.7359, -659301792.563987654, 0.00045));
        assert_eq!(Atom::distance(&atom1, &atom2), 0.0);
    }

    #[test]
    fn test_distance() {
        let atom1 = Atom::from_coords((1.0, -2.0, 3.0));
        let atom2 = Atom::from_coords((-4.0, 5.0, -6.0));
        assert_eq!(Atom::distance(&atom1, &atom2), (5.0_f64*5.0_f64 + 7.0_f64*7.0_f64 + 9.0_f64*9.0_f64).sqrt());
    }

    #[test]
    fn test_angle_zero() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (1.0, 1.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (2.0, 2.0, 2.0));
        assert_eq!(Atom::angle(&atom1, &atom2, &atom3), 0.0);
    }

    #[test]
    fn test_angle_pi() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-2.0, 2.0, -2.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (1.0, -1.0, 1.0));
        assert_eq!(Atom::angle(&atom1, &atom2, &atom3), std::f64::consts::PI);
    }

    #[test]
    fn test_angle_pi_2() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (1.0, 0.0, 0.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 1.0, 0.0));
        assert_eq!(Atom::angle(&atom1, &atom2, &atom3), std::f64::consts::PI / 2.0);
    }

    #[test]
    fn test_angle_pi_4() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (1.0, 1.0, 0.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (1.0, 1.0, 1.0));
        assert!(Atom::angle(&atom1, &atom2, &atom3) - std::f64::consts::FRAC_PI_4 < 1e-10);
    }

    #[test]
    fn test_angle_pi_4_second() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (1.0, 1.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (1.0, 1.0, 0.0));
        assert!(Atom::angle(&atom1, &atom2, &atom3) - std::f64::consts::FRAC_PI_4 < 1e-10);
    }

    #[test]
    fn test_out_of_plane_angle_zero() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-1.0, 0.0, 0.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 1.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (1.0, 0.0, 0.0));
        assert_eq!(Atom::out_of_plane_angle(&atom1, &atom2, &atom3, &atom4), 0.0);
    }

    #[test]
    fn test_out_of_plane_angle_pi_2() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 0.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 1.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (1.0, 0.0, 0.0));
        assert!((Atom::out_of_plane_angle(&atom1, &atom2, &atom3, &atom4) + std::f64::consts::FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_out_of_plane_angle_minus_pi_2() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 0.0, -1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 1.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (1.0, 0.0, 0.0));
        assert!((Atom::out_of_plane_angle(&atom1, &atom2, &atom3, &atom4) - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_out_of_plane_angle_pi_4() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 1.0, -1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (0.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (0.0, 1.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (1.0, 0.0, 0.0));
        assert!((Atom::out_of_plane_angle(&atom1, &atom2, &atom3, &atom4) - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
    }

    #[test]
    fn test_dihedral_angle_zero() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (2.0, 2.0, 2.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (1.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-1.0, 0.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (-2.0, 2.0, 2.0));
        assert_eq!(Atom::dihedral_angle(&atom1, &atom2, &atom3, &atom4), 0.0);
    }

    #[test]
    fn test_dihedral_angle_pi() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-1.0, 1.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (-1.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (2.0, 0.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (2.0, -1.0, -1.0));
        assert_eq!(Atom::dihedral_angle(&atom1, &atom2, &atom3, &atom4), std::f64::consts::PI);
    }

    #[test]
    fn test_dihedral_angle_pi_2() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-1.0, 1.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (-1.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (2.0, 0.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (2.0, -1.0, 1.0));
        assert!((Atom::dihedral_angle(&atom1, &atom2, &atom3, &atom4) - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_dihedral_angle_minus_pi_2() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-1.0, 1.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (-1.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (2.0, 0.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (2.0, 1.0, -1.0));
        assert!((Atom::dihedral_angle(&atom1, &atom2, &atom3, &atom4) + std::f64::consts::FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_dihedral_angle_pi_4() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-1.0, 1.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (-1.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (2.0, 0.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (2.0, 0.0, 1.0));
        assert!((Atom::dihedral_angle(&atom1, &atom2, &atom3, &atom4) - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
    }

    #[test]
    fn test_dihedral_angle_pi_4_second() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (2.0, 0.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (2.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-1.0, 0.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (-1.0, 1.0, 1.0));
        assert!((Atom::dihedral_angle(&atom1, &atom2, &atom3, &atom4) - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
    }

    #[test]
    fn test_dihedral_angle_3pi_4() {
        let atom1 = Atom::new(Some("H".into()), Some(1), Some(1.008), (-1.0, 1.0, 1.0));
        let atom2 = Atom::new(Some("O".into()), Some(8), Some(15.999), (-1.0, 0.0, 0.0));
        let atom3 = Atom::new(Some("H".into()), Some(1), Some(1.008), (2.0, 0.0, 0.0));
        let atom4 = Atom::new(Some("C".into()), Some(6), Some(12.011), (2.0, 0.0, -1.0));
        assert!((Atom::dihedral_angle(&atom1, &atom2, &atom3, &atom4) + 3.0 * std::f64::consts::FRAC_PI_4).abs() < 1e-10);
    }
}