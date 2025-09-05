use crate::error::{Error, Result};

/// Fixed spatial dimension (3D).
pub const DIM: usize = 3;

/// A gas particle with hard-sphere properties in D=3.
///
/// Fields:
/// - `id`: stable identifier
/// - `r`: position vector [x, y, z]
/// - `v`: velocity vector [vx, vy, vz]
/// - `radius`: hard-sphere radius (> 0)
/// - `mass`: particle mass (> 0)
/// - `collision_count`: incremented each time the particle participates in a realized event
#[derive(Debug, Clone)]
pub struct Particle {
    /// Stable particle identifier.
    pub id: u32,
    /// Position (x, y, z).
    pub r: [f64; DIM],
    /// Velocity (vx, vy, vz).
    pub v: [f64; DIM],
    /// Hard-sphere radius (> 0).
    pub radius: f64,
    /// Mass (> 0).
    pub mass: f64,
    /// Collision participation counter (for event invalidation).
    pub collision_count: u64,
}

impl Particle {
    /// Create a new particle after validating invariants.
    ///
    /// Errors:
    /// - `Error::InvalidParam` if `radius` or `mass` is non-positive or any component is NaN/inf.
    pub fn new(id: u32, r: [f64; DIM], v: [f64; DIM], radius: f64, mass: f64) -> Result<Self> {
        if !radius.is_finite() || radius <= 0.0 {
            return Err(Error::InvalidParam("radius must be finite and > 0".into()));
        }
        if !mass.is_finite() || mass <= 0.0 {
            return Err(Error::InvalidParam("mass must be finite and > 0".into()));
        }
        if !r.iter().all(|x| x.is_finite()) {
            return Err(Error::InvalidParam("position must be finite".into()));
        }
        if !v.iter().all(|x| x.is_finite()) {
            return Err(Error::InvalidParam("velocity must be finite".into()));
        }
        Ok(Self {
            id,
            r,
            v,
            radius,
            mass,
            collision_count: 0,
        })
    }

    /// Increment the collision counter (used for event invalidation).
    #[inline]
    pub fn bump_collision_count(&mut self) {
        self.collision_count = self.collision_count.saturating_add(1);
    }

    /// Returns the particle's kinetic energy: 1/2 m |v|^2.
    #[inline]
    pub fn kinetic_energy(&self) -> f64 {
        let vsq: f64 = self.v.iter().map(|&c| c * c).sum();
        0.5 * self.mass * vsq
    }

    /// Set position (validated as finite).
    pub fn set_position(&mut self, r: [f64; DIM]) -> Result<()> {
        if !r.iter().all(|x| x.is_finite()) {
            return Err(Error::InvalidParam("position must be finite".into()));
        }
        self.r = r;
        Ok(())
    }

    /// Set velocity (validated as finite).
    pub fn set_velocity(&mut self, v: [f64; DIM]) -> Result<()> {
        if !v.iter().all(|x| x.is_finite()) {
            return Err(Error::InvalidParam("velocity must be finite".into()));
        }
        self.v = v;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_particle_ok() -> Result<()> {
        let p = Particle::new(1, [0.0, 1.0, 2.0], [2.0, -3.0, 0.5], 0.5, 2.0)?;
        assert_eq!(p.id, 1);
        assert_eq!(p.r, [0.0, 1.0, 2.0]);
        assert_eq!(p.v, [2.0, -3.0, 0.5]);
        assert_eq!(p.radius, 0.5);
        assert_eq!(p.mass, 2.0);
        assert_eq!(p.collision_count, 0);
        Ok(())
    }

    #[test]
    fn invalid_radius_rejected() {
        let err = Particle::new(0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, 1.0).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("radius"));
    }

    #[test]
    fn invalid_mass_rejected() {
        let err = Particle::new(0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("mass"));
    }

    #[test]
    fn kinetic_energy_computed() -> Result<()> {
        // v = (3,4,0), |v|^2 = 25; KE = 0.5 * m * 25
        let p = Particle::new(7, [0.0, 0.0, 0.0], [3.0, 4.0, 0.0], 1.0, 2.0)?;
        assert!((p.kinetic_energy() - 25.0).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn bump_collision_count() -> Result<()> {
        let mut p = Particle::new(1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 1.0)?;
        assert_eq!(p.collision_count, 0);
        p.bump_collision_count();
        assert_eq!(p.collision_count, 1);
        Ok(())
    }
}
