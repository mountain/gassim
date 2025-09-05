use crate::core::particle::DIM;
use crate::core::{Event, EventKind, Particle};
use crate::error::{Error, Result};
use rand::{rng, rngs::StdRng, Rng, SeedableRng};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// Small numeric tolerance for time and geometric checks.
const EPS_TIME: f64 = 1e-12;

/// Simulation domain: static axis-aligned box in D=3 with adiabatic walls.
///
/// Walls are implicitly represented by six planes:
/// - For axis k in {0,1,2}, wall_id = 2*k (min wall at x_k = 0), wall_id = 2*k+1 (max wall at x_k = L_k).
#[derive(Debug)]
pub struct Simulation {
    time_now: f64,
    box_size: [f64; DIM],
    pub particles: Vec<Particle>,
    pq: BinaryHeap<Reverse<Event>>,
}

impl Simulation {
    /// Create a new simulation with `num_particles` hard spheres of identical `radius` and `mass`
    /// inside a static axis-aligned box with edge lengths `box_size` (each >= 2*radius).
    ///
    /// Particles are placed with simple rejection sampling to avoid initial overlap.
    /// Initial velocities are sampled uniformly in [-1, 1] for each component.
    pub fn new(
        num_particles: usize,
        box_size: [f64; DIM],
        radius: f64,
        mass: f64,
        seed: Option<u64>,
    ) -> Result<Self> {
        if num_particles == 0 {
            return Err(Error::InvalidParam("num_particles must be > 0".into()));
        }
        if !box_size.iter().all(|&l| l.is_finite() && l > 0.0) {
            return Err(Error::InvalidParam(
                "box_size components must be finite and > 0".into(),
            ));
        }
        if !radius.is_finite() || radius <= 0.0 {
            return Err(Error::InvalidParam("radius must be finite and > 0".into()));
        }
        if !mass.is_finite() || mass <= 0.0 {
            return Err(Error::InvalidParam("mass must be finite and > 0".into()));
        }
        for &l in &box_size {
            if l < 2.0 * radius {
                return Err(Error::InvalidParam(
                    "box_size must be at least 2 * radius in every dimension".into(),
                ));
            }
        }

        let mut rng: StdRng = match seed {
            Some(s) => SeedableRng::seed_from_u64(s),
            None => SeedableRng::seed_from_u64(rng().random()),
        };

        // Rejection sample non-overlapping initial positions
        let mut particles: Vec<Particle> = Vec::with_capacity(num_particles);
        let max_attempts = 1_000_000usize;
        for id in 0..(num_particles as u32) {
            let mut attempts = 0usize;
            let r = loop {
                if attempts >= max_attempts {
                    return Err(Error::InvalidParam(format!(
                        "failed to place particle {} without overlap; try fewer particles or smaller radius",
                        id
                    )));
                }
                attempts += 1;
                // Sample positions within [radius, L - radius]
                let mut r = [0.0_f64; DIM];
                for (k, r_k) in r.iter_mut().enumerate() {
                    let lo = radius;
                    let hi = box_size[k] - radius;
                    *r_k = rng.random_range(lo..=hi);
                }
                if !overlaps_existing(&particles, &r, radius) {
                    break r;
                }
            };

            // Random initial velocities in [-1, 1] for each component
            let mut v = [0.0_f64; DIM];
            v.iter_mut().for_each(|x| *x = rng.random_range(-1.0..=1.0));

            particles.push(Particle::new(id, r, v, radius, mass)?);
        }

        let mut sim = Self {
            time_now: 0.0,
            box_size,
            particles,
            pq: BinaryHeap::new(),
        };

        sim.schedule_initial_events()?;
        Ok(sim)
    }

    /// Returns current simulation time.
    pub fn time(&self) -> f64 {
        self.time_now
    }

    /// Number of particles.
    pub fn num_particles(&self) -> usize {
        self.particles.len()
    }

    /// Positions as a Vec of fixed-size arrays.
    pub fn positions(&self) -> Vec<[f64; DIM]> {
        self.particles.iter().map(|p| p.r).collect()
    }

    /// Velocities as a Vec of fixed-size arrays.
    pub fn velocities(&self) -> Vec<[f64; DIM]> {
        self.particles.iter().map(|p| p.v).collect()
    }

    /// Advance the simulation to `target_time` (must be ≥ current time).
    ///
    /// This implementation drifts all particles to each processed event time (simpler and robust
    /// for Phase 0 with small-to-medium N). It uses collision-count based invalidation to
    /// discard stale queue entries.
    pub fn advance_to(&mut self, target_time: f64) -> Result<()> {
        if !target_time.is_finite() {
            return Err(Error::InvalidParam("target_time must be finite".into()));
        }
        if target_time < self.time_now - EPS_TIME {
            return Err(Error::InvalidParam(
                "target_time cannot be earlier than current time".into(),
            ));
        }

        while self.time_now + EPS_TIME < target_time {
            let Some(Reverse(ev)) = self.pq.pop() else {
                // No more events; drift to target_time and stop
                self.drift_all(target_time)?;
                self.time_now = target_time;
                break;
            };

            let t_ev = ev.time_f64();

            if t_ev > target_time {
                // Next event occurs after the target; drift to target and stop
                self.drift_all(target_time)?;
                self.time_now = target_time;
                // Put the event back for future calls
                self.pq.push(Reverse(ev));
                break;
            }

            // Validate against current collision counts
            match ev.kind {
                EventKind::P2P { i, j } => {
                    let (ii, jj) = (i as usize, j as usize);
                    let (cci, ccj) = (
                        self.particles[ii].collision_count,
                        self.particles[jj].collision_count,
                    );
                    if !ev.is_valid(cci, Some(ccj)) {
                        continue;
                    }

                    // Drift all to event time and resolve collision
                    self.drift_all(t_ev)?;
                    self.time_now = t_ev;
                    self.resolve_p2p(ii, jj)?;

                    // Increment cc for involved
                    self.particles[ii].bump_collision_count();
                    self.particles[jj].bump_collision_count();

                    // Reschedule events involving i and j
                    self.reschedule_for_particle(ii)?;
                    self.reschedule_for_particle(jj)?;
                }
                EventKind::P2W { i, wall_id } => {
                    let ii = i as usize;
                    let cci = self.particles[ii].collision_count;
                    if !ev.is_valid(cci, None) {
                        continue;
                    }

                    self.drift_all(t_ev)?;
                    self.time_now = t_ev;
                    self.resolve_p2w(ii, wall_id)?;

                    self.particles[ii].bump_collision_count();

                    self.reschedule_for_particle(ii)?;
                }
            }
        }

        Ok(())
    }

    /// Compute total kinetic energy (diagnostic).
    pub fn kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }

    /// Rebuild the event queue from the current particle states and box.
    ///
    /// This should be called after externally modifying positions/velocities
    /// (e.g., via Python setters) to ensure event times are consistent.
    pub fn rebuild_event_queue(&mut self) -> Result<()> {
        self.pq.clear();
        self.schedule_initial_events()
    }

    // ============ Internal helpers ============

    fn schedule_initial_events(&mut self) -> Result<()> {
        // P2P for all pairs (naive)
        let n = self.particles.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if let Some(t) = self.predict_p2p_time(i, j) {
                    let ev = Event::new(
                        t,
                        EventKind::P2P {
                            i: self.particles[i].id,
                            j: self.particles[j].id,
                        },
                        self.particles[i].collision_count,
                        Some(self.particles[j].collision_count),
                    )?;
                    self.pq.push(Reverse(ev));
                }
            }
        }

        // P2W earliest per particle
        for i in 0..n {
            if let Some((t, wall_id)) = self.predict_p2w_time(i) {
                let ev = Event::new(
                    t,
                    EventKind::P2W {
                        i: self.particles[i].id,
                        wall_id,
                    },
                    self.particles[i].collision_count,
                    None,
                )?;
                self.pq.push(Reverse(ev));
            }
        }

        Ok(())
    }

    fn reschedule_for_particle(&mut self, i: usize) -> Result<()> {
        let n = self.particles.len();

        // New P2W for i
        if let Some((t, wall_id)) = self.predict_p2w_time(i) {
            let ev = Event::new(
                t,
                EventKind::P2W {
                    i: self.particles[i].id,
                    wall_id,
                },
                self.particles[i].collision_count,
                None,
            )?;
            self.pq.push(Reverse(ev));
        }

        // New P2P for i with all others
        for j in 0..n {
            if j == i {
                continue;
            }
            let (a, b) = if i < j { (i, j) } else { (j, i) };
            if let Some(t) = self.predict_p2p_time(a, b) {
                let ev = Event::new(
                    t,
                    EventKind::P2P {
                        i: self.particles[a].id,
                        j: self.particles[b].id,
                    },
                    self.particles[a].collision_count,
                    Some(self.particles[b].collision_count),
                )?;
                self.pq.push(Reverse(ev));
            }
        }

        Ok(())
    }

    /// Predict absolute collision time for pair (i, j) if any, else None.
    fn predict_p2p_time(&self, i: usize, j: usize) -> Option<f64> {
        let pi = &self.particles[i];
        let pj = &self.particles[j];

        let mut d = [0.0_f64; DIM];
        let mut u = [0.0_f64; DIM];
        for (k, u_k) in u.iter_mut().enumerate() {
            d[k] = pj.r[k] - pi.r[k];
            *u_k = pj.v[k] - pi.v[k];
        }
        let a = dot(&u, &u);
        if a <= EPS_TIME {
            return None; // No relative motion
        }
        let b = 2.0 * dot(&d, &u);
        let r_sum = pi.radius + pj.radius;
        let c = dot(&d, &d) - r_sum * r_sum;

        // Discriminant
        let disc = b * b - 4.0 * a * c;
        if disc <= 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();

        // Candidate earlier root
        let t_rel = (-b - sqrt_disc) / (2.0 * a);
        // Require approaching (b < 0) and positive time > EPS
        if b >= 0.0 || t_rel <= EPS_TIME {
            return None;
        }

        let t_abs = self.time_now + t_rel;
        if !t_abs.is_finite() {
            return None;
        }
        Some(t_abs)
    }

    /// Predict the earliest absolute particle-wall collision time for particle i.
    /// Returns (t_abs, wall_id).
    fn predict_p2w_time(&self, i: usize) -> Option<(f64, u32)> {
        let p = &self.particles[i];
        let mut best_t = f64::INFINITY;
        let mut best_wall: Option<u32> = None;

        for (k, ((&x, &v), &l)) in p.r.iter().zip(&p.v).zip(&self.box_size).enumerate() {
            if v < -EPS_TIME {
                // Approaching min wall at 0: contact when x + v t = radius
                let t_rel = (p.radius - x) / v; // v is negative
                if t_rel > EPS_TIME {
                    let t_abs = self.time_now + t_rel;
                    if t_abs < best_t {
                        best_t = t_abs;
                        best_wall = Some((2 * k) as u32);
                    }
                }
            }
            if v > EPS_TIME {
                // Approaching max wall at L: contact when x + v t = L - radius
                let t_rel = (l - p.radius - x) / v; // v is positive
                if t_rel > EPS_TIME {
                    let t_abs = self.time_now + t_rel;
                    if t_abs < best_t {
                        best_t = t_abs;
                        best_wall = Some((2 * k + 1) as u32);
                    }
                }
            }
        }

        best_wall.map(|w| (best_t, w))
    }

    /// Drift all particles to the specified absolute time by linear motion.
    fn drift_all(&mut self, to_time: f64) -> Result<()> {
        if to_time < self.time_now - EPS_TIME {
            return Err(Error::InvalidParam("cannot drift backwards in time".into()));
        }
        let dt = to_time - self.time_now;
        if dt.abs() <= EPS_TIME {
            return Ok(());
        }
        for p in &mut self.particles {
            for k in 0..DIM {
                let lo = p.radius;
                let hi = self.box_size[k] - p.radius;

                // Linear drift
                p.r[k] += p.v[k] * dt;

                // Safety clamp: if position is outside bounds, clamp to bounds.
                // This is a fallback for numerical drift or stale events.
                // The proper wall collision handling should be via P2W events.
                if p.r[k] < lo {
                    p.r[k] = lo;
                } else if p.r[k] > hi {
                    p.r[k] = hi;
                }
            }
        }
        Ok(())
    }

    /// Resolve an elastic hard-sphere collision between particles i and j.
    fn resolve_p2p(&mut self, i: usize, j: usize) -> Result<()> {
        // Unit normal at contact from i -> j
        let (ri, rj) = (self.particles[i].r, self.particles[j].r);
        let mut n = [0.0_f64; DIM];
        for ((nk, &rjk), &rik) in n.iter_mut().zip(rj.iter()).zip(ri.iter()) {
            *nk = rjk - rik;
        }
        let dist = dot(&n, &n).sqrt();
        if dist <= EPS_TIME {
            // Nearly coincident centers - numerical pathology
            return Err(Error::MathError(
                "degenerate contact normal in P2P collision".into(),
            ));
        }
        for nk in &mut n {
            *nk /= dist;
        }

        let (mi, mj) = (self.particles[i].mass, self.particles[j].mass);
        let mut u = [0.0_f64; DIM];
        for ((uk, &vj), &vi) in u
            .iter_mut()
            .zip(self.particles[j].v.iter())
            .zip(self.particles[i].v.iter())
        {
            *uk = vj - vi;
        }
        let u_n = dot(&u, &n);

        // Update velocities (only normal component changes)
        let fi = (2.0 * mj / (mi + mj)) * u_n;
        let fj = (2.0 * mi / (mi + mj)) * u_n;
        for (k, &nk) in n.iter().enumerate() {
            self.particles[i].v[k] += fi * nk;
            self.particles[j].v[k] -= fj * nk;
        }
        Ok(())
    }

    /// Resolve a particle-wall collision by specular reflection on the hit axis.
    fn resolve_p2w(&mut self, i: usize, wall_id: u32) -> Result<()> {
        let (axis, is_max) = wall_axis_side(wall_id);
        // Mirror reflection: flip the normal component
        self.particles[i].v[axis] = -self.particles[i].v[axis];

        // Snap position to the exact contact plane
        let lo = self.particles[i].radius;
        let hi = self.box_size[axis] - self.particles[i].radius;
        if is_max {
            self.particles[i].r[axis] = hi;
        } else {
            self.particles[i].r[axis] = lo;
        }
        Ok(())
    }
}

// ============ Utility helpers ============

#[inline]
fn dot(a: &[f64; DIM], b: &[f64; DIM]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

#[inline]
fn wall_axis_side(wall_id: u32) -> (usize, bool) {
    let axis = (wall_id / 2) as usize;
    let is_max = wall_id % 2 == 1;
    (axis, is_max)
}

fn overlaps_existing(existing: &[Particle], r: &[f64; DIM], radius: f64) -> bool {
    if existing.is_empty() {
        return false;
    }
    let min_sq = (2.0 * radius) * (2.0 * radius);
    for p in existing {
        let mut d = [0.0_f64; DIM];
        for ((dk, &rk), &pk) in d.iter_mut().zip(r.iter()).zip(p.r.iter()) {
            *dk = rk - pk;
        }
        let dsq = dot(&d, &d);
        if dsq < min_sq {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn make_small_sim_ok() -> Result<()> {
        let mut sim = Simulation::new(4, [10.0, 10.0, 10.0], 0.2, 1.0, Some(1234))?;
        assert_eq!(sim.num_particles(), 4);
        assert!(sim.kinetic_energy().is_finite());
        // Advance a little to trigger some events or drift
        sim.advance_to(0.1)?;
        assert!(sim.time() >= 0.1 - 1e-12);
        Ok(())
    }

    #[test]
    fn p2w_prediction_basic() -> Result<()> {
        // One particle heading toward x-min wall
        let mut sim = Simulation::new(1, [5.0, 5.0, 5.0], 0.5, 1.0, Some(2))?;
        sim.particles[0].r = [1.0, 2.5, 2.5];
        sim.particles[0].v = [-1.0, 0.0, 0.0];
        let (t, wall) = sim.predict_p2w_time(0).expect("should hit wall");
        let (axis, is_max) = wall_axis_side(wall);
        assert_eq!(axis, 0);
        assert!(!is_max);
        // Contact when x reaches radius = 0.5, starting at 1.0 moving -1.0: dt = 0.5
        assert!((t - 0.5).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn p2p_prediction_basic() -> Result<()> {
        // Two particles along x-axis moving toward each other
        let mut sim = Simulation::new(2, [10.0, 10.0, 10.0], 0.2, 1.0, Some(3))?;
        sim.particles[0].r = [3.0, 5.0, 5.0];
        sim.particles[1].r = [7.0, 5.0, 5.0];
        sim.particles[0].v = [1.0, 0.0, 0.0];
        sim.particles[1].v = [-1.0, 0.0, 0.0];
        // Radii sum = 0.4, distance = 4.0, need to reduce to 0.4 => gap to close = 3.6, relative speed = 2, t = 1.8
        let t = sim.predict_p2p_time(0, 1).expect("should collide");
        assert!((t - 1.8).abs() < 1e-12);
        Ok(())
    }
}
