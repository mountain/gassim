use crate::core::particle::DIM;
use crate::core::{Event, EventKind, Particle};
use crate::error::{Error, Result};
use rand::{rng, rngs::StdRng, Rng, SeedableRng};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// Small numeric tolerance for time and geometric checks.
const EPS_TIME: f64 = 1e-12;
const N_WALLS: usize = DIM * 2;

/// Kinematics for an axis-aligned planar wall with piecewise-constant velocity.
/// Position evolves as: pos(t) = base_pos + u * (t - base_time)
#[derive(Debug, Clone, Copy)]
struct WallKinematics {
    base_pos: f64,
    base_time: f64,
    u: f64,
}

impl WallKinematics {
    #[inline]
    fn pos_at(&self, t: f64) -> f64 {
        self.base_pos + self.u * (t - self.base_time)
    }
}

/// Simulation domain: axis-aligned box in D=3 with adiabatic walls (Phase 0) and
/// optional moving piston (Phase 1).
///
/// Walls are represented by six planes (for axis k in {0,1,2}):
/// - wall_id = 2*k     : min wall at coordinate along axis k
/// - wall_id = 2*k + 1 : max wall at coordinate along axis k
#[derive(Debug)]
pub struct Simulation {
    time_now: f64,
    box_size: [f64; DIM],
    pub particles: Vec<Particle>,
    pq: BinaryHeap<Reverse<Event>>,

    // Phase 1: moving piston/walls
    wall_min: [WallKinematics; DIM],
    wall_max: [WallKinematics; DIM],

    // Measurements (Phase 1)
    work_total: f64,                           // cumulative mechanical work on gas
    pressure_events: Vec<(f64, f64)>,          // (time, |impulse|) on piston wall
    pressure_events_all: Vec<(f64, f64, u32)>, // (time, |impulse|, wall_id) for all walls
    piston_wall_id: Option<u32>, // which wall is considered the piston (for pressure history)

    // Phase 1.5: simple histories for windowed diagnostics
    energy_history: Vec<(f64, f64)>, // (time, kinetic_energy)
    work_history: Vec<(f64, f64)>,   // (time, work_total snapshot)

    // Phase 2: thermal walls and heat accounting
    // thermal_walls[wall_id] = Some((T_bath, accommodation)) or None if adiabatic/specular
    thermal_walls: [Option<(f64, f64)>; N_WALLS],
    // Dedicated RNG stream per wall for reproducibility of thermal sampling
    rng_walls: Vec<StdRng>,
    // Heat and per-wall accounting
    heat_total: f64,
    heat_events: Vec<(f64, f64, u32)>, // (time, dQ, wall_id)
    // Per-wall accumulators
    work_by_wall: [f64; N_WALLS],
    heat_by_wall: [f64; N_WALLS],
    impulse_by_wall: [f64; N_WALLS],
}

impl Simulation {
    /// Create a new simulation with `num_particles` hard spheres of identical `radius` and `mass`
    /// inside an axis-aligned box with edge lengths `box_size` (each >= 2*radius).
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

        // Initialize static walls at t=0 with u=0
        let mut wall_min = [WallKinematics {
            base_pos: 0.0,
            base_time: 0.0,
            u: 0.0,
        }; DIM];
        let mut wall_max = [WallKinematics {
            base_pos: 0.0,
            base_time: 0.0,
            u: 0.0,
        }; DIM];
        for k in 0..DIM {
            wall_min[k].base_pos = 0.0;
            wall_max[k].base_pos = box_size[k];
        }

        // Phase 2: set up per-wall RNG streams for thermal sampling
        let mut rng_walls: Vec<StdRng> = Vec::with_capacity(N_WALLS);
        for _w in 0..N_WALLS {
            let s: u64 = rng.random();
            rng_walls.push(SeedableRng::seed_from_u64(s));
        }

        let mut sim = Self {
            time_now: 0.0,
            box_size,
            particles,
            pq: BinaryHeap::new(),
            wall_min,
            wall_max,
            work_total: 0.0,
            pressure_events: Vec::new(),
            pressure_events_all: Vec::new(),
            piston_wall_id: None,
            energy_history: Vec::new(),
            work_history: Vec::new(),

            // Phase 2 init
            thermal_walls: [None; N_WALLS],
            rng_walls,
            heat_total: 0.0,
            heat_events: Vec::new(),
            work_by_wall: [0.0; N_WALLS],
            heat_by_wall: [0.0; N_WALLS],
            impulse_by_wall: [0.0; N_WALLS],
        };
        // Initial snapshots at t=0
        let u0 = sim.kinetic_energy();
        sim.energy_history.push((sim.time_now, u0));
        sim.work_history.push((sim.time_now, sim.work_total));

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

    /// Phase 1: Set a wall's (piston) normal velocity. Velocity is piecewise-constant from the current time.
    ///
    /// - `wall_id` in [0, 5]: 2*k for min wall, 2*k+1 for max wall along axis k.
    /// - `velocity` is the wall coordinate rate along the axis.
    ///
    /// This call re-bases the wall position at the current time and updates the event queue.
    pub fn set_piston_velocity(&mut self, wall_id: u32, velocity: f64) -> Result<()> {
        let (axis, is_max) = wall_axis_side(wall_id);
        if axis >= DIM {
            return Err(Error::InvalidParam("invalid wall_id".into()));
        }
        // Re-base at current time
        if is_max {
            let cur_pos = self.wall_max[axis].pos_at(self.time_now);
            self.wall_max[axis].base_pos = cur_pos;
            self.wall_max[axis].base_time = self.time_now;
            self.wall_max[axis].u = velocity;
        } else {
            let cur_pos = self.wall_min[axis].pos_at(self.time_now);
            self.wall_min[axis].base_pos = cur_pos;
            self.wall_min[axis].base_time = self.time_now;
            self.wall_min[axis].u = velocity;
        }
        // Mark this wall as the piston for pressure history
        self.piston_wall_id = Some(wall_id);
        // Rebuild events as wall kinematics changed
        self.rebuild_event_queue()
    }

    /// Phase 1: Total mechanical work done on gas by the moving piston/walls so far.
    pub fn work_total(&self) -> f64 {
        self.work_total
    }

    /// Phase 1: Return a copy of pressure events recorded on the piston wall as (time, |impulse|).
    pub fn pressure_events(&self) -> Vec<(f64, f64)> {
        self.pressure_events.clone()
    }

    /// Pressure impulse events for all walls as (time, |impulse|, wall_id).
    pub fn pressure_events_all(&self) -> Vec<(f64, f64, u32)> {
        self.pressure_events_all.clone()
    }

    /// Phase 2: Configure a thermal wall with bath temperature T and accommodation α ∈ \[0,1].
    /// - wall_id in [0, 2*DIM-1]
    /// - T must be > 0 and finite
    /// - accommodation α in \[0,1]
    pub fn set_thermal_wall(
        &mut self,
        wall_id: u32,
        t_bath: f64,
        accommodation: f64,
    ) -> Result<()> {
        let idx = wall_id as usize;
        if idx >= N_WALLS {
            return Err(crate::error::Error::InvalidParam("invalid wall_id".into()));
        }
        if !t_bath.is_finite() || t_bath <= 0.0 {
            return Err(crate::error::Error::InvalidParam(
                "thermal wall T must be finite and > 0".into(),
            ));
        }
        if !accommodation.is_finite() || !(0.0..=1.0).contains(&accommodation) {
            return Err(crate::error::Error::InvalidParam(
                "accommodation must be in [0,1]".into(),
            ));
        }
        self.thermal_walls[idx] = Some((t_bath, accommodation));
        Ok(())
    }

    /// Phase 2: Return total accumulated heat Q (sum of ΔQ over thermal interactions).
    pub fn heat_total(&self) -> f64 {
        self.heat_total
    }

    /// Phase 2: Return a copy of heat events as (time, dQ, wall_id).
    pub fn heat_events(&self) -> Vec<(f64, f64, u32)> {
        self.heat_events.clone()
    }

    /// Per-wall total mechanical work since start.
    pub fn work_by_wall(&self) -> Vec<f64> {
        self.work_by_wall.to_vec()
    }

    /// Per-wall total pure heat since start.
    pub fn heat_by_wall(&self) -> Vec<f64> {
        self.heat_by_wall.to_vec()
    }

    /// Per-wall accumulated absolute normal impulse |Δp_n| (for pressure proxy).
    pub fn impulse_by_wall(&self) -> Vec<f64> {
        self.impulse_by_wall.to_vec()
    }

    /// Phase 1.5: Windowed mechanical pressure estimate on a wall using impulse history.
    ///
    /// Computes P ≈ (Σ |Δp|) / (A * Δt) over events within [time_now - window, time_now],
    /// where A is the wall area (product of the two orthogonal box lengths).
    ///
    /// Limitations:
    /// - Currently only supports the active piston wall (if set). Requesting a different wall_id
    ///   returns InvalidParam until per-wall histories are implemented.
    pub fn mechanical_pressure(&self, window: f64, wall_id: u32) -> Result<f64> {
        if !window.is_finite() || window <= 0.0 {
            return Err(crate::error::Error::InvalidParam(
                "window must be a positive finite float".into(),
            ));
        }
        let Some(piston_id) = self.piston_wall_id else {
            return Err(crate::error::Error::InvalidParam(
                "no piston wall configured; call set_piston_velocity first".into(),
            ));
        };
        if wall_id != piston_id {
            return Err(crate::error::Error::InvalidParam(
                "mechanical_pressure currently supports only the active piston wall_id".into(),
            ));
        }
        let t0 = self.time_now - window;
        let sum_impulse: f64 = self
            .pressure_events
            .iter()
            .filter(|(t, _)| *t >= t0)
            .map(|(_, imp)| *imp)
            .sum();
        if sum_impulse == 0.0 {
            return Err(crate::error::Error::DiagnosticInsufficientSamples(
                "no piston impulse events in window",
            ));
        }
        let area = self.area_for_wall(wall_id)?;
        Ok(sum_impulse / (area * window))
    }

    /// Phase 1.5/2: Return total mechanical work and total heat.
    pub fn work_heat(&self) -> (f64, f64) {
        (self.work_total, self.heat_total)
    }

    /// Phase 1.5: Compute drift-removed temperature tensor diagonals and scalar T from current velocities.
    ///
    /// Returns (Txx, Tyy, Tzz, T_scalar, count). Requires at least 2 selected samples.
    /// If `mask` is provided, it must have length N and select particles included in the statistics.
    pub fn temperature_tensor(&self, mask: Option<&[bool]>) -> Result<(f64, f64, f64, f64, usize)> {
        const D: usize = DIM;
        let n = self.particles.len();
        if n == 0 {
            return Err(crate::error::Error::InvalidParam("no particles".into()));
        }
        if let Some(m) = mask {
            if m.len() != n {
                return Err(crate::error::Error::InvalidParam(
                    "mask length must equal number of particles".into(),
                ));
            }
        }
        // Selection
        let mut count: usize = 0;
        let mut u = [0.0_f64; D];
        for (idx, p) in self.particles.iter().enumerate() {
            let selected = match mask {
                Some(m) => m[idx],
                None => true,
            };
            if selected {
                count += 1;
                for (k, uk) in u.iter_mut().enumerate() {
                    *uk += p.v[k];
                }
            }
        }
        if count < 2 {
            return Err(crate::error::Error::DiagnosticInsufficientSamples(
                "need at least 2 samples for temperature tensor",
            ));
        }
        for uk in u.iter_mut() {
            *uk /= count as f64;
        }

        // Accumulate centered second moments
        let m_ref = self.particles[0].mass;
        let mut sum_cc_diag = [0.0_f64; D];
        let mut sum_cc_trace = 0.0_f64;
        for (idx, p) in self.particles.iter().enumerate() {
            let selected = match mask {
                Some(m) => m[idx],
                None => true,
            };
            if !selected {
                continue;
            }
            let mut c = [0.0_f64; D];
            for k in 0..D {
                c[k] = p.v[k] - u[k];
            }
            for k in 0..D {
                let ck2 = c[k] * c[k];
                sum_cc_diag[k] += ck2;
                sum_cc_trace += ck2;
            }
        }
        // Temperature tensor components and scalar T (k_B = 1 units)
        let t_xx = m_ref * sum_cc_diag[0] / (count as f64);
        let t_yy = m_ref * sum_cc_diag[1] / (count as f64);
        let t_zz = m_ref * sum_cc_diag[2] / (count as f64);
        let t_scalar = (m_ref * sum_cc_trace) / ((DIM as f64) * (count as f64));
        Ok((t_xx, t_yy, t_zz, t_scalar, count))
    }

    /// Phase 1.5: Windowed −α/T rate using First Law closure (Phase 1 adiabatic approximation).
    ///
    /// Definition here uses the residual rate r = −(ΔU − ΔW) / (T̄ Δt)
    /// over the window [time_now − window, time_now]. Heat Q=0 in Phase 1.
    pub fn alpha_over_t(&self, window: f64) -> Result<f64> {
        if !window.is_finite() || window <= 0.0 {
            return Err(crate::error::Error::InvalidParam(
                "window must be a positive finite float".into(),
            ));
        }
        let t0 = self.time_now - window;
        let u0 = self.value_at_or_before(&self.energy_history, t0).ok_or(
            crate::error::Error::DiagnosticInsufficientSamples("insufficient U history for window"),
        )?;
        let w0 = self.value_at_or_before(&self.work_history, t0).ok_or(
            crate::error::Error::DiagnosticInsufficientSamples("insufficient W history for window"),
        )?;
        let u_now = self.kinetic_energy();
        let w_now = self.work_total;
        let delta_u = u_now - u0;
        let delta_w = w_now - w0;

        // Use current scalar temperature as proxy T̄ (k_B = 1 units)
        let (_tx, _ty, _tz, t_scalar, cnt) = self.temperature_tensor(None)?;
        if t_scalar <= 0.0 || cnt < 2 {
            return Err(crate::error::Error::DiagnosticInsufficientSamples(
                "insufficient temperature statistics",
            ));
        }
        Ok(-(delta_u - delta_w) / (t_scalar * window))
    }

    /// Search helper: last history value at or before t0.
    fn value_at_or_before(&self, hist: &[(f64, f64)], t0: f64) -> Option<f64> {
        // Linear scan backward (hist is append-only; small for tests)
        for (t, v) in hist.iter().rev() {
            if *t <= t0 {
                return Some(*v);
            }
        }
        None
    }

    /// Phase 1.5: Build histograms of velocity components for selected axes.
    ///
    /// Returns for each axis a pair (edges, counts) where:
    /// - edges has length bins+1
    /// - counts has length bins
    pub fn velocity_histogram(
        &self,
        axes: &[usize],
        bins: usize,
        range: Option<(f64, f64)>,
        mask: Option<&[bool]>,
    ) -> Result<Vec<(Vec<f64>, Vec<f64>)>> {
        if axes.is_empty() {
            return Err(crate::error::Error::InvalidParam(
                "axes must be non-empty".into(),
            ));
        }
        for &ax in axes {
            if ax >= DIM {
                return Err(crate::error::Error::InvalidParam(
                    "axis index out of bounds".into(),
                ));
            }
        }
        if bins < 1 {
            return Err(crate::error::Error::InvalidParam(
                "bins must be >= 1".into(),
            ));
        }
        let n = self.particles.len();
        if let Some(m) = mask {
            if m.len() != n {
                return Err(crate::error::Error::InvalidParam(
                    "mask length must equal number of particles".into(),
                ));
            }
        }
        // Collect per-axis min/max if range not provided
        let mut results = Vec::with_capacity(axes.len());
        for &ax in axes {
            let mut vmin = f64::INFINITY;
            let mut vmax = f64::NEG_INFINITY;
            let mut count_sel = 0usize;
            for (idx, p) in self.particles.iter().enumerate() {
                let selected = match mask {
                    Some(m) => m[idx],
                    None => true,
                };
                if !selected {
                    continue;
                }
                count_sel += 1;
                let v = p.v[ax];
                if v < vmin {
                    vmin = v;
                }
                if v > vmax {
                    vmax = v;
                }
            }
            if count_sel == 0 {
                return Err(crate::error::Error::DiagnosticInsufficientSamples(
                    "no samples selected for histogram",
                ));
            }
            let (lo, hi) = match range {
                Some((a, b)) => {
                    if !a.is_finite() || !b.is_finite() || a >= b {
                        return Err(crate::error::Error::InvalidParam(
                            "invalid histogram range".into(),
                        ));
                    }
                    (a, b)
                }
                None => {
                    if vmax <= vmin {
                        // Degenerate; widen a bit around vmin
                        (vmin - 1.0, vmax + 1.0)
                    } else {
                        (vmin, vmax)
                    }
                }
            };
            let width = (hi - lo) / (bins as f64);
            let mut edges = Vec::with_capacity(bins + 1);
            for j in 0..=bins {
                edges.push(lo + (j as f64) * width);
            }
            let mut counts = vec![0.0f64; bins];
            for (idx, p) in self.particles.iter().enumerate() {
                let selected = match mask {
                    Some(m) => m[idx],
                    None => true,
                };
                if !selected {
                    continue;
                }
                let v = p.v[ax];
                let mut bin = ((v - lo) / width).floor() as isize;
                if bin < 0 {
                    bin = 0;
                }
                if bin as usize >= bins {
                    bin = (bins as isize) - 1;
                }
                counts[bin as usize] += 1.0;
            }
            results.push((edges, counts));
        }
        Ok(results)
    }

    /// Phase 1.5: KL divergence between drift-removed component velocities and Gaussian model.
    ///
    /// For each requested axis, compute D_KL(P||Q) where:
    /// - P is the empirical histogram of centered component c_k = v_k - <v_k>
    /// - Q is Normal(0, σ²) with σ² = T/m, T the scalar temperature from temperature_tensor.
    ///   Returns the mean D_KL over axes.
    pub fn kl_stats(
        &self,
        bins: usize,
        range: Option<(f64, f64)>,
        axes: &[usize],
        mask: Option<&[bool]>,
    ) -> Result<f64> {
        if axes.is_empty() {
            return Err(crate::error::Error::InvalidParam(
                "axes must be non-empty".into(),
            ));
        }
        for &ax in axes {
            if ax >= DIM {
                return Err(crate::error::Error::InvalidParam(
                    "axis index out of bounds".into(),
                ));
            }
        }
        if bins < 5 {
            return Err(crate::error::Error::InvalidParam(
                "bins must be >= 5".into(),
            ));
        }
        // Selection and mean velocity
        let n = self.particles.len();
        if let Some(m) = mask {
            if m.len() != n {
                return Err(crate::error::Error::InvalidParam(
                    "mask length must equal number of particles".into(),
                ));
            }
        }
        let mut count = 0usize;
        let mut u = [0.0f64; DIM];
        for (idx, p) in self.particles.iter().enumerate() {
            let selected = match mask {
                Some(m) => m[idx],
                None => true,
            };
            if !selected {
                continue;
            }
            count += 1;
            for (k, uk) in u.iter_mut().enumerate() {
                *uk += p.v[k];
            }
        }
        if count < 5 {
            return Err(crate::error::Error::DiagnosticInsufficientSamples(
                "need at least 5 samples for KL",
            ));
        }
        for uk in u.iter_mut() {
            *uk /= count as f64;
        }
        // Scalar temperature and sigma^2 = T/m
        let (_tx, _ty, _tz, t_scalar, _c) = self.temperature_tensor(mask)?;
        let m_ref = self.particles[0].mass;
        let sigma2 = t_scalar / m_ref;
        if sigma2 <= 0.0 {
            return Err(crate::error::Error::DiagnosticInsufficientSamples(
                "non-positive variance for KL",
            ));
        }
        let sigma = sigma2.sqrt();

        let mut dkl_sum = 0.0f64;
        let two_pi = std::f64::consts::PI * 2.0;
        for &ax in axes {
            // Choose histogram range
            let (lo, hi) = match range {
                Some((a, b)) => (a, b),
                None => {
                    // symmetric around 0 using ±5σ of centered c
                    let r = 5.0 * sigma;
                    (-r, r)
                }
            };
            let width = (hi - lo) / (bins as f64);
            let mut counts = vec![0.0f64; bins];
            let mut total = 0.0f64;
            for (idx, p) in self.particles.iter().enumerate() {
                let selected = match mask {
                    Some(m) => m[idx],
                    None => true,
                };
                if !selected {
                    continue;
                }
                let c = p.v[ax] - u[ax];
                let mut bin = ((c - lo) / width).floor() as isize;
                if bin < 0 {
                    bin = 0;
                }
                if bin as usize >= bins {
                    bin = (bins as isize) - 1;
                }
                counts[bin as usize] += 1.0;
                total += 1.0;
            }
            if total <= 0.0 {
                return Err(crate::error::Error::DiagnosticInsufficientSamples(
                    "no samples in KL window",
                ));
            }
            // D_KL = Σ p_j ln(p_j / q_j), approximate q_j by φ(center)*Δx
            let mut dkl = 0.0f64;
            for (j, cnt) in counts.iter().enumerate().take(bins) {
                let p_j = *cnt / total;
                if p_j == 0.0 {
                    continue;
                }
                let center = lo + (j as f64 + 0.5) * width;
                let q_density =
                    (1.0 / (sigma * (two_pi).sqrt())) * (-center * center / (2.0 * sigma2)).exp();
                let mut q_j = q_density * width;
                if q_j < 1e-300 {
                    q_j = 1e-300;
                }
                dkl += p_j * (p_j / q_j).ln();
            }
            dkl_sum += dkl;
        }

        Ok(dkl_sum / (axes.len() as f64))
    }

    /// Phase 1.5: append current (time, U) and (time, W) to histories.
    #[inline]
    fn record_snapshot(&mut self) {
        let u = self.kinetic_energy();
        self.energy_history.push((self.time_now, u));
        self.work_history.push((self.time_now, self.work_total));
        // Optional: cap history growth in future with down-sampling/compaction
    }

    /// Compute wall area for the given wall id (product of the two orthogonal box lengths).
    fn area_for_wall(&self, wall_id: u32) -> Result<f64> {
        let (axis, _is_max) = wall_axis_side(wall_id);
        if axis >= DIM {
            return Err(crate::error::Error::InvalidParam("invalid wall_id".into()));
        }
        let mut area = 1.0_f64;
        for k in 0..DIM {
            if k != axis {
                area *= self.box_size[k];
            }
        }
        Ok(area)
    }

    /// Advance the simulation to `target_time` (must be ≥ current time).
    ///
    /// This implementation drifts all particles to each processed event time (full drift).
    /// It uses collision-count based invalidation to discard stale queue entries.
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
                // Snapshot at target_time
                self.record_snapshot();
                break;
            };

            let t_ev = ev.time_f64();

            if t_ev > target_time {
                // Next event occurs after the target; drift to target and stop
                self.drift_all(target_time)?;
                self.time_now = target_time;
                // Snapshot at target_time
                self.record_snapshot();
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
                    // Snapshot after event resolution
                    self.record_snapshot();

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
                    // Snapshot after event resolution/work accumulation
                    self.record_snapshot();

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

    /// Rebuild the event queue from the current particle states and wall kinematics.
    ///
    /// This should be called after externally modifying positions/velocities
    /// or wall kinematics to ensure event times are consistent.
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
    ///
    /// Moving-wall kinematics are handled via relative velocities:
    /// - Min wall: (x - pos_min_now) + (v_x - u_min) t = radius
    /// - Max wall: (pos_max_now - x) + (u_max - v_x) t = radius
    fn predict_p2w_time(&self, i: usize) -> Option<(f64, u32)> {
        let p = &self.particles[i];
        let mut best_t = f64::INFINITY;
        let mut best_wall: Option<u32> = None;

        for k in 0..DIM {
            // Min wall on axis k
            let pos_min = self.wall_min[k].pos_at(self.time_now);
            let u_min = self.wall_min[k].u;
            let s0_min = p.r[k] - pos_min;
            let rel_min = p.v[k] - u_min;
            if rel_min < -EPS_TIME {
                let t_rel = (p.radius - s0_min) / rel_min; // rel_min negative
                if t_rel > EPS_TIME {
                    let t_abs = self.time_now + t_rel;
                    if t_abs < best_t {
                        best_t = t_abs;
                        best_wall = Some((2 * k) as u32);
                    }
                }
            }

            // Max wall on axis k
            let pos_max = self.wall_max[k].pos_at(self.time_now);
            let u_max = self.wall_max[k].u;
            let s0_max = pos_max - p.r[k];
            let rel_max = u_max - p.v[k];
            if rel_max < -EPS_TIME {
                let t_rel = (p.radius - s0_max) / rel_max; // rel_max negative
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

    /// Drift all particles to the specified absolute time by linear motion,
    /// then safety-clamp to the dynamic domain at that time.
    fn drift_all(&mut self, to_time: f64) -> Result<()> {
        if to_time < self.time_now - EPS_TIME {
            return Err(Error::InvalidParam("cannot drift backwards in time".into()));
        }
        let dt = to_time - self.time_now;
        if dt.abs() <= EPS_TIME {
            return Ok(());
        }

        // Pre-compute dynamic wall positions at to_time
        let mut lo_bound = [0.0_f64; DIM];
        let mut hi_bound = [0.0_f64; DIM];
        for k in 0..DIM {
            lo_bound[k] = self.wall_min[k].pos_at(to_time);
            hi_bound[k] = self.wall_max[k].pos_at(to_time);
            // Keep bounds sane within original box in case of extreme velocities
            if lo_bound[k] < 0.0 {
                lo_bound[k] = 0.0;
            }
            if hi_bound[k] > self.box_size[k] {
                hi_bound[k] = self.box_size[k];
            }
        }

        for p in &mut self.particles {
            for k in 0..DIM {
                // Linear drift
                p.r[k] += p.v[k] * dt;

                // Safety clamp to current dynamic domain
                let lo = lo_bound[k] + p.radius;
                let hi = hi_bound[k] - p.radius;
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

    /// Resolve a particle-wall collision with a possibly moving and possibly thermal wall.
    ///
    /// - If wall is adiabatic/specular: reflect the normal component in the moving-wall frame.
    /// - If wall is thermal: with probability α, resample outgoing velocity from a half-space Maxwellian
    ///   at temperature T in the wall frame; otherwise specular reflect. Accumulate heat ΔQ.
    fn resolve_p2w(&mut self, i: usize, wall_id: u32) -> Result<()> {
        let (axis, is_max) = wall_axis_side(wall_id);

        // Wall velocity and position at the current time (event time)
        let (u_w, pos_now) = if is_max {
            (
                self.wall_max[axis].u,
                self.wall_max[axis].pos_at(self.time_now),
            )
        } else {
            (
                self.wall_min[axis].u,
                self.wall_min[axis].pos_at(self.time_now),
            )
        };

        // Copy pre-collision velocity for ΔQ accounting
        let v_before_vec = self.particles[i].v;
        let v_before_n = v_before_vec[axis];
        // Pre-collision speed squared for energy accounting
        let vsq_before: f64 = v_before_vec.iter().map(|&c| c * c).sum();

        let mut used_specular = false;

        if let Some((t_bath, accommodation)) = self.thermal_walls[wall_id as usize] {
            // Decide thermal resampling vs specular reflection
            let u: f64 = self.rng_walls[wall_id as usize].random_range(0.0..1.0);
            if u < accommodation {
                // Thermal resampling in wall frame
                let m = self.particles[i].mass;
                let sigma2 = t_bath / m;
                if !sigma2.is_finite() || sigma2 <= 0.0 {
                    return Err(Error::InvalidParam(
                        "thermal wall requires positive finite T and mass".into(),
                    ));
                }
                let sigma = sigma2.sqrt();

                // Sample tangential components: Normal(0, σ²) via Box–Muller (two independent)
                let (z0, z1) = sample_standard_normal_pair(&mut self.rng_walls[wall_id as usize]);
                // Sample normal magnitude (into gas): Rayleigh(σ); sign depends on wall side
                let r_ray = sample_rayleigh(&mut self.rng_walls[wall_id as usize], sigma);
                let sign_n = if is_max { -1.0 } else { 1.0 }; // into domain
                let c_n_plus = sign_n * r_ray;

                // Build outgoing velocity in lab frame
                let mut v_after_vec = [0.0_f64; DIM];
                // Normal component with wall velocity offset
                v_after_vec[axis] = c_n_plus + u_w;
                // Tangential axes mapping (the two axes != axis)
                let mut t_idx = 0usize;
                for (k, vk) in v_after_vec.iter_mut().enumerate() {
                    if k == axis {
                        continue;
                    }
                    let z = if t_idx == 0 { z0 } else { z1 };
                    *vk = sigma * z; // wall tangential velocity is zero
                    t_idx += 1;
                }

                // Assign velocity
                self.particles[i].v = v_after_vec;
            } else {
                // Specular branch (moving mirror)
                used_specular = true;
            }
        } else {
            // Not a thermal wall
            used_specular = true;
        }

        if used_specular {
            // Moving mirror reflection along normal
            let v_after = 2.0 * u_w - v_before_n;
            self.particles[i].v[axis] = v_after;
            // Tangential components unchanged
        }

        // Snap position to the exact contact plane accounting for radius
        if is_max {
            self.particles[i].r[axis] = pos_now - self.particles[i].radius;
        } else {
            self.particles[i].r[axis] = pos_now + self.particles[i].radius;
        }

        // Measurements: impulse, mechanical work, and pure heat split
        let m = self.particles[i].mass;
        let v_after_n = self.particles[i].v[axis];
        let delta_p = m * (v_after_n - v_before_n); // signed impulse on particle along normal

        // Per-wall impulse accumulation (pressure proxy)
        self.impulse_by_wall[wall_id as usize] += delta_p.abs();

        if let Some(piston_id) = self.piston_wall_id {
            if piston_id == wall_id {
                self.pressure_events.push((self.time_now, delta_p.abs()));
            }
        }
        // Record impulse event for all walls
        self.pressure_events_all
            .push((self.time_now, delta_p.abs(), wall_id));

        // Energy/work/heat accounting
        let vsq_after: f64 = self.particles[i].v.iter().map(|&c| c * c).sum();
        let d_e = 0.5 * m * (vsq_after - vsq_before);
        let d_w = u_w * delta_p; // mechanical work (wall -> gas)
        self.work_total += d_w;
        self.work_by_wall[wall_id as usize] += d_w;

        if self.thermal_walls[wall_id as usize].is_some() {
            let d_q = d_e - d_w; // pure heat
            self.heat_total += d_q;
            self.heat_by_wall[wall_id as usize] += d_q;
            self.heat_events.push((self.time_now, d_q, wall_id));
            debug_assert!((d_e - (d_w + d_q)).abs() <= 1e-12);
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

fn sample_standard_normal_pair<R: Rng + ?Sized>(rng: &mut R) -> (f64, f64) {
    // Box–Muller transform
    let mut u1: f64 = rng.random_range(0.0..1.0);
    // Guard against log(0)
    if u1 < 1e-16 {
        u1 = 1e-16;
    }
    let u2: f64 = rng.random_range(0.0..1.0);
    let r = (-2.0 * u1.ln()).sqrt();
    let theta = 2.0 * std::f64::consts::PI * u2;
    (r * theta.cos(), r * theta.sin())
}

fn sample_rayleigh<R: Rng + ?Sized>(rng: &mut R, sigma: f64) -> f64 {
    let mut u: f64 = rng.random_range(0.0..1.0);
    if u < 1e-16 {
        u = 1e-16;
    }
    sigma * (-2.0 * u.ln()).sqrt()
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

    #[test]
    fn piston_smoke_and_work_accumulates() -> Result<()> {
        let mut sim = Simulation::new(16, [10.0, 10.0, 10.0], 0.2, 1.0, Some(123))?;
        sim.advance_to(0.5)?;
        // Set x-max wall (axis 0 max => wall_id = 1) to move inward slowly
        sim.set_piston_velocity(1, -0.05)?;
        let w0 = sim.work_total();
        sim.advance_to(5.0)?;
        let w1 = sim.work_total();
        assert!(w1.is_finite());
        // Expect some collisions with piston; work magnitude should grow or equal
        assert!(w1 != w0 || !sim.pressure_events().is_empty());
        Ok(())
    }
}
