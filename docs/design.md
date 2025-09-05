# GasSim Design Document

This document specifies the architecture, core data structures, event-driven algorithms, performance and robustness strategies, Python bindings and packaging, testing and validation plan, and phased extensions for GasSim. It serves as the authoritative technical blueprint for building the simulator.

This document aligns with:
- Roadmap and milestones in docs/plan.md
- Engineering standards in AGENT.md
- Current scaffolding in Cargo.toml, src/lib.rs, and pyproject.toml

Note: The current Rust module exposes a minimal PyO3 function for scaffolding only. The design below describes the planned core engine and Python API for Phases 0–2.

## 1. Goals and Scope

- Primary goals
  - Implement a physically faithful and computationally robust Event-Driven Molecular Dynamics (EDMD) simulator for a hard-sphere gas (NVE in Phase 0).
  - Keep compute-intensive physics in Rust for performance and memory safety; provide a clean Python API (via PyO3) for configuration, control, and analysis.
  - Extend to a moving piston (mechanical work) and thermal walls (heat exchange) to reproduce canonical non-equilibrium thermodynamic processes and support statistical physics validations.

- Non-goals (current)
  - GPU or distributed execution.
  - Complex intermolecular potentials (focus on hard-sphere collisions).
  - Arbitrary complex geometries (start with axis-aligned boxes).
  - Built-in visualization front-end (handled from Python, e.g., Matplotlib).

- Quality and maintainability
  - Follow AGENT.md: rustfmt, clippy (warns as errors), unit and integration tests, complete Rust doc comments, and custom error handling without unwrap/expect in application logic.
  - Modular, testable design; descriptive naming; single-responsibility functions.

## 2. High-Level Architecture

- Rust core (compute engine)
  - Event prediction: particle–particle (P2P), particle–wall (P2W) for static walls, moving piston walls, and thermal walls.
  - Event queue: priority queue over event time (min-heap behavior).
  - Main event loop: deterministic time advancement and event resolution.
  - Collision resolution and state updates (elastic hard-sphere).
  - Measurements: kinetic energy, momentum impulse (pressure proxy), mechanical work, and heat.

- Python bindings (control and analysis)
  - Expose class `GasSim` to Python with initialization, advancement to target time, state access as NumPy arrays, and accessors for observables and histories.
  - Use PyO3 and maturin for packaging and distribution.

- Build and packaging
  - Cargo.toml: `crate-type = ["cdylib"]`, `pyo3 = "0.25"`.
  - pyproject.toml: maturin backend, `pyo3/extension-module` feature, Python ≥ 3.11, CPython/PyPy classifiers.

### Thermal contact and local thermometry (ADD)

- Core engine supports two thermal boundary primitives:
  1) Thermal Wall (Maxwell accommodation at a box face; piston may simultaneously be moving and thermal).
  2) ThermoPlane (an infinitesimally thin internal plane for in‑bulk contact thermometry; transmissive Maxwell “patch” with probability α).

- Observability plane/channel for heat events dQ parallels the work/impulse channel; enables contact‑temperature readout via zero‑heat‑flux condition and closes the δ‑circuit with dU = δQ − p dV. The contact 1‑form dU + p dV − T dS maps to AEG’s α = da − μ du − λ a dv under the dictionary (S,V;U,T,p) = (v, −u; a, λ a, μ), guaranteeing structural consistency.

## 3. Recommended Project Layout

- src/
  - lib.rs (PyO3 module entry; binds Rust API to Python)
  - core/
    - particle.rs
    - wall.rs
    - event.rs
    - sim.rs (Simulation state and event loop)
    - pq.rs (priority queue wrapper)
    - measure.rs (measurements/observables)
    - math.rs (geometry and numerics helpers)
  - error.rs (custom error enum via thiserror)
  - ffi/ (Rust↔Python bridges and zero-copy array exports)
- tests/ (Rust integration tests)
- docs/
  - design.md (this document)
  - plan.md (roadmap and phases)

Future crates (planned):
- thiserror (custom error types)
- numpy (zero-copy NumPy array exports)
- rand (initialization and thermal wall sampling)
- parry (broad-phase acceleration; planned)
- criterion (benchmarks; optional)
- log/env_logger (logging)

## 4. Physical Model and Assumptions

- Hard-sphere gas with elastic collisions (no potential energy).
- Particle attributes: position r, velocity v, radius a, mass m (initially uniform a and m).
- Dimensions: initial implementation targets 2D; design remains 3D-ready (D ∈ {2,3}).
- Boundaries by phase:
  - Phase 0: static axis-aligned box with adiabatic walls.
  - Phase 1: introduce a moving piston (wall with prescribed velocity).
  - Phase 2: introduce thermal walls (heat bath behavior).
- Units: allow dimensionless or SI-consistent units; configuration records time/length/mass scales.

## 5. Core Data Structures

- Particle
  - id: u32
  - r: VecD (f64[D])
  - v: VecD (f64[D])
  - radius: f64
  - mass: f64
  - collision_count: u64 (incremented after each realized event; used to invalidate stale events)

- Wall
  - id: u32
  - kind: enum WallKind { Adiabatic, Moving { u: VecD }, HeatBath { T: f64 } }
  - n: VecD (unit outward normal)
  - pos: f64 (plane position along normal; e.g., x = pos)

- Event
  - time: f64
  - kind: enum EventKind { P2P { i: u32, j: u32 }, P2W { i: u32, wall_id: u32 } }
  - cc_i: u64, cc_j: Option<u64> (collision_count snapshots for invalidation)
  - Priority: sorted by time ascending; use `std::cmp::Reverse` or a custom Ord wrapper over BinaryHeap.

- Simulation
  - time_now: f64
  - particles: Vec<Particle>
  - walls: Vec<Wall>
  - pq: BinaryHeap<Reverse<Event>> (or a dedicated min-heap type)
  - meas: Measurements (accumulators and histories)
  - params: SimParams (dimension, tolerances, RNG seed, etc.)
  - last_update_time per particle (for lazy drift; optional)

- Measurements
  - work_total: f64
  - heat_total: f64
  - impulse_accum_piston: f64 (for pressure estimate)
  - histories: optional ring buffers/time series (windowed or down-sampled)
  - derived on-demand: kinetic energy, temperature proxy, etc.

### Additional/extended types (ADD)

- Thermal boundary kinds:

```rust
enum ThermalBcKind {
    // Classical wall (boundary of domain)
    Wall { T: f64, accommodation: f64 },             // accommodation ∈ [0,1], 1=Maxwell, 0=specular
    // Internal transmissive patch for contact thermometry
    Plane { T: f64, accommodation: f64, normal: VecD, point: VecD, area: f64 },
}
```

- Boundary container:

```rust
struct Boundary {
    id: u32,
    kind: ThermalBcKind,
    // geometry: plane x·n = pos for faces; for Plane, (normal, point, finite area or infinite flag)
    n: VecD,
    pos: f64,      // for axis-aligned faces
}
```

- EventKind addition (covers both thermal Wall and Plane):

```rust
// ADD to existing EventKind:
P2Thermal { i: u32, bnd_id: u32 }
```

- Measurements extensions:

```rust
struct HeatEvent { t: f64, dQ: f64, bnd_id: u32 }
struct PressureEvent { t: f64, impulse_abs: f64, wall_id: u32 }
struct WindowCfg { dt_window: f64, max_len: usize }
struct Moments { /* D known at compile time in Rust module scope */
    sum_c: VecD,
    // temperature tensor accumulators
    // in 2D this would be 2x2, in 3D 3x3; represent as flat array or small fixed matrix
    // shown here conceptually:
    // sum_cc[j][k] = Σ c_j c_k over samples
}
```

Accumulators:
- heat_total: f64; heat_history: RingBuffer<HeatEvent>
- work_total: f64; pressure_history: RingBuffer<PressureEvent>
- optional sliding-window velocity moments for anisotropy metrics and KL diagnostics (see §10).

## 6. Event Prediction and Resolution

### 6.1 Particle–Particle (P2P) collision time

Given i, j:
- Relative displacement d = r_j − r_i
- Relative velocity u = v_j − v_i
- Contact at |d + u t| = R = a_i + a_j
- Quadratic: a t^2 + b t + c = 0
  - a = u·u
  - b = 2 d·u
  - c = d·d − R^2
- Conditions:
  - If a ≈ 0 (|u|≈0): no relative motion → no collision.
  - Discriminant Δ = b^2 − 4ac; if Δ < 0: no real root → no collision.
  - Candidate root t = (-b - sqrt(Δ)) / (2a); require t > ε and b < 0 (approaching).

Resolution (elastic):
- At contact, unit normal n = (r_j - r_i) / |r_j - r_i|
- Standard elastic update (conserves momentum and kinetic energy):
  - Let u_n = (v_j - v_i)·n
  - With masses (m_i, m_j), update the normal components accordingly (center-of-mass frame reflection).
  - A convenient form:
    - v_i' = v_i + (2 m_j / (m_i + m_j)) u_n n
    - v_j' = v_j - (2 m_i / (m_i + m_j)) u_n n
  - Only the normal components change; tangential components are preserved.

### 6.2 Particle–Wall (P2W) collision time — static wall

- Signed distance to plane s0 = n·r - pos
- Normal velocity component u_n = n·v
- Contact condition: s0 + u_n t = ± radius (choose sign for the nearer face)
- Require that u_n points toward the wall and t > ε; otherwise no collision.

Resolution (mirror reflection):
- v' = v - 2 (v·n) n

### 6.3 Particle–Wall (P2W) collision time — moving wall (piston)

- Wall velocity u_w (only normal component matters physically)
- In wall frame: u_n_rel = (v - u_w)·n
- Contact: s0 + u_n_rel t = ± radius → t = (±radius - s0) / u_n_rel, valid if u_n_rel < 0 and t > ε

Resolution:
- Reflect in wall frame, then transform back:
  - v_rel' = v_rel - 2 (v_rel·n) n, where v_rel = v - u_w
  - v' = v_rel' + u_w

### 6.4 Thermal wall (heat bath)

- On collision with a thermal wall at temperature T_b:
  - Resample outgoing velocity from a half-space Maxwellian distribution corresponding to T_b (only velocities pointing into the domain).
  - Normal component sampled from the flux-weighted distribution; tangential components from Maxwellian at T_b.
- Accumulate heat: ΔE = E_after - E_before.

#### 6.4.1 Detailed spec — Maxwell accommodation (ADD)

- Collision rule (accommodation α): with probability α, resample in the wall frame; with probability 1−α, do specular reflection.
  - Let wall normal n (into gas), wall velocity U (piston allowed). Pre‑collision: c⁻ = v⁻ − U.
  - Resampling at bath temperature T_b (m = mass, k_B = 1 in code units):
    - Tangential components: c⁺_{t} ~ Normal(0, σ²), with σ² = T_b / m (1 tangential component in 2D; 2 in 3D).
    - Normal (into gas): c⁺_n ~ Rayleigh(σ) (flux‑weighted half‑space Maxwell).
  - Transform back: v⁺ = c⁺ + U.
- Heat accounting: ΔQ = (m/2)(|v⁺|² − |v⁻|²) (positive when wall → gas). Push (t, ΔQ, bnd_id) into heat_history; accumulate heat_total.
- Notes: The Rayleigh sampler can use inverse‑transform r = σ sqrt(−2 ln(1−u)).

### 6.5 Event invalidation

- Each particle maintains `collision_count` (cc) incremented after any realized event it participates in.
- Events store cc snapshots at scheduling time.
- On pop: if a participant’s current cc ≠ stored cc, the event is invalid; discard it.

### 6.6 ThermoPlane (internal contact thermometer) (ADD)

- Geometric predicate: free‑flight line–plane intersection (like P2W), but the plane is transmissive.
- On crossing, with probability α perform transmissive Maxwellization (same sampler as thermal wall but keep the sign of the normal so the particle continues across the plane); otherwise pass through unchanged.
- Log each modified crossing as a heat event with the same ΔQ definition.
- This enables in‑bulk contact thermometry by tuning T_b so that the window‑averaged ĖQ → 0 on that plane; see Python control loop in §11.

## 7. Event Loop and Time Advancement

Initialization:
- For each particle, predict the earliest P2P contacts (Phase 0 may do all pairs) and P2W contacts, and push the soonest ones into the queue.

Advancement (`advance_to(target_time)`):
- While `time_now < target_time`:
  - Pop event e with the minimal time.
  - If `e.time > target_time`: drift all particles to `target_time` (or lazily drift involved particles), set `time_now = target_time`, break.
  - Validate cc snapshots; if invalid → discard and continue.
  - Drift the involved entities to `e.time` (lazy drift: only particles touched by e).
  - Resolve the collision (update velocities; apply measurements).
  - Increment collision_count for involved particles.
  - Reschedule events for involved particles against relevant candidates and walls; push the new earliest events.

Drift strategies:
- Lazy/on-demand: store `last_update_time` per particle; update only touched particles to the current time to reduce global work.
- Full drift (optional, for clarity at small N): update all particle positions each event.

Determinism:
- Provide tie-break rules when event times are equal within tolerance (e.g., prefer P2P over P2W, then by id ordering) to ensure reproducible trajectories.

Pseudocode (simplified):
```
loop:
  e = pq.pop_min()
  if e.time > T: drift_all_to(T); break
  if !valid(e): continue
  drift_involved_to(e.time)
  resolve(e)
  bump_collision_counts(e.involved)
  measure.on_event(e)
  reschedule_for(e.involved)
```

### Measurement hooks (ADD)

- After every resolved event:
  - If P2W on a piston wall: accumulate impulse and mechanical work.
  - If P2Thermal: accumulate ΔQ.
  - Update sliding‑window velocity moments for anisotropy/KL if enabled (O(1) amortized per event with Welford‑style updates).
- Provide a low‑overhead window compactor (down‑sample) to keep histories bounded.

## 8. Numerical Robustness and Determinism

- Tolerances: use combined absolute/relative ε (e.g., 1e-12) for time and geometric checks to avoid duplicate events at t ≈ 0 and grazing contacts.
- Normalization: guard against near-zero-length normals (e.g., if |d| is tiny).
- Stable math: factor common terms, avoid catastrophic cancellation where possible; centralize sensitive operations in math.rs.
- Priority ordering: stable tie-breaking by event kind and ids for determinism.
- Bounds checking: ensure positions remain inside the domain after drift (with tolerance) and adjust if necessary.

### Thermal randomness & reproducibility (ADD)

- Thermal sampling uses a dedicated RNG stream per boundary id; seeds recorded in metadata for replay.
- Tie‑break precedence: P2P > P2W > P2Thermal on equal time within ε, then by (i, j, wall_id/bnd_id) index for determinism.

## 9. Performance Design

- Complexity:
  - Naive all-pairs prediction is O(N^2 log N) with a priority queue; acceptable for functional validation at small/medium N in Phase 0.
  - Introduce acceleration later:
    - Cell lists / uniform grids to reduce candidate pairs.
    - Broad-phase via Parry (planned) for generalization.

- Data layout:
  - Start with AoS (Vec<Particle>). For large N, consider migrating hot paths to SoA to improve cache locality.

- Parallelism:
  - EDMD is mostly sequential due to causality in the event loop.
  - Parallelize measurements/aggregation where possible; consider batched diagnostics.

- Python GIL:
  - Release GIL during long-running `advance_to` using `Python::allow_threads` to avoid blocking Python.

### Selective diagnostics (ADD)

- Temperature‑tensor and KL diagnostics are opt‑in and can be computed on a spatial mask (bulk vs near‑wall), or on a coarse slab grid to cap cost.
- For large N, export periodic moments instead of full velocities; Python builds tensors and histograms.

## 10. Measurements and Observables

- Kinetic energy: E_k = Σ ½ m |v|^2 (no potential energy for hard spheres).
- Pressure and work (Phase 1):
  - On piston collisions, accumulate normal momentum change Δp.
  - Estimate instantaneous pressure P_inst ≈ Δp / (A Δt) where A is piston area (length in 2D).
  - Mechanical work W accumulates via force·displacement or via impulse with piston kinematics.
- Heat (Phase 2):
  - On thermal wall collisions, ΔE contributes to the cumulative heat Q.
- Histories:
  - Maintain windowed/down-sampled time series for pressure, work, heat as requested by the API.

### 10.1 Temperature definitions (operational) (ADD)

- Drift‑removed kinetic temperature (scalar):
  - u = ⟨v⟩, c_i = v_i − u, T = (m / (D N)) Σ_i ||c_i||².
- Temperature tensor:
  - T_{jk} = (m / N) Σ_i c_{i,j} c_{i,k}, and T = (1/D) tr(T).
- Anisotropy indicators:
  - η_∞ = max_j |T_{jj}/T − 1|, and η_F = ||T − T I||_F / (T √D).

### 10.2 Pressure and work (clarification) (ADD)

- Mechanical pressure at piston: use work–volume slope p_mech = −ΔW / ΔV or windowed impulse/area/time (already present).
- For EOS checks, compare bulk p_th ≈ (1/D) tr(P) ≈ n k_B T vs wall‑normal measure; differences quantify near‑wall anisotropy.

### 10.3 Heat and contact temperature (ADD)

- Heat flux per boundary:
  - ĖQ_b ≈ Σ_{k in window} ΔQ_k / Δt.
- Contact thermometry (operational): a Python controller adjusts boundary T_b until ĖQ_b → 0 → returns contact temperature T_contact for that site/plane.

### 10.4 Distribution diagnostics (KL/KS) (ADD)

- 1D (per component) KL: histogram P of c_x vs Q = Normal(0, σ²) with σ² = T/m:
  - D_KL(P‖Q) = Σ_j p_j ln(p_j / q_j).
- 3D speed version: compare radial histogram to Maxwell speed PDF.
- Report time series of KL (bulk and near‑wall).

### 10.5 δ‑Circuit measurement (non‑equilibrium entropy accounting) (ADD)

- Baseline (ideal monatomic gas, k_B = 1):
  - dS_IG = (1/T) dU + (p_th/T) dV, and ΔS_IG = N [ (3/2) ln(T₂/T₁) + ln(V₂/V₁) ].
- δ‑circuit 1‑form:
  - α_IG = dU − (−p_mech dV + T dS_IG).
  - The vertical residual −α_IG / T is the running entropy‑production measure; integrate over time/windows.
- AEG consistency: under the contactomorphism (S,V;U,T,p) = (v, −u; a, λ a, μ), α_TD = dU + p dV − T dS pulls back to α = da − μ du − λ a dv, so the δ‑circuit is geometrically faithful.

## 11. Python API (PyO3)

Module name: `gassim`. Expose a `GasSim` class.

Phase 0 (minimal):
- Initialization:
  - `GasSim(num_particles: int, box_size: Sequence[float], radius: float = 1.0, mass: float = 1.0, dim: int = 2, seed: int | None = None, ...)`
  - Optional initial velocity distribution (uniform/normal) and target temperature; may be staged in later minor releases.
- Advancement:
  - `advance_to(target_time: float) -> None` (releases GIL internally).
- State access:
  - `get_positions() -> np.ndarray[float64]  # shape (N, D)`
  - `get_velocities() -> np.ndarray[float64] # shape (N, D)`

Phase 1 additions:
- Moving wall configuration at init or via method:
  - `set_piston(wall_id: int, velocity: float | Sequence[float]) -> None` (initially piecewise-constant velocity).
- Observables:
  - `get_work_done() -> float`
  - `get_pressure_history(window: float | None = None) -> np.ndarray`

Phase 2 additions:
- Thermal walls:
  - `set_thermal_wall(wall_id: int, T: float) -> None`
- Observables:
  - `get_heat_flow() -> float`

Implementation notes:
- Zero-copy arrays via `numpy` crate to expose internal buffers where safe; else copy.
- Error mapping: Rust `Error` → `PyErr` with clear messages for invalid parameters, out-of-bounds, or state errors.

### Phase 1.5 diagnostics (ADD)

- `get_temperature_tensor(mask: Optional[np.ndarray]) -> (Txx, Tyy, Tzz, T_scalar, count)`
- `get_velocity_histogram(axes=("x","y","z"), bins=80, range=None, mask=None) -> dict`
- `get_kl_stats(mask=None) -> float` (returns current window KL for chosen mode)

### Thermal contact / thermometry (Phase 2) (ADD)

- `set_thermal_wall(wall_id: int, T: float, accommodation: float = 1.0) -> None`
- `set_thermo_plane(plane_id: int, point: Sequence[float], normal: Sequence[float], T: float, accommodation: float = 0.02, area: Optional[float] = None) -> None`
- `get_heat_history(boundary_id: int, window: float) -> np.ndarray  # [t, dQ]`
- Controller pattern (Python‑side): `find_contact_temperature(boundary_id, tol=..., max_iter=...)` adjusts T to drive ⟨ĖQ⟩ → 0.

### δ‑circuit helpers (ADD)

- `get_mechanical_pressure(window: float, wall_id: int) -> float`
- `get_entropy_IG(time_series: bool = False) -> (S_vals or ΔS)`
- `get_alpha_over_T(window: float) -> float  # windowed −α/T`
- `get_work_heat() -> (W_total, Q_total)`

Implementation notes: all data are windowed; masks allow bulk vs near‑wall selection.

## 12. Error Handling and Logging

- Errors:
  - Define a comprehensive `Error` enum using `thiserror` (e.g., InvalidParam, MathError, NoCollision, OutOfBounds, NotImplemented).
  - All fallible functions return `Result<T, Error>`; no unwrap/expect in application logic.

- Logging:
  - Integrate `log` facade with optional `env_logger` for debug builds.
  - Consider Python-side toggles to control verbosity.

- Documentation:
  - All public items have Rust doc comments (///) covering purpose, parameters, returns, and error variants, enabling complete `cargo doc`.

### Diagnostics guardrails (ADD)

- KL/anisotropy modules validate binning and sample counts; return a dedicated DiagnosticInsufficientSamples error if below threshold.

## 13. Testing and Validation

Rust unit tests (in-module #[cfg(test)]):
- event.rs: P2P/P2W collision time edge cases (no collision, grazing, equal speed, t≈0); resolution conserves momentum and kinetic energy.
- sim.rs: event invalidation, ordering determinism under tie/tolerance, and long-run stability of totals (e.g., no particle escapes due to drift error).

Rust integration tests (in tests/):
- Small simulations advancing to a given time; assert:
  - Energy conservation within tolerance over long NVE runs (Phase 0).
  - Velocity distribution approaches Maxwell–Boltzmann at equilibrium (statistical checks on histogram/mean/variance).
  - Isotropy at equilibrium: <v_x^2> ≈ <v_y^2> (and ≈ <v_z^2> in 3D).

Python-level tests (optional, via maturin develop):
- Import `gassim`, construct `GasSim`, run `advance_to`, inspect arrays and basic invariants.
- Phase 1: compression/expansion with piston; check quasi-static PV^γ constancy (γ = 5/3 for monatomic hard spheres) as piston speed → 0; for finite rates, entropy change ΔS > 0 qualitatively.
- Phase 2: validate First Law ΔE = W + Q; steady-state conduction with T_L ≠ T_R shows heat flux J_q and negative temperature gradient with preliminary Fourier behavior.

CI gates:
- `cargo fmt --check`
- `cargo clippy -- -D warnings`
- `cargo test`
- `cargo doc --no-deps`

### Additional unit/integration tests (ADD)

1) Thermal sampler tests: moments of sampled c⁺ match Maxwell half‑space within statistical error; Rayleigh law for c_n.
2) First Law closure (Phase 2): ΔU = Σ_b Q_b + W over runs (tolerance relative to U).
3) Contact temperature: single thermal wall at T_b relaxes to T_b; thermo‑plane controller recovers T_contact with |T_contact − T_kin| bulk median ≤ 5%.
4) δ‑circuit “must‑pass” (Phase 1/2):
   - Same‑end‑volume compression ensemble: all runs share V_end.
   - Quasi‑static limit: ⟨−α/T⟩ increases monotonically with |piston speed|, → 0 as speed → 0.
   - Work–energy residual small (|ΔU − W| in Phase 1 adiabatic; |ΔU − W − Q| in Phase 2).
   - Baseline correlation: deviations in T and wall pressure vs adiabatic/EOS baselines positively correlate with ⟨−α/T⟩.

### Diagnostic thresholds (defaults) (ADD)

- Bulk anisotropy: η_F < 0.1 typical; warn if exceeded persistently.
- KL (bulk): D_KL < 1e−2 typical for near‑equilibrium; warn on sustained spikes.

## 14. Reproducibility and Randomness

- Seed management:
  - Initialization accepts an explicit RNG seed; record it in simulation metadata for reproducibility.
  - Thermal wall resampling uses the same RNG stream (thread-safe as needed).

- Metadata:
  - Expose `get_metadata()` returning version, build info, parameters, dimension, and seed.

### Per‑boundary RNG seeds and controllers (ADD)

- Record per‑boundary RNG seeds and controller states (e.g., find_contact_temperature iterations) in metadata for exact reruns.

## 15. Versioning and Packaging

- Semantic versioning:
  - 0.x may include breaking changes; Python API stability is preferred where feasible.
- Packaging:
  - Built and published via maturin. Python ≥ 3.11; CPython and PyPy supported.
- Module/package name: `gassim`.

## 16. Milestones and Alignment with docs/plan.md

- Phase 0 (Foundational EDMD, NVE):
  - Implement P2P and static P2W predictions/resolutions, event loop, and measurements for energy.
  - Validate energy conservation, Maxwell–Boltzmann distribution, isotropy.
  - Expose minimal Python API (init, advance_to, get_positions/velocities).

- Phase 1 (Mechanical work — moving piston):
  - Moving wall kinematics integrated into predictions/resolutions.
  - Impulse accumulation for pressure; compute work; pressure history.
  - Validation: quasi-static PV^γ constancy; ΔS > 0 for finite-rate processes.

- Phase 2 (Heat exchange — thermal wall):
  - Thermal wall resampling and heat accumulation.
  - Validation: First Law ΔE = W + Q; steady-state conduction; preliminary Fourier’s law.

### Phase 1.5 (Diagnostics & δ‑circuit) (ADD)

- Deliver anisotropy/KL diagnostics, δ‑circuit accounting, and “must‑pass” CI gate.
- Scripts: `examples/phase1_adiabatic_alpha.py`, `examples/diagnostics_qc.py`.

### Phase 2 details (Thermal contact & thermometry) (ADD)

- Implement Thermal Wall and ThermoPlane; expose heat histories; First‑Law regression.
- Example: piston compression with one thermal wall; conduction test between T_L ≠ T_R.

## 17. Risks and Mitigations

- Numerical pathologies (grazing collisions, t≈0):
  - Unified tolerances; tie-breaking; guard against repeated immediate events; normalize with care.
- Event storms and stale queues:
  - Periodic partial or full rescheduling; maintain counters and heuristics to rebuild.
- API evolution:
  - Track changes in this design doc; provide migration notes when altering Python/Rust APIs.
- Performance at scale:
  - Plan for cell lists and broad-phase acceleration (Parry); evolve towards SoA if necessary.

## 18. Future Enhancements

- Local field analysis:
  - Spatial binning (slabs/cells) for density ρ(x), velocity u(x), temperature T(x).
- Complex phenomena:
  - Shock waves from high-speed piston; shear flows via sliding walls or Lees–Edwards BCs.
- Parallelism:
  - Parallel measurement aggregation; multi-instance ensembles for statistics.
- Visualization and tooling:
  - Convenience exporters to NumPy/Parquet/NetCDF; example notebooks for analyses.

### Additional local fields and geometry (ADD)

- Slab/cell binning to export ρ(r), u(r), T(r) and contact‑temperature field from sparse thermo‑planes (kriging/IDW blend with kinetic‑temperature map).
- Variable accommodation, moving thermo‑planes (Lagrangian probes), and specular/roughness models.

## 19. Appendix: Notation

- r, v: position and velocity vectors (ℝ^D)
- a, m: particle radius and mass
- n: wall outward unit normal
- pos: wall plane position along n (e.g., x = pos)
- R: sum of radii at contact, R = a_i + a_j
- ε: numerical tolerance (absolute/relative)
- Δp, W, Q: impulse, work, and heat accumulators

## 20. Appendix: Maxwell sampling cookbook (ADD)

- Tangential components: Box–Muller or Ziggurat for Normal(0, σ²).
- Positive normal (flux‑weighted): sample R ~ Rayleigh(σ), set c⁺_n = R.
- Specular mixing: with probability 1−α, reflect c⁺_n = −c⁻_n, keep tangentials.
- Heat increment: ΔQ = (m/2)(||v⁺||² − ||v⁻||²).
- Contact controller (pseudo‑code):

```python
# every Δt window
qdot = sum(dQ) / Δt
if abs(qdot) < tol:
    done
T_b -= kappa * qdot  # simple integral control; or bracket root by bisection
```

### Notes on Theory Alignment (for the Methods section of the paper)

- The simulator’s δ‑circuit is the numerical instantiation of the contact 1‑form α_TD = dU + p dV − T dS; under the AEG dictionary (S,V;U,T,p) = (v, −u; a, λ a, μ), it pulls back to α = da − μ du − λ a dv. The horizontal part δa = μ du + λ a dv records reversible bookkeeping; the vertical residual α gives the non‑commutative “area/entropy‑production” measure −α/T. This justifies using macroscopic U, V, T, p to assess irreversibility, independent of microstate enumerations.
