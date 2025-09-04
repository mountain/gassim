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

### 6.5 Event invalidation

- Each particle maintains `collision_count` (cc) incremented after any realized event it participates in.
- Events store cc snapshots at scheduling time.
- On pop: if a participant’s current cc ≠ stored cc, the event is invalid; discard it.

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

## 8. Numerical Robustness and Determinism

- Tolerances: use combined absolute/relative ε (e.g., 1e-12) for time and geometric checks to avoid duplicate events at t ≈ 0 and grazing contacts.
- Normalization: guard against near-zero-length normals (e.g., if |d| is tiny).
- Stable math: factor common terms, avoid catastrophic cancellation where possible; centralize sensitive operations in math.rs.
- Priority ordering: stable tie-breaking by event kind and ids for determinism.
- Bounds checking: ensure positions remain inside the domain after drift (with tolerance) and adjust if necessary.

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

## 12. Error Handling and Logging

- Errors:
  - Define a comprehensive `Error` enum using `thiserror` (e.g., InvalidParam, MathError, NoCollision, OutOfBounds, NotImplemented).
  - All fallible functions return `Result<T, Error>`; no unwrap/expect in application logic.

- Logging:
  - Integrate `log` facade with optional `env_logger` for debug builds.
  - Consider Python-side toggles to control verbosity.

- Documentation:
  - All public items have Rust doc comments (///) covering purpose, parameters, returns, and error variants, enabling complete `cargo doc`.

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
- Phase 1: compression/expansion with piston; check quasi-static PV^γ constancy (γ=5/3 for monatomic hard spheres) as piston speed→0; for finite rates, entropy change ΔS > 0 qualitatively.
- Phase 2: validate First Law ΔE = W + Q; steady-state conduction with T_L ≠ T_R shows heat flux J_q and negative temperature gradient with preliminary Fourier behavior.

CI gates:
- `cargo fmt --check`
- `cargo clippy -- -D warnings`
- `cargo test`
- `cargo doc --no-deps`

## 14. Reproducibility and Randomness

- Seed management:
  - Initialization accepts an explicit RNG seed; record it in simulation metadata for reproducibility.
  - Thermal wall resampling uses the same RNG stream (thread-safe as needed).

- Metadata:
  - Expose `get_metadata()` returning version, build info, parameters, dimension, and seed.

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

## 19. Appendix: Notation

- r, v: position and velocity vectors (ℝ^D)
- a, m: particle radius and mass
- n: wall outward unit normal
- pos: wall plane position along n (e.g., x = pos)
- R: sum of radii at contact, R = a_i + a_j
- ε: numerical tolerance (absolute/relative)
- Δp, W, Q: impulse, work, and heat accumulators
