# GasSim Design Document

This document specifies the architecture, core data structures, event-driven algorithms, performance and robustness strategies, Python bindings and packaging, testing and validation plan, and phased extensions for GasSim. It serves as the authoritative technical blueprint for building and evolving the simulator.

This document aligns with:
- Roadmap and milestones in docs/plan.md
- Engineering standards in AGENT.md
- Current implementation in Cargo.toml, src/lib.rs, and pyproject.toml

Note: The current Rust module implements Phases 0, 1, 1.5, and 2, exposing a full `GasSim` class via PyO3 with a validated EDMD core, moving piston, diagnostics, and thermal walls.

## 1. Goals and Scope

- Primary goals
  - Implement a physically faithful and computationally robust Event-Driven Molecular Dynamics (EDMD) simulator for a hard-sphere gas (Phase 0: NVE).
  - Keep compute-intensive physics in Rust for performance and memory safety; provide a clean Python API (via PyO3) for configuration, control, and analysis.
  - Extend to a moving piston (mechanical work) and thermal walls (heat exchange) to reproduce canonical non-equilibrium thermodynamic processes and support statistical physics validations.

- Non-goals (current)
  - GPU or distributed execution.
  - Complex intermolecular potentials (focus on hard-sphere collisions).
  - Arbitrary complex geometries (use axis-aligned boxes).
  - Built-in visualization front-end (handled from Python, e.g., Matplotlib).

- Quality and maintainability
  - Follow AGENT.md: rustfmt, clippy (warns as errors), unit and integration tests, complete Rust doc comments, and custom error handling (no unwrap/expect in app logic).
  - Modular, testable design; descriptive naming; single-responsibility functions.

## 2. High-Level Architecture

- Rust core (compute engine)
  - Event prediction: particle–particle (P2P), particle–wall (P2W) for static and moving walls; thermal walls with Maxwell accommodation.
  - Event queue: priority queue over event time (BinaryHeap with min-time ordering via Reverse).
  - Main event loop: deterministic time advancement, event invalidation via collision counters, and robust resolution.
  - Collision resolution and state updates (elastic hard-sphere).
  - Measurements: kinetic energy, momentum impulse (pressure proxy), mechanical work, heat; simple histories.

- Python bindings (control and analysis)
  - Expose class `GasSim` to Python with initialization, advancement to target time, state access as NumPy arrays, and accessors for observables and histories.
  - Use PyO3 and maturin for packaging and distribution; release the GIL during long compute.

- Build and packaging
  - Cargo.toml:
    - `crate-type = ["cdylib","rlib"]`
    - `pyo3 = "0.26"` with features `["extension-module","abi3-py313"]`
    - `numpy = "0.26"`, `thiserror = "2"`, `rand = "0.9"`, `ordered-float = "5"`
  - pyproject.toml:
    - maturin backend, Python `>=3.13`, classifiers for CPython and PyPy, dynamic version.
    - `[tool.maturin] features = ["pyo3/extension-module"]`

### Thermal contact and local thermometry (ADD)

- Core engine supports a thermal wall primitive (implemented) and a future ThermoPlane primitive (planned):
  1) Thermal Wall (Maxwell accommodation at a box face; piston may simultaneously be moving and thermal).
  2) ThermoPlane (an infinitesimally thin internal plane for in‑bulk contact thermometry; transmissive Maxwell “patch” with probability α).

- Observability plane/channel for heat events dQ parallels the work/impulse channel; enables contact‑temperature readout via zero‑heat‑flux condition and closes the δ‑circuit with dU = δQ − p dV. The contact 1‑form dU + p dV − T dS maps to AEG’s α = da − μ du − λ a dv under the dictionary (S,V;U,T,p) = (v, −u; a, λ a, μ), guaranteeing structural consistency.

## 3. Project Layout

- src/
  - lib.rs (PyO3 module entry; binds Rust API to Python)
  - core/
    - particle.rs (Particle type and DIM)
    - event.rs (Event, EventKind, deterministic ordering)
    - sim.rs (Simulation state, event loop, physics, measurements)
    - mod.rs (module wiring)
  - error.rs (custom error enum via thiserror)
- tests/ (Rust integration tests)
  - phase0.rs, phase1.rs, phase1_5.rs, phase2.rs
- docs/
  - design.md (this document)
  - plan.md (roadmap and phases)
  - build.md, publish.md, supplement.md

Planned modules (future):
- measure.rs (structured measurements/observables)
- math.rs (geometry/numerics helpers)
- pq.rs (explicit min-heap wrapper)
- wall.rs (generalized boundary definitions)

## 4. Physical Model and Assumptions

- Hard-sphere gas with elastic collisions (no potential energy).
- Particle attributes: position r, velocity v, radius a, mass m (uniform a and m in current implementation).
- Dimensions: current implementation is fixed 3D (DIM = 3).
- Boundaries by phase:
  - Phase 0: static axis-aligned box with adiabatic walls.
  - Phase 1: moving piston walls with piecewise‑constant normal velocity.
  - Phase 2: thermal walls (heat bath behavior with Maxwell accommodation).
- Units: allow dimensionless or SI-consistent units; configuration records time/length/mass scales in metadata (implicit in parameters).

## 5. Core Data Structures

- Particle
  - id: u32
  - r: [f64; 3]
  - v: [f64; 3]
  - radius: f64
  - mass: f64
  - collision_count: u64 (incremented after each realized event; used to invalidate stale events)

- Event
  - time: f64 (NotNan)
  - kind: enum EventKind { P2P { i: u32, j: u32 }, P2W { i: u32, wall_id: u32 } }
  - cc_i: u64, cc_j: Option<u64> (collision_count snapshots for invalidation)
  - Deterministic ordering: by time, then kind (P2P < P2W), then cc tie-breaker.

- Simulation
  - time_now: f64
  - box_size: [f64; 3] (edge lengths)
  - particles: Vec<Particle>
  - pq: BinaryHeap<Reverse<Event>>
  - Moving walls (Phase 1):
    - wall_min[3], wall_max[3]: kinematics (base_pos, base_time, u)
  - Measurements:
    - work_total: f64; work_by_wall[6]
    - pressure_events: Vec<(t, |Δp|)> (active piston)
    - pressure_events_all: Vec<(t, |Δp|, wall_id)>
    - impulse_by_wall[6]
    - energy_history: Vec<(t, U)>
    - work_history: Vec<(t, W)>
  - Thermal (Phase 2):
    - thermal_walls[6]: Option<(T_bath, accommodation α)>
    - rng_walls: per-wall RNG streams for reproducibility
    - heat_total: f64; heat_by_wall[6]; heat_events: Vec<(t, dQ, wall_id)>

## 6. Event Prediction and Resolution

### 6.1 Particle–Particle (P2P) collision time

- Relative displacement d = r_j − r_i; relative velocity u = v_j − v_i.
- Contact at |d + u t| = R = a_i + a_j.
- Quadratic a t² + b t + c = 0 with:
  - a = u·u, b = 2 d·u, c = d·d − R².
- Conditions:
  - a ≈ 0: no relative motion → no collision.
  - Discriminant Δ < 0: no collision.
  - Candidate earlier root t = (−b − √Δ) / (2a); require b < 0 (approaching) and t > ε.

Resolution (elastic):
- Unit normal n = (r_j − r_i)/|r_j − r_i|.
- Normal relative speed u_n = (v_j − v_i)·n.
- Update normal components only (conserve momentum and kinetic energy):
  - v_i' = v_i + (2 m_j / (m_i + m_j)) u_n n
  - v_j' = v_j − (2 m_i / (m_i + m_j)) u_n n

### 6.2 Particle–Wall (P2W) collision time — static or moving wall (piston)

- Walls indexed: wall_id = 2*k (min), 2*k+1 (max) along axis k in {0,1,2}.
- Moving-wall frame:
  - Min wall: (x − pos_min_now) + (v_x − u_min) t = radius (require v_x − u_min < 0)
  - Max wall: (pos_max_now − x) + (u_max − v_x) t = radius (require u_max − v_x < 0)
- Choose smallest positive t > ε.

Resolution:
- Specular reflection in moving wall frame; transform back:
  - v'_n = 2 u_w − v_n (normal component), tangential unchanged.
- Thermal branch (Phase 2): see §6.4.

### 6.3 Thermal wall (heat bath, Phase 2)

- On collision with a thermal wall at temperature T_b with accommodation α:
  - With probability α, resample outgoing velocity in wall frame from half-space Maxwellian:
    - Tangential components ~ Normal(0, σ²) with σ² = T_b/m
    - Normal (into gas): Rayleigh(σ) with sign into the domain
  - Otherwise, specular reflect in wall frame.
  - Transform back to lab frame by adding wall velocity.
- Heat accounting: dE from speed change; dW = u_w Δp; pure heat dQ = dE − dW.
- Accumulate totals and event histories.

### 6.4 Event invalidation

- Each particle maintains `collision_count`; incremented after realized events.
- Events store cc snapshots at scheduling time.
- On pop: if a participant’s current cc ≠ stored cc, discard (stale).

## 7. Event Loop and Time Advancement

Initialization:
- Rejection-sample non-overlapping initial positions; random velocities in [-1, 1].
- Schedule initial P2P (all pairs, naive) and earliest P2W events.

Advancement (`advance_to(target_time)`):
- While `time_now < target_time`:
  - Pop earliest event; if its time > target_time, drift to target and stop (reinsert event).
  - Validate cc snapshots; discard if invalid.
  - Drift all particles to event time (full drift strategy).
  - Resolve event:
    - P2P: update velocities elastically.
    - P2W: moving mirror or thermal resampling; record impulse/work/heat.
  - Record snapshot: append (t, U) and (t, W).
  - Increment collision_count of involved particles.
  - Reschedule for involved particle(s).

Determinism:
- Priority ordering: by time, then EventKind (P2P < P2W), then cc tie-break.
- EPS_TIME = 1e-12 for time/geometric checks and positivity thresholds.

## 8. Numerical Robustness and Determinism

- Tolerances: combined absolute ε (1e-12) to avoid duplicate/grazing pathologies.
- Normalization: guard against degenerate P2P contact normals.
- Stable math: centralize sensitive operations; avoid catastrophic cancellation.
- Bounds checking: after drift, clamp positions inside the dynamic domain with tolerance.
- Randomness/reproducibility:
  - Explicit seed in `Simulation::new(..., seed)`.
  - Per-wall RNG streams for thermal sampling; deterministic tie-breaking in the queue.

## 9. Performance Design

- Complexity:
  - Naive all-pairs P2P scheduling (O(N² log N) with PQ). Acceptable for functional validation at small/medium N.
- Planned acceleration:
  - Cell lists / uniform grids to prune P2P candidates.
  - Broad-phase via Parry in the future.
- Data layout:
  - AoS (Vec<Particle>) currently; consider SoA for hotspots if scaling up.
- Parallelism:
  - EDMD loop is inherently sequential; parallelize diagnostics/aggregation as feasible.
- Python GIL:
  - Long-running `advance_to` releases GIL via `Python::detach` (PyO3 0.26).

## 10. Measurements and Observables

- Kinetic energy: E_k = Σ ½ m |v|².
- Pressure/work (Phase 1):
  - On piston wall contacts, record |Δp|; accumulate per-wall impulse.
  - Mechanical work: dW = u_w Δp; accumulate `work_total` and `work_by_wall[wall_id]`.
  - Histories:
    - `pressure_events`: (t, |Δp|) for active piston wall.
    - `pressure_events_all`: (t, |Δp|, wall_id) for all walls.
  - Windowed pressure estimate: P ≈ (Σ |Δp|)/(A Δt) over a window for the active piston wall.
- Heat (Phase 2):
  - Pure heat dQ = dE − dW at thermal events; `heat_total`, `heat_by_wall[wall_id]`, `heat_events`.
- Diagnostics (Phase 1.5):
  - Temperature tensor (drift-removed) → (Txx, Tyy, Tzz, T_scalar).
  - Velocity histograms and KL divergence vs Gaussian model (per-axis, centered).
  - δ‑circuit residual −α/T (windowed): −(ΔU − ΔW)/(T̄ Δt) for adiabatic windows.

## 11. Python API (PyO3)

Module name: `gassim`. Expose a `GasSim` class.

Initialization:
- `GasSim(num_particles: int, box_size: Sequence[float], radius: float = 1.0, mass: float = 1.0, dim: int = 3, seed: int | None = None)`
  - Current build supports `dim=3` only; box_size length must equal 3.

Advancement:
- `advance_to(target_time: float) -> None` (releases GIL during compute).

State access:
- `get_positions() -> np.ndarray[float64]  # shape (N, 3)`
- `get_velocities() -> np.ndarray[float64] # shape (N, 3)`
- `set_positions(positions: np.ndarray) -> None` (rebuilds event queue)
- `set_velocities(velocities: np.ndarray) -> None` (rebuilds event queue)

Phase 1 (moving piston):
- `set_piston(wall_id: int, velocity: float) -> None`
- `get_work_done() -> float`
- `get_pressure_history(window: float | None = None, wall_id: int | None = None) -> np.ndarray  # shape (M,2)`
- `get_mechanical_pressure(window: float, wall_id: int) -> float`

Phase 1.5 (diagnostics & δ‑circuit):
- `get_temperature_tensor(mask: Optional[np.ndarray]) -> (Txx, Tyy, Tzz, T_scalar, count)`
- `get_velocity_histogram(axes: Optional[List[int]] = None, bins: int = 80, range: Optional[Tuple[float,float]] = None, mask: Optional[np.ndarray] = None) -> Dict[int, Dict[str, np.ndarray]]`
- `get_kl_stats(bins: int = 80, range: Optional[Tuple[float,float]] = None, axes: Optional[List[int]] = None, mask: Optional[np.ndarray] = None) -> float`
- `get_alpha_over_t(window: float) -> float`
- `get_work_heat() -> Tuple[float, float]`

Phase 2 (thermal walls):
- `set_thermal_wall(wall_id: int, T: float, accommodation: float = 1.0) -> None`
- `get_heat_flow() -> float`
- `get_heat_history(window: Optional[float] = None, wall_id: Optional[int] = None) -> np.ndarray  # shape (M,3)`
- `get_work_by_wall() -> np.ndarray[float64] # shape (6,)`
- `get_heat_by_wall() -> np.ndarray[float64] # shape (6,)`
- `get_impulse_by_wall() -> np.ndarray[float64] # shape (6,)`

Implementation notes:
- Zero-copy arrays via `numpy` crate where applicable; otherwise copy on export.
- Error mapping: Rust `Error` → `PyErr` with clear messages for invalid parameters or state errors.

### Thermal contact / thermometry (Phase 2+) (ADD)

- Planned:
  - `set_thermo_plane(plane_id: int, point: Sequence[float], normal: Sequence[float], T: float, accommodation: float = 0.02, area: Optional[float] = None) -> None`
  - `get_heat_history(boundary_id: int, window: float) -> np.ndarray  # [t, dQ]`
  - Controller: `find_contact_temperature(boundary_id, tol, max_iter)` drives ⟨ĖQ⟩ → 0.

### δ‑circuit helpers (ADD)

- `get_mechanical_pressure(window: float, wall_id: int) -> float` (extended to arbitrary walls)
- `get_entropy_IG(time_series: bool = False) -> (S_vals or ΔS)`
- `get_alpha_over_T(window: float) -> float  # windowed −α/T`
- `get_work_heat() -> (W_total, Q_total)`

## 12. Error Handling and Logging

- Errors:
  - `thiserror`-based `Error` enum (e.g., InvalidParam, MathError, NoCollision, OutOfBounds, NotImplemented, DiagnosticInsufficientSamples).
  - All fallible functions return `Result<T, Error>`; no unwrap/expect in application logic.

- Logging:
  - Integrate `log` facade with optional `env_logger` for debug builds (planned toggles from Python).

- Documentation:
  - Public items have Rust doc comments; `cargo doc` builds API docs.

## 13. Testing and Validation

Rust unit tests:
- event.rs: deterministic ordering and cc validation.
- particle.rs: invariants and kinetic energy.

Rust integration tests (tests/):
- Phase 0: energy conservation (long NVE runs), isotropy after mixing.
- Phase 1: moving piston smoke test; impulse/work sanity.
- Phase 1.5: temperature tensor, histograms, KL diagnostic, windowed α/T.
- Phase 2: First Law closure ΔU = W + Q; dual thermal walls heat histories; moving thermal piston closure over windows.

CI gates (recommended):
- `cargo fmt --check`
- `cargo clippy -- -D warnings`
- `cargo test`
- `cargo doc --no-deps`

Additional tests (ADD):
- Thermal sampler moments vs theory (half-space Maxwell, Rayleigh for normal).
- Contact temperature recovery via thermo-plane controller.
- δ‑circuit “must‑pass” ensemble checks across piston speeds and end volumes.

## 14. Reproducibility and Randomness

- Seed management:
  - Initialization accepts an explicit RNG seed; record in metadata (parameters) for reproducibility.
  - Thermal wall resampling uses per-wall RNG streams; deterministic within a given configuration.

- Metadata:
  - Expose via Python-side helpers (planned) version, build info, parameters, dimension, and seed.

## 15. Versioning and Packaging

- Semantic versioning:
  - 0.x may include breaking changes; attempt to keep Python API stable where feasible.
- Packaging:
  - Built and published via maturin.
  - Python ≥ 3.13 (abi3-py313); CPython and PyPy supported per classifiers.

## 16. Milestones and Alignment with docs/plan.md

- Phase 0 (Foundational EDMD, NVE): Completed and validated.
- Phase 1 (Mechanical work — moving piston): Completed; impulse/work histories; pressure estimate.
- Phase 1.5 (Diagnostics & δ‑circuit): Completed; temperature tensor, histograms/KL, α/T.
- Phase 2 (Heat exchange — thermal wall): Completed; First Law closure and histories.

Planned:
- Local field analysis (slab/cell binning) and complex flows (shock, shear).
- Thermo-planes and contact thermometry; extended δ‑circuit helpers.

## 17. Risks and Mitigations

- Numerical pathologies (grazing collisions, t≈0):
  - Tolerances, tie-breaking, guard against repeated immediate events; normalized contact normals.
- Event storms and stale queues:
  - Partial/full rescheduling heuristics (planned); current collision-count invalidation.
- API evolution:
  - Track changes in this design doc; provide migration notes for Python/Rust APIs.
- Performance at scale:
  - Introduce cell lists/broad-phase; consider SoA layouts.

## 18. Future Enhancements

- Local field analysis:
  - Spatial binning for ρ(r), u(r), T(r); export compact diagnostics and histories.
- Complex phenomena:
  - Shock waves via high-speed piston; shear flows via sliding walls or Lees–Edwards BCs.
- Parallelism:
  - Parallel measurement aggregation; multi-instance ensembles.
- Visualization and tooling:
  - Exporters to NumPy/Parquet/NetCDF; example notebooks.

### Additional local fields and geometry (ADD)

- Binning and contact‑temperature field estimation (kriging/IDW over sparse thermo‑planes).
- Variable accommodation, moving thermo‑planes (Lagrangian probes), and specular/roughness models.

## 19. Appendix: Notation

- r, v: position and velocity vectors (ℝ³)
- a, m: particle radius and mass
- n: wall outward unit normal
- R: sum of radii at contact, R = a_i + a_j
- ε: numerical tolerance (absolute/relative)
- Δp, W, Q: impulse, work, and heat accumulators

## 20. Appendix: Maxwell sampling cookbook (ADD)

- Tangential components: Box–Muller for Normal(0, σ²).
- Positive normal (flux‑weighted): sample R ~ Rayleigh(σ), set c⁺_n = ±R into the domain.
- Specular mixing: with probability 1−α, reflect c⁺_n = −c⁻_n, keep tangentials.
- Heat increment: ΔQ = (m/2)(||v⁺||² − ||v⁻||²).
- Controller (Python‑side):

```python
# every Δt window
qdot = sum(dQ) / Δt
if abs(qdot) < tol:
    done
T_b -= kappa * qdot  # integral control; or bracket root by bisection
```

### Notes on Theory Alignment

- The simulator’s δ‑circuit is the numerical instantiation of the contact 1‑form α_TD = dU + p dV − T dS; under the AEG dictionary (S,V;U,T,p) = (v, −u; a, λ a, μ), it pulls back to α = da − μ du − λ a dv. The horizontal part δa = μ du + λ a dv records reversible bookkeeping; the vertical residual α gives the non‑commutative “area/entropy‑production” measure −α/T, enabling macroscopic irreversibility assessment independent of microstate enumeration.
