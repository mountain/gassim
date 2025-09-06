# Project Implementation Plan

This plan tracks GasSim from concept to a robust research tool, with phased milestones, key tasks, and verifiable outcomes.

---

## Phase 0: Foundational Architecture — A Validated Hard-Sphere Gas Simulator (NVE Ensemble)

Status: Completed (2025-09-05)

### Goal
Build a physically faithful and computationally robust event-driven simulation core. This serves as the foundation for all subsequent features.

### Key Tasks

#### Rust Core Implementation
- Data Structures:
  - `Particle` (position, velocity, radius, mass, collision_count for invalidation).
  - `Event` (time, event type, collision_count snapshots for invalidation).
- Event Prediction:
  - Analytical collision times for particle–particle (P2P) and particle–wall (P2W) contacts (static walls).
- Event Loop:
  - Priority queue based on `BinaryHeap` with min-time ordering via `Reverse`.
  - Main loop that handles invalidation via collision counters, deterministic progression, and safe time advancement.

#### Python Interface (PyO3)
- Wrap Rust `Simulation` in a Python `GasSim` class.
- Core interfaces:
  - `__new__(num_particles, box_size, radius=1.0, mass=1.0, dim=3, seed=None)`
  - `advance_to(target_time)` (releases the GIL internally)
  - `get_positions()`, `get_velocities()` return NumPy arrays

### Validation Criteria
- Energy Conservation: total kinetic energy constant (within FP error) over long simulations.
- Distribution Check: approach to Maxwell–Boltzmann distribution at equilibrium (qualitative).
- Isotropy: ⟨vx²⟩ ≈ ⟨vy²⟩ ≈ ⟨vz²⟩ after mixing.

### Current Status (2025-09-05)
- Rust core implemented:
  - Event-driven EDMD engine with P2P/P2W prediction and collision invalidation via `collision_count`.
  - Elastic P2P resolution; specular P2W reflection.
  - Full drift to event time; numeric tolerances and bounds safety.
  - Initialization: non-overlapping placement + random velocities; kinetic energy diagnostic.
- Python API (PyO3 0.26, abi3-py313):
  - `GasSim.__new__(num_particles, box_size, radius, mass, dim, seed)`
  - `advance_to(target_time)` (GIL released)
  - `get_positions()`, `get_velocities()`
- Validation:
  - Energy conservation (NVE) passes.
  - Isotropy check passes after mixing.
- Tooling:
  - Toolchain: Python ≥ 3.13; `pyo3 = "0.26"` with `["extension-module","abi3-py313"]`; `numpy = "0.26"`.
  - Build guide: `docs/build.md`
- Code references:
  - Core engine: `src/core/sim.rs`
  - Python bindings: `src/lib.rs`
  - Tests: `tests/phase0.rs`

---

## Phase 1: Introducing Mechanical Work — The Moving Piston

Status: Completed (2025-09-06)

### Goal
Implement externally driven moving walls (piston) enabling mechanical work on the system.

### Key Tasks (as implemented)

#### Rust Core Upgrade
- Dynamic Boundaries:
  - Six axis-aligned walls represented by per-wall kinematics with piecewise-constant velocity (piston model).
- Event Prediction:
  - P2W times computed in the moving-wall frame for both min and max faces along each axis.
- Collision Resolution:
  - Moving mirror reflection in the wall frame; transform back to lab frame.

#### Measurement Module (Rust)
- Impulse and Work:
  - Record per-wall absolute normal impulse |Δp| at each P2W event.
  - Accumulate total work `W` and `work_by_wall[wall_id]`.
  - Maintain `(time, |impulse|)` history for the active piston wall, plus a unified `(time, |impulse|, wall_id)` history across all walls.

#### Python Interface Expansion
- `set_piston(wall_id: int, velocity: float)` sets piecewise-constant piston velocity.
- `get_work_done() -> float`
- `get_pressure_history(window: Optional[float]=None, wall_id: Optional[int]=None) -> np.ndarray`
- `get_mechanical_pressure(window: float, wall_id: int) -> float` (windowed estimate on the active piston wall)

### Validation (2025-09-06)
- `tests/phase1.rs`: moving piston smoke test (work accumulation and/or non-empty impulse history).
- `tests/phase1_5.rs`: mechanical pressure window test (probabilistic; requires impulses in window if pressure unavailable).

---

## Phase 1.5: Diagnostics and δ‑Circuit (Operational Metrics)

Status: Completed (2025-09-06)

### Goal
Provide operational diagnostics for equilibrium quality, anisotropy, and a windowed δ‑circuit residual −α/T.

### Implemented Features
- Temperature tensor:
  - `get_temperature_tensor(mask=None) -> (Txx, Tyy, Tzz, T_scalar, count)` via Rust `temperature_tensor`.
- Velocity histograms and KL diagnostics:
  - `get_velocity_histogram(axes=None, bins=80, range=None, mask=None) -> dict`
  - `get_kl_stats(bins=80, range=None, axes=None, mask=None) -> float` (mean D_KL over requested axes).
- δ‑circuit residual:
  - `get_alpha_over_t(window: float) -> float` using −(ΔU − ΔW)/(T̄ Δt) over the window (Phase 1 adiabatic approximation).
- Work/heat summary:
  - `get_work_heat() -> (W_total, Q_total)`

### Validation (2025-09-06)
- `tests/phase1_5.rs`: temperature tensor finiteness/positivity; histogram shape checks; KL non-negativity; α/T finiteness.

---

## Phase 2: Introducing Heat Exchange — The Thermal Wall

Status: Completed (2025-09-06)

### Goal
Implement thermal walls (heat baths) to support isothermal processes, heat conduction, and First Law validation.

### Key Tasks (as implemented)

#### Rust Core Upgrade
- Thermal walls:
  - Per-wall configuration `(T_bath, accommodation α)` with dedicated RNG streams for reproducibility.
  - Collision: with probability α, resample outgoing velocity from a half-space Maxwellian at `T_bath` in the wall frame (tangential ~ Normal(0, σ²), normal ~ Rayleigh(σ), σ² = T/m); else specular reflect.
- Heat accounting:
  - Event-level split: `dQ = dE − dW`; accumulate `heat_total`, `heat_by_wall[wall_id]`, and append `(time, dQ, wall_id)` to heat history.

#### Python Interface Expansion
- `set_thermal_wall(wall_id: int, T: float, accommodation: float = 1.0) -> None`
- `get_heat_flow() -> float` (total Q)
- `get_heat_history(window: Optional[float]=None, wall_id: Optional[int]=None) -> np.ndarray`
- Per-wall aggregations:
  - `get_work_by_wall() -> np.ndarray`
  - `get_heat_by_wall() -> np.ndarray`
  - `get_impulse_by_wall() -> np.ndarray`

### Validation (2025-09-06)
- `tests/phase2.rs`:
  - Thermal wall heat accumulation with First Law closure (W≈0): `|ΔU − (W+Q)|` below tolerance.
  - Dual thermal walls produce non-empty heat histories with sane event values.
  - Moving thermal piston: First Law closes over windows; heat events include the moving wall id.

---

## Project Status Snapshot (2025-09-06)

- Implemented and validated:
  - Phase 0 (NVE EDMD core, 3D).
  - Phase 1 (moving piston/walls, work/pressure measurement).
  - Phase 1.5 (diagnostics: temperature tensor, histograms, KL, −α/T).
  - Phase 2 (thermal walls with accommodation, heat accounting, First Law).
- Python API (PyO3 0.26, abi3-py313, Python ≥ 3.13) provides:
  - Core: init, advance, get/set positions/velocities, state accessors.
  - Work/pressure: per-wall impulse history, work totals, windowed pressure.
  - Heat: per-wall heat totals and event histories.
  - Diagnostics: temperature tensor, histograms, KL, δ‑circuit residual.
- Build/Tooling:
  - `Cargo.toml`: `pyo3 = "0.26"` with `["extension-module","abi3-py313"]`, `numpy = "0.26"`, `thiserror = "2"`, `rand = "0.9"`.
  - `pyproject.toml`: Python ≥ 3.13; maturin backend.
  - Build guide: `docs/build.md`.
- Tests:
  - `cargo test` validates Phases 0, 1, 1.5, 2.

---

## Phase 3 & Future Outlook

Building on the foundation of the first two phases, expand to more advanced analyses and more complex physical phenomena.

### Local Field Analysis
- Spatial binning (slab/cell) in Rust to measure local fields:
  - Density ρ(x), Velocity u(x), Temperature T(x)

### Complex Phenomena Simulation
- Shock Waves:
  - High-speed piston (Mach M > 1) and local field analysis to resolve shock structure.
- Shear Flow:
  - Lees–Edwards or sliding wall boundary conditions to measure viscosity η.

### Entropy Production and Contact Thermometry
- Python post-processing for entropy production density from local fields.
- Thermo-planes (internal transmissive Maxwell patches) for in-bulk contact thermometry (planned).
- Controller to drive heat flux to zero to estimate contact temperature.
