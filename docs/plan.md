# Project Implementation Plan

This plan aims to build GasSim from a concept into a powerful research tool.

The project is divided into several progressive phases, each with clear goals, key tasks, and verifiable outcomes.

---

## Phase 0: Foundational Architecture — A Validated Hard-Sphere Gas Simulator (NVE Ensemble)

Status: Completed (2025-09-05)

### Goal
Build a physically faithful and computationally robust event-driven simulation core. This serves as the foundation for all subsequent features.

### Key Tasks

#### Rust Core Implementation
- Data Structures:
  - Define `Particle` (position, velocity, radius, mass).
  - Define `Event` (time, event type, and a `collision_count` for invalidation).
- Event Prediction:
  - Implement precise analytical calculations for collision times between particles and between particles and fixed walls.
- Event Loop:
  - Implement a priority queue based on `BinaryHeap`.
  - Build the main `tick()` loop, ensuring correct event handling, time progression, and management of invalidated events.

#### Python Interface (PyO3)
- Wrap the Rust `Simulation` struct into a `GasSim` Python class.
- Provide basic interfaces:
  - `__new__(num_particles, box_size, ...)`: Initialize the simulation environment.
  - `advance_to(time)`: Evolve the simulation to a specific point in time.
  - `get_positions()`, `get_velocities()`: Return the system state as NumPy arrays.

### Validation Criteria
- Energy Conservation: Total kinetic energy remains constant (within floating-point error) over long simulations.
- Distribution Check: At equilibrium, the particle velocity distribution conforms to the Maxwell–Boltzmann distribution.
- Isotropy: At equilibrium, average kinetic energy is equal in all directions (⟨vx²⟩ ≈ ⟨vy²⟩ ≈ ⟨vz²⟩).

### Current Status (2025-09-05)
- Status: Completed.
- Rust core implemented:
  - Event-driven MD engine with priority queue, deterministic ordering, collision invalidation via `collision_count`.
  - Analytical P2P and static P2W collision prediction; elastic P2P resolution; specular P2W reflection.
  - Lazy drift to event time; numeric tolerances and bounds safety.
  - Initialization via non-overlapping placement and random velocities; kinetic energy diagnostic.
- Python API (PyO3 ≥ 0.26, abi3-py313) implemented:
  - `GasSim.__new__(num_particles, box_size, radius, mass, dim, seed)`
  - `advance_to(target_time)` (releases GIL via `Python::detach`)
  - `get_positions()`, `get_velocities()` returning NumPy arrays (numpy crate 0.26).
- Validation:
  - Energy conservation test passes (NVE).
  - Isotropy check passes after mixing.
  - Unit tests for particles, events, and simulation math pass.
- Tooling/infra:
  - Toolchain upgraded to Python ≥ 3.13; `pyo3 = "0.26"` with `["extension-module","abi3-py313"]`; `numpy = "0.26"`.
  - No environment variables required for build/test.
  - CI (GitHub Actions) runs fmt, clippy (warnings as errors), tests, and docs across Ubuntu/macOS/Windows.
- Code references:
  - Core engine: `src/core/sim.rs`
  - Python bindings: `src/lib.rs`
  - Integration tests: `tests/phase0.rs`
  - Build guide: `docs/BUILD.md`
---

## Phase 1: Introducing Mechanical Work — The Moving Piston

### Goal
Implement an externally driven piston with a predefined trajectory, enabling the simulator to perform work on the system.

### Key Tasks

#### Rust Core Upgrade
- Dynamic Boundaries:
  - Refactor `box_size` into a configurable list of walls.
  - Each `Wall` object contains geometry (normal vector `n`), position, and wall velocity `u_w(t)` (initially piecewise constant).
- Update Event Prediction:
  - Implement collision time with a moving wall: `t_hit = (L0 - x0) / (v_x - u_w)`.
- Update Collision Resolution:
  - Implement the velocity update for collisions with a moving wall (moving mirror reflection).

#### Measurement Module (Rust)
- Impulse Accumulation:
  - In the `tick()` loop, record the impulse Δp transferred to the piston upon each particle collision.
- Work and Pressure Calculation:
  - Within a specified time window Δt, compute instantaneous pressure `P_inst` and cumulative work `W` done on the gas by the piston based on accumulated impulse.

#### Python Interface Expansion
- Allow defining the wall’s motion trajectory during `GasSim` initialization (e.g., specify one wall as a piston with velocity `u_w`).
- Add methods:
  - `get_work_done()`
  - `get_pressure_history()`

### First Target Experiment & Validation
- Finite-Rate Adiabatic Compression/Expansion:
  - Set up a piston moving at constant velocity with all walls adiabatic.
- Validation:
  - In the quasi-static limit (`u_w → 0`), verify that `P V^γ = const` holds (for a monatomic hard-sphere gas, `γ = 5/3`).
  - For rapid, finite-rate compression, verify that the system’s entropy change `ΔS > 0` after the process.

---

## Phase 2: Introducing Heat Exchange — The Thermal Wall

### Goal
Implement a thermal wall (heat bath) capable of exchanging energy with the system, enabling the simulation of isothermal and other thermodynamic processes.

### Key Tasks

#### Rust Core Upgrade
- Wall Types:
  - Add a wall-type field (e.g., `Adiabatic`, `HeatBath`) and a temperature attribute to the `Wall` struct.
- Thermal Reflection Logic:
  - Implement thermal wall collisions: when a particle hits a thermal wall, outgoing velocity is resampled from a half-space Maxwellian distribution corresponding to the wall temperature `T_b` instead of mirror reflection.

#### Measurement Module (Rust)
- Heat Flow Calculation:
  - On each thermal wall collision, compute kinetic energy change `ΔE = E_after - E_before` and accumulate into total heat flow `Q`.

#### Python Interface Expansion
- Allow one or more walls to be set as thermal baths with specified temperatures from Python.
- Add a `get_heat_flow()` method.

### Target Experiments & Validation
- Isothermal Compression/Expansion:
  - Set the wall opposite the piston as a thermal bath and drive the piston to compress the gas.
  - Validation:
    - Record total `ΔE`, `W`, and `Q` and verify the First Law of Thermodynamics: `ΔE = W + Q`.
    - Accumulate data for later verification of the Jarzynski equality.
- Steady-State Heat Conduction:
  - Set left and right walls as thermal baths with different temperatures, `T_L` and `T_R`.
  - Validation:
    - Once a non-equilibrium steady state is reached, measure heat flux `J_q` through the system.
    - Using spatial binning, measure the temperature gradient `∇T` and perform preliminary verification of Fourier’s Law, `J_q ∝ -∇T`.

---

## Phase 3 & Future Outlook

Building on the foundation of the first two phases, expand to more advanced analyses and more complex physical phenomena.

### Local Field Analysis
- Implement spatial binning (Slab Analysis) in the Rust core to measure local fields:
  - Density `ρ(x)`
  - Velocity `u(x)`
  - Temperature `T(x)`

### Complex Phenomena Simulation
- Shock Waves:
  - Drive the piston at high speeds (Mach number `M > 1`) and use local field analysis to observe shock structure.
- Shear Flow:
  - Implement Lees–Edwards or sliding wall boundary conditions to measure fluid shear viscosity `η`.

### Entropy Production Calculation
- Utilize Python to post-process local field data returned from Rust.
- Estimate entropy production density using continuum mechanics formulas, connecting microscopic simulation with macroscopic irreversibility.
