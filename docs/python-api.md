# GasSim Python API: Rationale, Design, and Examples

This document explains how to use the GasSim Python API exposed by the Rust event‑driven molecular dynamics (EDMD) engine. It covers the rationale behind key design choices, the current API surface, error semantics, performance characteristics, and end‑to‑end examples aligned with the latest implementation.

- Module/package name: `gassim` (PyO3)
- Current capabilities: 3D hard‑sphere EDMD core, moving piston (mechanical work), diagnostics (temperature tensor, histograms, KL, δ‑circuit residual), thermal walls (pure heat exchange) with First‑Law‑consistent accounting.

See also:
- docs/plan.md for phased roadmap and status
- docs/design.md for engine architecture and data structures
- examples/ for runnable Python examples

---

## 1) Rationale and Design Principles

- Keep physics hot paths in Rust; expose a thin, clean Python layer for configuration, control, and analysis.
- Deterministic event ordering and robust invalidation via collision counters.
- Windowed, low‑overhead measurement channels (impulse/pressure, work, heat, diagnostics).
- Reproducibility: explicit seed for initialization and per‑wall RNG streams for thermal sampling.
- Python GIL released during long `advance_to` calls for responsiveness in Python environments.
- Simple, explicit APIs with clear shapes/types for NumPy data.

Current simplifications:
- Fixed 3D (DIM=3) to focus on stability and validations.
- Axis‑aligned rectangular box; six walls indexed 0..5 (see §3).
- Piecewise‑constant piston velocities and per‑wall thermal baths (Maxwell accommodation).

---

## 2) Installation and Requirements

- Python: >= 3.13 (abi3‑py313)
- Build backend: maturin
- Crates: pyo3 = "0.26" with features ["extension‑module","abi3‑py313"], numpy = "0.26", thiserror = "2", rand = "0.9", ordered‑float = "5"

Typical local dev flow:
- `pip install maturin`
- `maturin develop` (in project root) to build and install `gassim` into the active environment

---

## 3) Coordinate System and Wall Indexing

- 3D box: axes k ∈ {0,1,2} correspond to x, y, z.
- Wall indices: for axis k,
  - 2k   = min face (x_min, y_min, z_min)
  - 2k+1 = max face (x_max, y_max, z_max)
- Valid `wall_id` range: 0..5
- Piston may be configured on any wall (`set_piston`); thermal baths can be set per wall (`set_thermal_wall`).

---

## 4) Class and Methods

Class: `gassim.GasSim`

Initialization:
- `GasSim(num_particles: int, box_size: Sequence[float], radius: float = 1.0, mass: float = 1.0, dim: int = 3, seed: int | None = None)`
  - dim must be 3 (current build). `box_size` must be length 3 with each L > 2*radius.
  - Particles initialized non‑overlapping by rejection sampling; velocities uniform in [-1,1] per component.
  - If `seed` is provided, initialization and per‑wall RNG streams are reproducible.

Time advancement:
- `advance_to(target_time: float) -> None`
  - Advances to an absolute time (monotone, ≥ current time). Releases the GIL internally.

State access:
- `get_positions() -> np.ndarray[float64]` shape (N, 3)
- `get_velocities() -> np.ndarray[float64]` shape (N, 3)
- `set_positions(positions: np.ndarray[float64]) -> None` shape (N, 3)
  - Rebuilds event queue; caller must ensure positions are finite, within domain, and non‑overlapping.
- `set_velocities(velocities: np.ndarray[float64]) -> None` shape (N, 3)
  - Rebuilds event queue; values must be finite.

Mechanical work and pressure:
- `set_piston(wall_id: int, velocity: float) -> None`
  - Sets piecewise‑constant normal velocity from the current time; re‑bases wall position at call time. Marks this wall as the active piston for `pressure_events`.
- `get_work_done() -> float`
  - Total mechanical work W on gas since start (sum over all walls).
- `get_pressure_history(window: float | None = None, wall_id: int | None = None) -> np.ndarray`
  - If `wall_id` is None: returns piston history, shape (M, 2) rows [time, |impulse|].
  - If `wall_id` is set: returns all‑walls history filtered by the given wall, same shape.
  - `window` limits to events with time ≥ (now−window).
- `get_mechanical_pressure(window: float, wall_id: int) -> float`
  - Windowed estimate using impulse sum over (area * window). Currently limited to the active piston wall.

Thermal walls and heat:
- `set_thermal_wall(wall_id: int, T: float, accommodation: float = 1.0) -> None`
  - Maxwell accommodation α ∈ [0,1]; α=1 means pure Maxwellian resampling in the wall frame, else mix with specular reflection.
- `get_heat_flow() -> float`
  - Total accumulated pure heat Q (sum of dQ = dE − dW over thermal events).
- `get_heat_history(window: float | None = None, wall_id: int | None = None) -> np.ndarray`
  - Returns shape (M, 3) rows [time, dQ, wall_id], optionally filtered by window/wall.

Per‑wall aggregations:
- `get_work_by_wall() -> np.ndarray[float64]` shape (6,)
- `get_heat_by_wall() -> np.ndarray[float64]` shape (6,)
- `get_impulse_by_wall() -> np.ndarray[float64]` shape (6,)

Diagnostics (Phase 1.5):
- `get_temperature_tensor(mask: Optional[np.ndarray[bool]]) -> (Txx, Tyy, Tzz, T_scalar, count)`
  - Drift‑removed; returns per‑axis components and scalar T, plus sample count.
- `get_velocity_histogram(axes: Optional[List[int]] = None, bins: int = 80, range: Optional[Tuple[float,float]] = None, mask: Optional[np.ndarray[bool]] = None) -> Dict[int, Dict[str, np.ndarray]]`
  - Returns per‑axis dict: `{axis_index: {"edges": np.ndarray, "counts": np.ndarray}}`.
- `get_kl_stats(bins: int = 80, range: Optional[Tuple[float,float]] = None, axes: Optional[List[int]] = None, mask: Optional[np.ndarray[bool]] = None) -> float`
  - Mean KL divergence between centered component histograms and Gaussian model (σ² = T/m).
- `get_alpha_over_t(window: float) -> float`
  - Windowed δ‑circuit residual −α/T ≈ −(ΔU − ΔW)/(T̄ Δt) for adiabatic windows.

---

## 5) Error Semantics

- All Rust errors are mapped to Python `ValueError` with descriptive messages.
- Common cases:
  - Invalid parameters (e.g., `wall_id` out of range, non‑positive T, accommodation not in [0,1]).
  - Time monotonicity violation for `advance_to`.
  - Shape/type mismatches for arrays in set/get operations.
  - DiagnosticInsufficientSamples (e.g., not enough events or samples in the selected window or mask) → `ValueError` with a specific message.
- Best practice: design runs to ensure enough collisions and use windows that contain adequate data.

---

## 6) Performance Notes

- `advance_to` releases the GIL internally for long computations.
- Event loop uses collision‑count invalidation and deterministic ordering (P2P < P2W at equal time).
- Histories are append‑only vectors; for large/long runs consider down‑sampling externally or advancing in segments.

---

## 7) Reproducibility

- `seed` in `GasSim.__init__` seeds initialization and per‑wall RNG streams (thermal).
- For thermal walls, stochastic resampling is still reproducible given the same seed and sequence of events.

---

## 8) Examples

All examples assume:
```python
import numpy as np
from gassim import GasSim
```

1) NVE baseline (Phase 0)
```python
# Initialize a 3D box, 64 particles
sim = GasSim(
    num_particles=64,
    box_size=[20.0, 20.0, 20.0],
    radius=0.2,
    mass=1.0,
    dim=3,            # must be 3 in this build
    seed=12345,
)

# Advance in steps to allow many events
for t in np.linspace(0.5, 50.0, 50):
    sim.advance_to(float(t))

# Inspect positions/velocities
X = sim.get_positions()     # shape (64, 3)
V = sim.get_velocities()    # shape (64, 3)

# Basic invariants (user-side)
ke = 0.5 * 1.0 * np.sum(V**2)   # kinetic energy
print("KE ~", ke)
```

2) Adiabatic piston compression (Phase 1)
```python
sim = GasSim(96, [30.0, 30.0, 30.0], radius=0.25, mass=1.0, dim=3, seed=98765)

# Mix with static walls
sim.advance_to(0.5)

# Make x_max (wall_id=1) the piston, moving inward slowly
sim.set_piston(1, velocity=-0.05)

# Advance
sim.advance_to(5.0)

W = sim.get_work_done()
print("Work done on gas:", W)

# Pressure history on piston
hist = sim.get_pressure_history(window=2.0)  # shape (M,2) [t, |Δp|]
if hist.shape[0] > 0:
    P = sim.get_mechanical_pressure(window=2.0, wall_id=1)
    print("Windowed mechanical pressure (wall 1):", P)
```

3) Thermal wall (Phase 2) and First Law closure
```python
sim = GasSim(64, [20.0, 20.0, 20.0], radius=0.2, mass=1.0, dim=3, seed=424242)

# Configure x_min wall (0) as thermal at T=1.0, full accommodation
sim.set_thermal_wall(0, T=1.0, accommodation=1.0)

# Record U0 and (W,Q) snapshots
def kinetic_energy(V, m=1.0): return 0.5 * m * np.sum(V**2)
U0 = kinetic_energy(sim.get_velocities())

W0, Q0 = sim.get_work_heat()
sim.advance_to(5.0)
W1, Q1 = sim.get_work_heat()
U1 = kinetic_energy(sim.get_velocities())

dU = U1 - U0
dW = W1 - W0
dQ = Q1 - Q0
print("First-Law residual |ΔU - (ΔW + ΔQ)| =", abs(dU - (dW + dQ)))

# Heat history (array of [t, dQ, wall_id])
events = sim.get_heat_history()
print("Heat events:", events.shape[0])
```

4) Diagnostics and δ‑circuit residual (Phase 1.5)
```python
sim = GasSim(80, [20, 20, 20], 0.2, 1.0, 3, seed=24680)
sim.advance_to(10.0)

# Temperature tensor
tx, ty, tz, T, count = sim.get_temperature_tensor()
print("Txx, Tyy, Tzz, T =", tx, ty, tz, T, "count =", count)

# Velocity histograms and KL divergence
hists = sim.get_velocity_histogram(axes=[0,1,2], bins=64)
dkl = sim.get_kl_stats(bins=64, axes=[0,1,2])
print("D_KL(components vs Gaussian) =", dkl)

# δ‑circuit residual over a window (adiabatic piston window suggested)
sim.set_piston(1, velocity=-0.03)
sim.advance_to(3.0)
val = sim.get_alpha_over_t(window=1.0)
print("Windowed −α/T =", val)
```

5) Per‑wall aggregation helpers
```python
sim = GasSim(64, [20, 20, 20], 0.2, 1.0, 3, seed=1122)
sim.set_piston(1, -0.02)         # x_max piston
sim.set_thermal_wall(0, 1.0, 1.0) # x_min thermal
sim.advance_to(4.0)

work_by_wall   = sim.get_work_by_wall()     # shape (6,)
heat_by_wall   = sim.get_heat_by_wall()     # shape (6,)
impulse_by_wall= sim.get_impulse_by_wall()  # shape (6,)

# Indices: [x_min, x_max, y_min, y_max, z_min, z_max]
print("W per wall:", work_by_wall)
print("Q per wall:", heat_by_wall)
print("|Δp| per wall:", impulse_by_wall)
```

6) Setting positions/velocities (advanced)
```python
sim = GasSim(8, [10.0, 10.0, 10.0], 0.4, 1.0, 3, seed=7)

# Overwrite state (ensure finite values, within domain, non-overlapping!)
X = sim.get_positions()
V = sim.get_velocities()
X[:, 0] = np.clip(X[:, 0], 0.4, 9.6)   # keep away from walls by radius
V[:, 2] *= 0.5                         # scale z-velocity

sim.set_positions(X)  # rebuilds event queue
sim.set_velocities(V) # rebuilds event queue

sim.advance_to(1.0)
```

---

## 9) Tips and Gotchas

- Dimension: must pass `dim=3`. `box_size` must be length 3 and each L > 2*radius.
- Windows: ensure `window` fits within the simulated time span and contains enough events; otherwise some diagnostics may raise `ValueError` for insufficient samples.
- Piston pressure: `get_mechanical_pressure` is currently limited to the active piston wall (`wall_id` passed to `set_piston`).
- Histories:
  - `get_pressure_history(window=None)` without `wall_id` returns the piston history (shape (M,2)).
  - Per‑wall impulse totals are available via `get_impulse_by_wall()` (shape (6,)).
- Thermal walls: accommodation ∈ [0,1]. With α=1, the wall behaves as an ideal Maxwellian bath in its frame; the normal component is Rayleigh‑distributed (flux‑weighted).
- Reproducibility: provide a `seed`. For thermal runs, results are still stochastic but reproducible given identical seeds and event sequences.
- Large runs: consider chunking `advance_to` calls and periodically extracting/compacting histories on the Python side.

---

## 10) Theory Alignment and First Law

- Mechanical work increment: dW = u_w Δp (wall frame→lab frame consistent accounting).
- Pure heat increment at thermal events: dQ = dE − dW.
- First Law closure:
  - Over any window [t0, t1], ΔU = Σ(dW) + Σ(dQ) holds (validated in tests, including moving thermal piston).
- δ‑circuit residual (operational): −α/T ≈ −(ΔU − ΔW)/(T̄ Δt) over adiabatic windows (Phase 1.5 diagnostic).

---

## 11) Cross‑References

- examples/half_box_relax.py
- examples/piston_compression_collect.py
- examples/isothermal_compression.py
- docs/design.md (architecture)
- docs/plan.md (phases and validations)
