# GasSim: Event-Driven Hard-Sphere Gas Simulator

GasSim is a high-performance event-driven molecular dynamics (EDMD) simulator for a hard-sphere gas. The compute-intensive physics core is implemented in Rust for performance and safety, with a clean Python API (via PyO3 + maturin) for configuration, control, and analysis.

Phase 0 delivers a validated NVE (adiabatic, static box) simulator with:
- Analytical prediction of particle–particle (P2P) and particle–wall (P2W, static) collision times
- Exact elastic collision resolution (hard spheres)
- Deterministic event queue with invalidation via per-particle collision counts
- Minimal Python API: initialize, advance in time, access positions/velocities
- Validation tests: energy conservation and approximate isotropy

For roadmap details, see docs/plan.md and docs/design.md.

Phase 1 adds:
- Moving piston (axis-aligned moving wall with piecewise-constant velocity)
- Mechanical work accumulation on the gas
- Pressure proxy via impulse events on the piston wall
- Python API additions: set_piston(wall_id, velocity), get_work_done(), get_pressure_history(window=None)

## Install and Build (Local Development)

Requirements:
- Rust toolchain (edition 2021)
- Python >= 3.13 (CPython or PyPy)
- maturin >= 1.9

Steps:
1) (Optional, on macOS) ensure a recent Xcode Command Line Tools installation.
2) Install maturin:
   pip install maturin
3) Build and install the Python extension in editable mode:
   maturin develop -m pyproject.toml
This compiles the Rust crate and installs the `gassim` Python module into your current environment.

Format, lint, test (Rust):
```
cargo fmt --all --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo doc --no-deps
```

## Quick Start (Python)

Build and install the extension locally, then run a minimal simulation:

```
pip install maturin
maturin develop -m pyproject.toml

python - <<'PY'
from gassim import GasSim
sim = GasSim(num_particles=64, box_size=[20.0, 20.0, 20.0], radius=0.5, mass=1.0, dim=3, seed=42)
sim.advance_to(1.0)
print("positions:", sim.get_positions().shape, "velocities:", sim.get_velocities().shape)
PY
```

## Example: Piston Compression (Phase 1)

A simple experiment that drives a piston (the x-max wall) inward at a slow, constant speed and queries work and a pressure proxy:

```
python - <<'PY'
from gassim import GasSim
sim = GasSim(num_particles=128, box_size=[20.0, 20.0, 20.0], radius=0.2, mass=1.0, dim=3, seed=1)

# Let system mix briefly with static walls
sim.advance_to(0.5)

# Set piston on x-max wall (axis 0, max side => wall_id = 1) to move inward
sim.set_piston(1, -0.05)

# Advance and query measurements
sim.advance_to(10.0)
print("work:", sim.get_work_done())

# Last 5 time units of piston impulse events (pressure proxy)
hist = sim.get_pressure_history(window=5.0)
print("pressure events (time, |impulse|) rows:", hist.shape)
PY
```

## Example: Half-Box Relaxation (Phase 0)

An end-to-end example that initializes particles in the left half of the box with a Maxwellian velocity distribution, evolves to equilibrium, and compares the empirical speed distribution to the Maxwell prediction:

```
python examples/half_box_relax.py
```

See docs/design.md and docs/plan.md for the architecture and roadmap details.
