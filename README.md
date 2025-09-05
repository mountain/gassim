# GasSim: Event-Driven Hard-Sphere Gas Simulator

GasSim is a high-performance event-driven molecular dynamics (EDMD) simulator for a hard-sphere gas. The compute-intensive physics core is implemented in Rust for performance and safety, with a clean Python API (via PyO3 + maturin) for configuration, control, and analysis.

Phase 0 delivers a validated NVE (adiabatic, static box) simulator with:
- Analytical prediction of particle–particle (P2P) and particle–wall (P2W, static) collision times
- Exact elastic collision resolution (hard spheres)
- Deterministic event queue with invalidation via per-particle collision counts
- Minimal Python API: initialize, advance in time, access positions/velocities
- Validation tests: energy conservation and approximate isotropy

For roadmap details, see docs/plan.md and docs/design.md.

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

## Example: Half-Box Relaxation (Phase 0)

An end-to-end example that initializes particles in the left half of the box with a Maxwellian velocity distribution, evolves to equilibrium, and compares the empirical speed distribution to the Maxwell prediction:

```
python examples/half_box_relax.py
```

See docs/design.md and docs/plan.md for the architecture and roadmap details.
