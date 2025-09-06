# GasSim: Event-Driven Hard-Sphere Gas Simulator

GasSim is a high-performance event-driven molecular dynamics (EDMD) simulator for a hard-sphere gas. The compute-intensive physics core is implemented in Rust for performance and safety, with a clean Python API (via PyO3) for configuration, control, and analysis.

Highlights:
- Event-driven engine with deterministic priority queue and collision invalidation
- Analytical prediction/resolution of particle–particle and particle–wall collisions
- Moving walls (piston) with mechanical work and pressure impulse measurements
- Thermal walls (heat bath) with Maxwell accommodation and First-Law-consistent heat accounting
- Diagnostics: temperature tensor, velocity histograms, KL divergence, windowed δ‑circuit residual
- Cross-platform Python package

For technical details, see:
- docs/python-api.md — Python API guide with rationale and examples
- docs/design.md — core architecture and data structures
- docs/plan.md — roadmap snapshot and validations

## Installation

Install from PyPI:
```
pip install gassim
```

Requirements:
- Python >= 3.13
- Platform: Linux, macOS, Windows (wheels provided where available)

If you want to build from source or develop locally, see docs/build.md.

## Quick Start (Python)

```
python - <<'PY'
from gassim import GasSim

# Create a 3D box with 64 particles
sim = GasSim(
    num_particles=64,
    box_size=[20.0, 20.0, 20.0],
    radius=0.3,
    mass=1.0,
    dim=3,           # current build supports 3D
    seed=42,
)

# Advance to time t=1.0
sim.advance_to(1.0)

# Access state
pos = sim.get_positions()    # shape (64, 3)
vel = sim.get_velocities()   # shape (64, 3)
print("positions:", pos.shape, "velocities:", vel.shape)

# Configure a moving piston on the x-max wall (wall_id = 1)
sim.advance_to(0.5)
sim.set_piston(1, -0.05)     # slow inward motion
sim.advance_to(2.0)
print("work done:", sim.get_work_done())

# Optional: recent piston impulse events [time, |impulse|]
import numpy as np
hist = np.asarray(sim.get_pressure_history(window=1.0))
print("recent pressure events rows:", hist.shape[0])
PY
```

## Examples

- Half-box relaxation (equilibration and Maxwell speed PDF):
  ```
  python examples/half_box_relax.py
  ```
- Piston compression with macro variable collection:
  ```
  python examples/piston_compression_collect.py
  ```
- Isothermal compression (thermal wall + piston, First Law check):
  ```
  python examples/isothermal_compression.py
  ```

## Development

Local build/testing (see docs/build.md for details):
```
# Rust checks
cargo fmt --all --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo doc --no-deps

# Python API tests
python -m pip install -U maturin pytest numpy
maturin develop -m pyproject.toml
pytest -q tests_py
```

## License

Licensed under the MIT License. See the LICENSE file for details.

Copyright (c) 2025 GasSim Contributors.
