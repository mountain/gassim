# GasSim: a gas simulation tool

This project is a high-performance particle simulation platform for research on non-equilibrium thermodynamic processes.

It employs event-driven molecular dynamics to accurately simulate a hard-sphere gas, providing a controllable and measurable microscopic laboratory for studying complex physical phenomena.

## Core Philosophy

This project combines the exceptional performance and memory safety of Rust with the flexibility, ease of use, and rich scientific computing ecosystem of Python.

### Computational Core (Rust)

All compute-intensive physical calculations—including event prediction, collision resolution, and state updates—are handled in a compiled Rust module to ensure maximum performance.

### Control and Analysis (Python)

Users interact with the simulation via Python scripts to configure parameters, control the simulation’s progress, run experiments, and perform data analysis and visualization, improving productivity and flexibility.

## Planned Features

This simulator aims to reproduce several canonical non-equilibrium thermodynamic processes, including:
- Mechanical processes: finite-rate adiabatic/isothermal compression and expansion
- Thermodynamic processes: steady-state heat conduction and heat exchange in isothermal conditions
- Fluid-dynamics phenomena: piston-driven shock waves and shear flows
- Statistical-physics validation: tests of non-equilibrium relations such as the Jarzynski equality and the Crooks fluctuation theorem

## Tech Stack

- Core engine: Rust
- Python bindings: PyO3
- Build and packaging: Maturin
- Collision detection (broad phase): Parry (planned)
- Scientific computing and visualization: NumPy, Matplotlib

## Getting Started (Planned)

### Install the Python package

```bash
pip install gassim
```

### Run a simulation script

```bash
python run_simulation.py
```

### Example: run_simulation.py

```python
from thermosim_ed import GasSim
import matplotlib.pyplot as plt

# 1. Initialize the simulator
sim = GasSim(num_particles=500, box_size=[10.0, 10.0])

# 2. Evolve the system
sim.advance_to(target_time=50.0)

# 3. Get results and visualize
positions = sim.get_positions()  # Returns a NumPy array
plt.scatter(positions[:, 0], positions[:, 1])
plt.show()
```

## Project Status

This project is under active development. For the detailed development roadmap, please see docs/plan.md.
