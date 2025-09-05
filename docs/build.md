# Build and Install (Local)

This guide shows two ways to build and install the `gassim` Python extension locally:

- Option A: maturin develop (requires a virtualenv or conda env)
- Option B: Build a wheel and install it with pip (works without a venv)

Either option installs the Rust extension module so Python can import `gassim`.

## Prerequisites

- Rust toolchain (edition 2021). Install via https://rustup.rs
- Python >= 3.13 (CPython or PyPy recommended)
- pip and setuptools in your target Python environment
- maturin >= 1.9

Install maturin:
```
python -m pip install -U maturin
```

## Option A. maturin develop (requires a virtualenv/conda)

1) Create and activate a virtual environment:
```
# venv example
python -m venv .venv
source .venv/bin/activate
# OR conda
# conda create -n gassim-py313 python=3.13 -y
# conda activate gassim-py313
```

2) Build and install into the active environment:
```
maturin develop -m pyproject.toml
```

3) Test import:
```
python -c "import gassim; from gassim import GasSim; print('GasSim ok')"
```

## Option B. Build a wheel and install it with pip (no venv required)

1) From the repository root (where Cargo.toml is), build a wheel:
```
maturin build --release -m pyproject.toml
```

This writes a `.whl` file into `target/wheels/`.

2) Install the freshly built wheel into your current Python environment:
```
python -m pip install --force-reinstall target/wheels/gassim-*.whl
```

3) Test import:
```
python -c "import gassim; from gassim import GasSim; print('GasSim ok')"
```

## Running the Example

After installing the extension (Option A or B), run:
```
python examples/half_box_relax.py
```

This will:
- Evolve a half-box initialization toward equilibrium
- Save periodic x–y position snapshots to `snapshots/positions_t*.png`
- Save the final speed-vs-Maxwell PDF plot to `snapshots/speed_distribution.png`

## Notes

- If you update Rust code (e.g., `src/core/sim.rs`), you need to reinstall the extension so Python loads the new binary:
  - Option A: re-run `maturin develop -m pyproject.toml` (with the same venv activated)
  - Option B: re-run `maturin build --release -m pyproject.toml` and `python -m pip install --force-reinstall target/wheels/gassim-*.whl`
- Ensure you rebuild/install into the exact Python environment that you use to run the example/script.
- If you hit `maturin develop` errors about not finding a virtualenv, use Option B to build+install a wheel, or create/activate a venv and retry.
