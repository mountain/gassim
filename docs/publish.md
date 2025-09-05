# Publishing GasSim (macOS, universal2) — Simple Guide

This guide shows how to publish GasSim wheels from a macOS laptop as a single universal2 wheel (arm64 + x86_64). It aims to be minimal and reliable for early releases.

If you only need to test locally, you can skip the PyPI steps and just build/install a wheel on your machine.

---

## 1) Prerequisites

- macOS with Xcode Command Line Tools installed (for Clang and SDKs)
  - xcode-select --install
- Python >= 3.13 in your publishing environment
- Rust toolchain installed (via rustup)
- maturin >= 1.9
  - pip install maturin

Optional (recommended):
- Have both Rust targets installed (arm64 and x86_64), so maturin can produce a universal2 wheel:
  - rustup target add aarch64-apple-darwin x86_64-apple-darwin

Note:
- GasSim already targets Python 3.13 via abi3, so one wheel per OS/arch covers all CPython 3.13.x patch versions.

---

## 2) Pre-Publish Checks (Local)

Run standard quality gates before cutting a release:

```
cargo fmt --all --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo doc --no-deps
```

If these pass, proceed.

---

## 3) Versioning

- Update the version in Cargo.toml (and ensure the same version is reflected for Python via pyproject metadata if/when needed).
- Commit changes and optionally tag:
  - git add -A
  - git commit -m "chore(release): v0.x.y"
  - git tag v0.x.y

(You can push the tag later if you maintain a consistent release flow.)

---

## 4) Build a universal2 wheel (Local)

To produce a single universal2 wheel:

```
maturin build --release --universal2 -m pyproject.toml
```

This will generate a .whl under target/wheels/.

Quick local install test (from the repo root):

```
pip install --force-reinstall target/wheels/gassim-*.whl
python -c "import gassim; from gassim import GasSim; print('GasSim ok')"
```

If this works, you’re ready to publish.

---

## 5) (Optional) Test on TestPyPI First

Create an account on https://test.pypi.org/ and generate an API token.

Publish to TestPyPI (replace the token as needed):

```
export MATURIN_PYPI_TOKEN="pypi-AgENdGVzdC5weXBpLm9yZwIk..."
maturin publish --release --universal2 -m pyproject.toml --repository-url https://test.pypi.org/legacy/
```

Install from TestPyPI in a clean venv to verify:

```
python -m venv .venv
source .venv/bin/activate
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple gassim
python -c "import gassim; from gassim import GasSim; print('GasSim from TestPyPI ok')"
deactivate
```

---

## 6) Publish to PyPI (Production)

1) Create an API token at https://pypi.org/ (Account Settings -> API tokens).

2) Export it as an environment variable (do NOT commit this):

```
export MATURIN_PYPI_TOKEN="pypi-AgENdH..."
```

3) Publish:

```
maturin publish --release --universal2 -m pyproject.toml
```

maturin will:
- Build the universal2 wheel
- Upload to PyPI using the token in MATURIN_PYPI_TOKEN

(If your environment uses a shell where `export` differs, set the variable appropriately.)

---

## 7) Post-Publish Sanity Check

In a fresh virtual environment:

```
python -m venv .venv
source .venv/bin/activate
pip install gassim
python - <<'PY'
from gassim import GasSim
sim = GasSim(num_particles=8, box_size=[10.0, 10.0, 10.0], radius=0.5, mass=1.0, dim=3, seed=1)
sim.advance_to(0.1)
print("positions:", sim.get_positions().shape, "velocities:", sim.get_velocities().shape)
PY
deactivate
```

If this runs, the published wheel is installable and functional.

---

## 8) Notes and Tips

- Universal2 builds: maturin `--universal2` automatically builds for both x86_64 and arm64 and combines into a single wheel. Having both Rust targets installed via `rustup target add` helps ensure the build succeeds on either Apple Silicon or Intel Macs.
- macOS deployment target: maturin typically chooses a reasonable default. If you need broader compatibility, you can export `MACOSX_DEPLOYMENT_TARGET`, e.g.:
  - export MACOSX_DEPLOYMENT_TARGET=11.0
- Skip sdist initially for simplicity. You can add `--sdist` later if you want to ship a source distribution as well:
  - maturin publish --release --universal2 --sdist -m pyproject.toml
- Credentials: keep your PyPI token in environment variables or a local keychain. Never commit tokens.
- Re-publish same version: PyPI will reject overwriting existing files. If you need to publish again, bump the version.

That’s it. This process stays on one macOS laptop and produces a single universal2 wheel for Python 3.13, keeping the early release flow simple and reliable.
