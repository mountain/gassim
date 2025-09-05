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

Since your maturin does not expose a `--universal2` flag, build per-architecture wheels and fuse them with `delocate`:

1) Ensure Rust targets are installed:
```
rustup target add aarch64-apple-darwin x86_64-apple-darwin
```

2) Build both wheels into a dedicated directory (wheelhouse). Run these from the repository root (where Cargo.toml is):
```
maturin build --release --target aarch64-apple-darwin -o wheelhouse
maturin build --release --target x86_64-apple-darwin  -o wheelhouse
```

3) Install delocate and fuse the two wheels into a single universal2 wheel:
```
python -m pip install -U delocate

# Identify the two platform-specific wheel filenames under wheelhouse/, then merge:
delocate-merge -w wheelhouse \
  wheelhouse/gassim-<ver>-cp313-abi3-macosx_11_0_arm64.whl \
  wheelhouse/gassim-<ver>-cp313-abi3-macosx_10_12_x86_64.whl
# Output: a new universal2 wheel is written into wheelhouse/ with an automatically determined name (macosx_universal2 tag).
```
Notes:
- Replace `<ver>` with the actual version in the filenames.
- The exact macOS platform tags may vary (e.g., `macosx_12_0_arm64`). Pass the two per-arch wheels you just built; `delocate-merge` writes a new universal2 wheel into `wheelhouse/` with an automatically determined name.

4) Quick local install test of the universal2 wheel (from the repo root):
```
pip install --force-reinstall wheelhouse/gassim-*-macosx_10_12_universal2.whl
python -c "import gassim; from gassim import GasSim; print('GasSim ok')"
```

---

## 5) (Optional) Test on TestPyPI First

Create an account on https://test.pypi.org/ and generate an API token.

Upload the fused universal2 wheel with `twine`:
```
python -m pip install -U twine
export TWINE_USERNAME="__token__"
export TWINE_PASSWORD="pypi-AgENdGVzdC5weXBpLm9yZwIk..."  # your TestPyPI token
twine upload --repository testpypi wheelhouse/gassim-*-macosx_10_12_universal2.whl
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

2) Upload the fused universal2 wheel with `twine`:
```
python -m pip install -U twine
export TWINE_USERNAME="__token__"
export TWINE_PASSWORD="pypi-AgENdH..."  # your PyPI token
twine upload wheelhouse/gassim-*-macosx_10_12_universal2.whl
```

This uploads the already-built universal2 wheel to PyPI without rebuilding.

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

- Universal2 builds: maturin does not provide a built-in `--universal2` flag on your setup. Use the two-arch build + `delocate-merge` method described above to produce a `macosx_universal2` wheel.
- macOS deployment target: maturin typically chooses a reasonable default. If you need broader compatibility, you can export `MACOSX_DEPLOYMENT_TARGET`, e.g.:
  - export MACOSX_DEPLOYMENT_TARGET=11.0
- Source distribution (sdist): optional at this stage. To produce one for upload with `twine`:
  - maturin sdist -m pyproject.toml
- Credentials: keep your PyPI token in environment variables or a local keychain. Never commit tokens.
- Re-publish same version: PyPI will reject overwriting existing files. If you need to publish again, bump the version.

That’s it. This process stays on one macOS laptop and produces a single universal2 wheel for Python 3.13, keeping the early release flow simple and reliable.

## 9) Troubleshooting

- Error: “the manifest-path must be a path to a Cargo.toml file”
  - Cause: running `maturin build` outside the repository root, or passing `-m pyproject.toml` to a subdirectory.
  - Fix: run commands from the repository root (where `Cargo.toml` lives) and omit `-m pyproject.toml` for `build`:
    - maturin build --release --target aarch64-apple-darwin -o wheelhouse
    - maturin build --release --target x86_64-apple-darwin  -o wheelhouse
  - Sanity check: `ls Cargo.toml` should list the manifest in your current directory.

- Error: cannot find Rust target triple
  - Fix: install targets with
    - rustup target add aarch64-apple-darwin x86_64-apple-darwin

- delocate-merge tag mismatch
  - The per-arch wheel tags may differ (e.g., `macosx_12_0_arm64`, `macosx_10_15_x86_64`). Use the exact filenames output by your `maturin build` runs when calling `delocate-merge`.
