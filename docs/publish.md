# Publishing GasSim (macOS) — Practical Guide

This guide shows how to publish GasSim wheels from a macOS laptop. It provides a simple, tested approach that works reliably for early releases.

---

## 1) Prerequisites

- macOS with Xcode Command Line Tools installed
  - xcode-select --install
- Python >= 3.13 (use `python3` command)
- Rust toolchain installed (via https://rustup.rs)
- maturin >= 1.9
  - python3 -m pip install -U maturin

Note: GasSim targets Python 3.13 via abi3, so one wheel covers all CPython 3.13.x patch versions.

---

## 2) Pre-Publish Checks (Local)

Run quality gates before cutting a release:

```
cargo fmt --all --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo doc --no-deps
```

All must pass before proceeding.

---

## 3) Versioning

Update the version in Cargo.toml, commit changes, and optionally tag:
```
# Edit Cargo.toml version field
git add -A
git commit -m "chore(release): v0.x.y"
git tag v0.x.y
```

---

## 4) Build and Test Wheel (Local)

1) Build a release wheel:
```
maturin build --release
```

The wheel will be written to `target/wheels/gassim-VERSION-cp313-abi3-macosx_ARCH.whl`.

2) Test the wheel locally (important!):
```
python3 -m pip install --force-reinstall target/wheels/gassim-*.whl
python3 -c "import gassim; from gassim import GasSim; print('GasSim ok')"
```

3) Run the example to ensure it works:
```
python3 examples/half_box_relax.py
```

If this completes without errors and prints results like "Final internal energy U = ...", the wheel is ready.

---

## 5) (Optional) Test on TestPyPI First

1) Create an account at https://test.pypi.org/ and generate an API token.

2) Upload to TestPyPI:
```
python3 -m pip install -U twine
export TWINE_USERNAME="__token__"
export TWINE_PASSWORD="pypi-AgENdGVzdC5weXBpLm9yZwIk..."  # your TestPyPI token
twine upload --repository testpypi target/wheels/gassim-*.whl
```

3) Test install from TestPyPI in a clean venv:
```
python3 -m venv test_venv
source test_venv/bin/activate
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple gassim
python3 -c "import gassim; print('TestPyPI install works')"
deactivate
rm -rf test_venv
```

---

## 6) Publish to PyPI (Production)

1) Create an API token at https://pypi.org/ (Account Settings -> API tokens).

2) Upload to PyPI:
```
export TWINE_USERNAME="__token__"
export TWINE_PASSWORD="pypi-AgENdH..."  # your PyPI production token
twine upload target/wheels/gassim-*.whl
```

---

## 7) Post-Publish Verification

Test in a fresh environment:
```
python3 -m venv verify_venv
source verify_venv/bin/activate
pip install gassim
python3 -c "from gassim import GasSim; sim = GasSim(num_particles=8, box_size=[10.0, 10.0, 10.0], radius=0.5, mass=1.0, dim=3, seed=1); sim.advance_to(0.1); print('Published wheel works!')"
deactivate
rm -rf verify_venv
```

---

## 8) Advanced: Universal2 Wheels (Optional)

For broader compatibility across Intel and Apple Silicon Macs:

1) Install additional Rust targets:
```
rustup target add aarch64-apple-darwin x86_64-apple-darwin
```

2) Build both architectures:
```
mkdir wheelhouse
maturin build --release --target aarch64-apple-darwin -o wheelhouse
maturin build --release --target x86_64-apple-darwin -o wheelhouse
```

3) Merge using delocate:
```
python3 -m pip install -U delocate
ls wheelhouse/  # Check actual filenames
delocate-merge -w wheelhouse wheelhouse/gassim-*-arm64.whl wheelhouse/gassim-*-x86_64.whl
```

4) Upload the universal2 wheel:
```
twine upload wheelhouse/gassim-*-universal2.whl
```

---

## 9) Important Notes

- **Always rebuild and reinstall after code changes**: Python will cache the old extension module. Use `--force-reinstall` to ensure the latest version is loaded.
- **Test thoroughly**: Run the example script after each wheel build/install to catch issues early.
- **Keep it simple**: Single-architecture wheels work fine for most use cases. Only use universal2 if you need broad compatibility.
- **Version bumps**: PyPI won't accept re-uploads of the same version. Bump the version in Cargo.toml for each release.

## 10) Troubleshooting

- **"out of bounds" errors**: Ensure you've rebuilt and reinstalled the wheel after recent code changes.
- **Import errors**: Check that you're using the correct Python environment where you installed the wheel.
- **Build errors**: Verify you're running maturin from the repository root (where Cargo.toml is located).
- **Platform tag mismatches**: Use `ls target/wheels/` or `ls wheelhouse/` to see actual filenames when using delocate-merge.

This process keeps releases simple and reliable while maintaining high quality through testing at each step.
