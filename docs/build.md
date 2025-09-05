# Build, Test, and Package Guide

This document describes how to build, lint, test, and generate docs for GasSim across environments.

The project is a Rust core library with a Python extension module exposed via PyO3 and built with Maturin.

## Prerequisites

- Rust (stable toolchain)
- Python >= 3.13 (CPython or PyPy)
- Maturin (for building the Python extension)
  - Install: `pip install maturin`

## Dependencies and Versions

- PyO3: 0.26
- numpy crate: 0.23

These are aligned in `Cargo.toml` to avoid FFI link conflicts. The crate is built as an extension module targeting the stable Python ABI (abi3) for Python 3.13.

## Important Environment Variables

None required for this project.

With PyO3 ≥ 0.26 and the `abi3-py313` feature enabled, PyO3 automatically handles newer interpreter versions (e.g., Python 3.13) without manual environment variables. When the active Python interpreter is newer than PyO3’s maximum supported version, PyO3 falls back to an abi3 build for the latest supported version.

## Standard Rust Workflow

This is the recommended sequence to ensure conformance with AGENT.md:

- Format:
  - `cargo fmt --all -- --check`
- Lint (fail on warnings):
  - `cargo clippy -- -D warnings`
- Test:
  - `cargo test`
- Docs:
  - `cargo doc --no-deps`

## Python Extension Build (Maturin)

- Develop (in-place build, importable from the venv):
  - `maturin develop`
- Build wheel:
  - `maturin build`
- Build and publish (example):
  - `maturin publish`

macOS note:
- For universal2 wheels on Apple Silicon/Intel, see Maturin’s docs (e.g., `--universal2` flag).
- Ensure you’re building inside a Python virtual environment matching your target.

## Example: Local CI Script

Below is a bash script you can paste into a local CI step or run manually.

```bash
#!/usr/bin/env bash
set -euo pipefail

echo "[1/4] Formatting check"
cargo fmt --all -- --check

echo "[2/4] Lint (clippy -D warnings)"
cargo clippy -- -D warnings

echo "[3/4] Tests"
cargo test

echo "[4/4] Docs"
cargo doc --no-deps

echo "Done."
```

## Troubleshooting

- Link conflicts or version mismatches:
  - Ensure `Cargo.toml` aligns PyO3 (`0.26`) and numpy crate (`0.23`).
  - Confirm the crate features include `extension-module` and `abi3-py313`.
- Importing from Python fails after `maturin develop`:
  - Verify you’re in the same Python environment where `maturin develop` was run.
  - On macOS, ensure correct architecture (e.g., not mixing arm64 vs x86_64).
- Dense initial configurations:
  - For very dense initial configurations (large N, large radius, small box), initialization can fail. Use fewer particles, a smaller radius, or a larger box.

## What We Run in CI (Equivalent)

- `cargo fmt --all -- --check`
- `cargo clippy -- -D warnings`
- `cargo test`
- `cargo doc --no-deps`
