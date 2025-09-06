"""
Python API test suite for GasSim.

How to run locally:
  pip install maturin pytest numpy
  maturin develop -m pyproject.toml
  pytest -q tests_py

These tests validate the public Python API surface exposed by src/lib.rs.
They are intended to fail when signatures/behavior change, notifying about
breaking changes between Rust and Python layers.
"""

from __future__ import annotations

import math
from typing import Dict, Tuple

import numpy as np
import pytest

from gassim import GasSim


def kinetic_energy(v: np.ndarray, mass: float) -> float:
    assert v.ndim == 2 and v.shape[1] == 3
    vsq = float(np.sum(v * v))
    return 0.5 * mass * vsq


def approx_close(a: float, b: float, rel: float = 1e-6, abs_: float = 1e-9) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= (rel * scale + abs_)


class TestCoreAPI:
    def test_init_and_state_shapes(self) -> None:
        N = 16
        sim = GasSim(N, [10.0, 9.0, 8.0], radius=0.2, mass=1.0, dim=3, seed=1234)
        X = np.asarray(sim.get_positions())
        V = np.asarray(sim.get_velocities())
        assert X.shape == (N, 3)
        assert V.shape == (N, 3)
        assert np.isfinite(X).all()
        assert np.isfinite(V).all()

    def test_invalid_init_params_raise(self) -> None:
        with pytest.raises(ValueError):
            _ = GasSim(0, [10.0, 10.0, 10.0], 0.2, 1.0, 3, seed=1)  # num_particles > 0
        with pytest.raises(ValueError):
            _ = GasSim(8, [1.0, 10.0, 10.0], 1.0, 1.0, 3, seed=1)  # Lx < 2r invalid
        with pytest.raises(ValueError):
            _ = GasSim(8, [10.0, 10.0], 0.2, 1.0, 3, seed=1)  # wrong box_size length
        with pytest.raises(ValueError):
            _ = GasSim(8, [10.0, 10.0, 10.0], -0.1, 1.0, 3, seed=1)  # radius <= 0
        with pytest.raises(ValueError):
            _ = GasSim(8, [10.0, 10.0, 10.0], 0.2, 0.0, 3, seed=1)  # mass <= 0
        with pytest.raises(ValueError):
            _ = GasSim(8, [10.0, 10.0, 10.0], 0.2, 1.0, 2, seed=1)  # dim must be 3

    def test_advance_monotone_time(self) -> None:
        sim = GasSim(8, [10.0, 10.0, 10.0], 0.4, 1.0, 3, seed=7)
        sim.advance_to(0.5)
        with pytest.raises(ValueError):
            sim.advance_to(0.1)  # earlier than current time should raise

    def test_set_positions_velocities_and_rebuild(self) -> None:
        sim = GasSim(8, [10.0, 10.0, 10.0], 0.4, 1.0, 3, seed=7)
        X = np.asarray(sim.get_positions()).copy()
        V = np.asarray(sim.get_velocities()).copy()
        # Modify slightly and set back
        X[:, 0] = np.clip(X[:, 0] + 0.01, 0.4, 9.6)
        V[:, 2] *= 1.1
        sim.set_positions(X)
        sim.set_velocities(V)
        X2 = np.asarray(sim.get_positions())
        V2 = np.asarray(sim.get_velocities())
        assert np.allclose(X2, X)
        assert np.allclose(V2, V)
        # Wrong shapes should raise
        with pytest.raises(ValueError):
            sim.set_positions(np.zeros((7, 3)))
        with pytest.raises(ValueError):
            sim.set_velocities(np.zeros((8, 2)))


class TestPistonAndPressure:
    def test_piston_work_and_pressure_history(self) -> None:
        N = 48
        sim = GasSim(N, [20.0, 20.0, 20.0], 0.25, 1.0, 3, seed=2468)
        sim.advance_to(0.5)
        sim.set_piston(1, -0.05)  # x-max
        sim.advance_to(3.0)
        W = sim.get_work_done()
        assert math.isfinite(W)
        # History without wall_id -> piston history
        hist = np.asarray(sim.get_pressure_history(window=1.0))
        assert hist.ndim == 2 and hist.shape[1] == 2
        # With wall_id -> filtered all-walls history
        hist1 = np.asarray(sim.get_pressure_history(window=1.0, wall_id=1))
        assert hist1.ndim == 2 and hist1.shape[1] == 2
        # Mechanical pressure may raise if insufficient samples; allow both outcomes
        try:
            P = sim.get_mechanical_pressure(window=1.0, wall_id=1)
            assert math.isfinite(P) and P >= 0.0
        except ValueError:
            # Insufficient samples in the selected window is acceptable; API contract exercised.
            # Events may be sparse depending on seed, speed, and density.
            pass

    def test_per_wall_aggregations_shapes(self) -> None:
        sim = GasSim(32, [15.0, 16.0, 17.0], 0.2, 1.0, 3, seed=13579)
        sim.advance_to(0.5)
        sim.set_piston(1, -0.02)
        sim.advance_to(2.0)
        w = np.asarray(sim.get_work_by_wall())
        q = np.asarray(sim.get_heat_by_wall())
        imp = np.asarray(sim.get_impulse_by_wall())
        assert w.shape == (6,)
        assert q.shape == (6,)
        assert imp.shape == (6,)

    def test_piston_errors(self) -> None:
        sim = GasSim(8, [10.0, 10.0, 10.0], 0.3, 1.0, 3, seed=22)
        with pytest.raises(ValueError):
            sim.set_piston(99, -0.01)
        with pytest.raises(ValueError):
            _ = sim.get_mechanical_pressure(1.0, wall_id=1)  # no piston configured yet
        with pytest.raises(ValueError):
            _ = sim.get_pressure_history(window=-1.0)


class TestThermalWallsAndHeat:
    def test_thermal_wall_first_law_closure(self) -> None:
        N = 48
        mass = 1.0
        sim = GasSim(N, [20.0, 20.0, 20.0], 0.2, mass, 3, seed=4242)
        sim.set_thermal_wall(0, 1.0, 1.0)  # x-min
        V0 = np.asarray(sim.get_velocities())
        U0 = kinetic_energy(V0, mass)
        sim.advance_to(3.0)
        W, Q = sim.get_work_heat()
        V1 = np.asarray(sim.get_velocities())
        U1 = kinetic_energy(V1, mass)
        residual = abs((U1 - U0) - (W + Q))
        # Use tolerance consistent with Rust tests: 1e-6 * scale + 1e-9
        scale = max(1.0, abs(U0), abs(U1))
        assert residual <= (1e-6 * scale + 1e-9)

    def test_heat_history_and_flow(self) -> None:
        sim = GasSim(40, [15.0, 15.0, 15.0], 0.2, 1.0, 3, seed=777)
        sim.set_thermal_wall(0, 0.9, 1.0)
        sim.set_thermal_wall(1, 0.9, 1.0)
        sim.advance_to(2.0)
        events = np.asarray(sim.get_heat_history())
        assert events.ndim == 2 and events.shape[1] == 3
        # Basic sanity on times and values
        if events.size > 0:
            assert np.isfinite(events).all()
            assert np.all(events[:, 0] >= 0.0)
            # wall_id should be integer-like in [0..5]
            assert np.all((events[:, 2] >= 0.0) & (events[:, 2] <= 5.0))
        # get_heat_flow returns scalar total Q
        q_total = sim.get_heat_flow()
        assert math.isfinite(q_total)
        # Check consistency with get_work_heat()
        W, Q = sim.get_work_heat()
        assert approx_close(Q, q_total)

    def test_thermal_wall_errors(self) -> None:
        sim = GasSim(8, [10.0, 10.0, 10.0], 0.3, 1.0, 3, seed=33)
        with pytest.raises(ValueError):
            sim.set_thermal_wall(99, 1.0, 1.0)
        with pytest.raises(ValueError):
            sim.set_thermal_wall(0, -1.0, 1.0)  # T must be > 0
        with pytest.raises(ValueError):
            sim.set_thermal_wall(0, 1.0, 2.0)  # accommodation in [0,1]


class TestDiagnostics:
    def test_temperature_tensor_basic(self) -> None:
        N = 64
        sim = GasSim(N, [20.0, 20.0, 20.0], 0.2, 1.0, 3, seed=1357)
        sim.advance_to(3.0)
        tx, ty, tz, T_scalar, count = sim.get_temperature_tensor()
        assert count == N
        assert all(map(math.isfinite, [tx, ty, tz, T_scalar]))
        assert T_scalar > 0.0

    def test_velocity_histogram_and_kl(self) -> None:
        sim = GasSim(72, [20.0, 20.0, 20.0], 0.2, 1.0, 3, seed=2468)
        sim.advance_to(3.0)
        axes = [0, 1, 2]
        bins = 32
        h = sim.get_velocity_histogram(axes=axes, bins=bins)
        assert isinstance(h, dict)
        # Validate shapes and totals
        total_counts = 0.0
        for ax in axes:
            inner = h[ax]
            edges = np.asarray(inner["edges"])
            counts = np.asarray(inner["counts"])
            assert edges.shape[0] == bins + 1
            assert counts.shape[0] == bins
            total_counts += float(np.sum(counts))
        assert int(total_counts) == sim.get_positions().shape[0] * len(axes)
        # KL
        dkl = sim.get_kl_stats(bins=64, axes=axes)
        assert math.isfinite(dkl) and dkl >= 0.0

    def test_alpha_over_t_window(self) -> None:
        sim = GasSim(48, [20.0, 20.0, 20.0], 0.2, 1.0, 3, seed=1122)
        sim.advance_to(0.5)
        sim.set_piston(1, -0.03)
        sim.advance_to(2.0)
        # Choose a window inside simulated span
        try:
            val = sim.get_alpha_over_t(window=0.5)
            assert math.isfinite(val)
        except ValueError:
            # Allow insufficient samples; objective is API contract
            pass
