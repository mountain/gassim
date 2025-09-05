#!/usr/bin/env python3
"""
Half-box relaxation example for GasSim (Phase 0)

- Initializes particles confined to the left half of a 3D box (no overlaps).
- Initializes velocities from a Maxwellian at temperature T0 (kB = 1 units).
- Evolves the system using the event-driven core to time T_end.
- Reports internal energy U, inferred temperature T, and ideal-gas pressure p_th = 2U/(3V).
- Compares empirical speed distribution with Maxwell speed PDF.

Usage:
  1) Build and install the extension locally:
     pip install maturin
     maturin develop -m pyproject.toml

  2) Run this script:
     python examples/half_box_relax.py
"""

import math
import sys
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt
import os

try:
    from gassim import GasSim
except Exception as e:
    print("Failed to import gassim. Did you run 'maturin develop'?", file=sys.stderr)
    raise


def rejection_sample_half_box(
    N: int, box_size: Tuple[float, float, float], radius: float, rng: np.random.Generator
) -> np.ndarray:
    """
    Return (N,3) positions in the left half of the box:
      x in [radius, Lx/2 - radius], y in [radius, Ly - radius], z in [radius, Lz - radius]
    ensuring hard-sphere non-overlap with minimal rejection sampling (O(N^2)).
    """
    Lx, Ly, Lz = box_size
    rmin = np.array([radius, radius, radius], dtype=np.float64)
    rmax = np.array([Lx / 2.0 - radius, Ly - radius, Lz - radius], dtype=np.float64)
    if np.any(rmax <= rmin):
        raise ValueError("Box too small for half-box placement at given radius")

    pos = np.zeros((N, 3), dtype=np.float64)
    min_sq = (2.0 * radius) ** 2
    for i in range(N):
        attempts = 0
        while True:
            attempts += 1
            if attempts > 1_000_000:
                raise RuntimeError("Failed to place particle without overlap")
            candidate = rmin + rng.random(3) * (rmax - rmin)
            # Check overlap only against already placed
            if i == 0:
                pos[i] = candidate
                break
            d2 = np.sum((pos[:i] - candidate) ** 2, axis=1)
            if np.all(d2 >= min_sq):
                pos[i] = candidate
                break
    return pos


def maxwellian_velocities(N: int, mass: float, T0: float, rng: np.random.Generator) -> np.ndarray:
    """
    Sample velocities from a 3D Maxwellian with kB = 1:
      Each component ~ N(0, sigma^2) with sigma^2 = T0 / mass.
    Remove center-of-mass drift to avoid net momentum.
    """
    sigma = math.sqrt(T0 / mass)
    v = rng.normal(loc=0.0, scale=sigma, size=(N, 3))
    v -= v.mean(axis=0, keepdims=True)
    return v


def main():
    # Simulation parameters
    N = 256
    Lx, Ly, Lz = 50.0, 40.0, 30.0
    radius = 0.5
    mass = 1.0
    T0 = 1.0  # initial temperature (kB=1 units)
    seed = 20250905
    box_size = (Lx, Ly, Lz)

    # Snapshot configuration
    out_dir = "snapshots"
    snapshot_stride = 20  # save a PNG every 20 steps
    os.makedirs(out_dir, exist_ok=True)

    # Construct simulator via Python API
    sim = GasSim(
        num_particles=N,
        box_size=[Lx, Ly, Lz],
        radius=radius,
        mass=mass,
        dim=3,
        seed=seed,
    )

    rng = np.random.default_rng(seed)

    # Overwrite initial positions to be left-half-only and non-overlapping
    positions = rejection_sample_half_box(N, box_size, radius, rng)
    sim.set_positions(positions)

    # Overwrite initial velocities with Maxwellian at T0 (remove COM drift)
    velocities = maxwellian_velocities(N, mass, T0, rng)
    sim.set_velocities(velocities)

    # Evolve the system to equilibrium
    T_end = 200.0
    steps = 400
    for s in range(1, steps + 1):
        t = T_end * s / steps
        sim.advance_to(t)

        # Save periodic position snapshots (x-y) as PNGs
        if s % snapshot_stride == 0:
            R_snap = sim.get_positions()
            fig_s, ax_s = plt.subplots(figsize=(5.0, 4.0))
            ax_s.scatter(R_snap[:, 0], R_snap[:, 1], s=6, alpha=0.5)
            ax_s.set_title(f"Positions (x-y) at t={t:.3f}")
            ax_s.set_xlabel("x")
            ax_s.set_ylabel("y")
            ax_s.set_xlim(radius, Lx - radius)
            ax_s.set_ylim(radius, Ly - radius)
            ax_s.set_aspect("equal", adjustable="box")
            fig_s.tight_layout()
            fig_s.savefig(os.path.join(out_dir, f"positions_t{t:08.3f}.png"), dpi=150)
            plt.close(fig_s)

    # Gather results
    R = sim.get_positions()
    V = sim.get_velocities()

    # Internal energy U, inferred temperature T (kB = 1): U = (3/2) N T for monatomic ideal gas
    U = 0.5 * mass * np.sum(V * V)
    T_inf = (2.0 * U) / (3.0 * N)

    # Ideal-gas pressure prediction p_th = 2U / (3V)
    V_box = Lx * Ly * Lz
    p_th = 2.0 * U / (3.0 * V_box)

    print(f"N = {N}, box = ({Lx}, {Ly}, {Lz}), radius = {radius}, mass = {mass}")
    print(f"Final internal energy U = {U:.6f}")
    print(f"Inferred temperature T = {T_inf:.6f}")
    print(f"Ideal-gas pressure prediction p_th = {p_th:.6f}")

    # Empirical speed distribution vs. Maxwell speed PDF
    speeds = np.linalg.norm(V, axis=1)
    s_max = np.percentile(speeds, 99.5)
    bins = np.linspace(0.0, s_max, 60)
    hist, edges = np.histogram(speeds, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Maxwell speed distribution in 3D:
    # f(v) = 4π (m / (2π T))^(3/2) v^2 exp(- m v^2 / (2 T))
    m_over_2T = mass / (2.0 * T_inf)
    pref = 4.0 * math.pi * (m_over_2T / math.pi) ** 1.5
    f_maxwell = pref * (centers**2) * np.exp(-m_over_2T * (centers**2))

    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

    ax[0].scatter(R[:, 0], R[:, 1], s=6, alpha=0.5)
    ax[0].set_title("Particle positions (x-y)")
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")
    ax[0].set_xlim(radius, Lx - radius)
    ax[0].set_ylim(radius, Ly - radius)
    ax[0].set_aspect("equal", adjustable="box")

    ax[1].plot(centers, hist, label="Empirical speed PDF", drawstyle="steps-mid")
    ax[1].plot(centers, f_maxwell, label="Maxwell speed PDF", linestyle="--")
    ax[1].set_xlabel("speed")
    ax[1].set_ylabel("PDF")
    ax[1].set_title("Speed distribution")
    ax[1].legend()

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "speed_distribution.png"), dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
