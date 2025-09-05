"""
Piston Compression with Macro Variable Collection (Phase 1)

This example:
- Initializes a 3D hard-sphere gas
- Sets a moving piston (x-max wall) with constant inward velocity
- Advances the simulation in time, collecting macroscopic variables:
  * Cumulative mechanical work W(t)
  * Pressure proxy via impulse events on the piston wall
  * Temperature proxy from kinetic energy (assuming k_B = 1 in code units)

Run:
  python examples/piston_compression_collect.py

Requirements:
  - maturin develop -m pyproject.toml (to build/install gassim)
  - numpy (and optionally matplotlib to visualize)
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np

from gassim import GasSim


@dataclass
class Config:
    num_particles: int = 9128
    box_size: Tuple[float, float, float] = (20.0, 20.0, 20.0)  # (Lx, Ly, Lz)
    radius: float = 0.2
    mass: float = 1.0
    dim: int = 3  # fixed in this build
    seed: int = 123

    # Piston (x-max wall => wall_id=1)
    piston_wall_id: int = 1
    piston_velocity: float = -0.05  # inward (negative along +x outward normal)

    # Time control
    t_warmup: float = 0.5  # let system mix with static walls first
    t_final: float = 200.0
    dt: float = 0.2

    # Output
    csv_path: str = "piston_run.csv"


def kinetic_temperature_proxy(vel: np.ndarray, mass: float, kB: float = 1.0) -> float:
    """
    Compute a temperature proxy from velocities using equipartition for an ideal monatomic gas:
      E_k = (D/2) N k_B T
    Return T in code units (k_B=1 by default).
    """
    assert vel.ndim == 2 and vel.shape[1] == 3, "expect velocities shape (N, 3)"
    # total kinetic energy (sum over particles): 0.5 * m * |v|^2
    vsq = np.sum(vel * vel, axis=1)  # (N,)
    Ek = 0.5 * mass * float(np.sum(vsq))
    N, D = vel.shape
    T = (2.0 * Ek) / (D * N * kB)
    return T


def pressure_proxy_from_impulses(impulses: np.ndarray, area: float, window: float) -> float:
    """
    Estimate instantaneous pressure from impulse events over a time window:
      P_inst ≈ (sum |Δp|) / (A * Δt)
    impulses array has shape (M, 2): [time, |impulse|]
    """
    if impulses.size == 0:
        return 0.0
    sum_imp = float(np.sum(impulses[:, 1]))
    if area <= 0 or window <= 0:
        return 0.0
    return sum_imp / (area * window)


def main(cfg: Config) -> None:
    # Construct simulator
    sim = GasSim(
        num_particles=cfg.num_particles,
        box_size=list(cfg.box_size),
        radius=cfg.radius,
        mass=cfg.mass,
        dim=cfg.dim,
        seed=cfg.seed,
    )

    # Derived quantities
    Lx, Ly, Lz = cfg.box_size
    piston_area = Ly * Lz  # area of x-max wall in 3D
    print(f"Initialized: N={cfg.num_particles}, box={cfg.box_size}, area_piston={piston_area}")

    # Warmup with static walls
    sim.advance_to(cfg.t_warmup)
    print(f"Warmup to t={cfg.t_warmup:g}, W={sim.get_work_done():.6g}")

    # Enable piston (x-max wall) with constant inward velocity
    sim.set_piston(cfg.piston_wall_id, cfg.piston_velocity)
    print(
        f"Set piston: wall_id={cfg.piston_wall_id}, velocity={cfg.piston_velocity:+.6g} (piecewise-constant)"
    )

    # Time loop with sampling
    times: List[float] = []
    works: List[float] = []
    temps: List[float] = []
    presses: List[float] = []

    t = cfg.t_warmup
    while t < cfg.t_final - 1e-15:
        t_next = min(cfg.t_final, t + cfg.dt)
        sim.advance_to(t_next)

        # Measurements
        W = sim.get_work_done()
        # Pressure proxy using impulses in last window (Δt)
        impulses = sim.get_pressure_history(window=(t_next - t))
        P = pressure_proxy_from_impulses(impulses, piston_area, (t_next - t))

        # Temperature proxy from velocities (k_B=1)
        v = sim.get_velocities()
        T = kinetic_temperature_proxy(v, cfg.mass, kB=1.0)

        times.append(t_next)
        works.append(W)
        temps.append(T)
        presses.append(P)

        print(
            f"t={t_next:8.3f}  W={W:14.6e}  P_est={P:12.6e}  T_proxy={T:10.6f}  events={impulses.shape[0]}"
        )

        t = t_next

    # Save to CSV
    with open(cfg.csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "work", "pressure_est", "temperature_proxy"])
        for row in zip(times, works, presses, temps):
            w.writerow([f"{row[0]:.9g}", f"{row[1]:.9g}", f"{row[2]:.9g}", f"{row[3]:.9g}"])
    print(f"Wrote CSV: {cfg.csv_path} (rows={len(times)})")

    # Optional: quick plot (requires matplotlib)
    try:
        import matplotlib.pyplot as plt  # type: ignore

        fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
        ax[0].plot(times, works, label="Work W(t)")
        ax[0].set_ylabel("Work")
        ax[0].legend()

        ax[1].plot(times, presses, label="Pressure proxy", color="tab:red")
        ax[1].set_ylabel("Pressure (proxy)")
        ax[1].legend()

        ax[2].plot(times, temps, label="Temperature proxy", color="tab:green")
        ax[2].set_ylabel("T (proxy)")
        ax[2].set_xlabel("time")
        ax[2].legend()

        fig.tight_layout()
        plt.savefig("piston_run.png", dpi=150)
        print("Wrote figure: piston_run.png")
        plt.show()
    except Exception:
        pass


if __name__ == "__main__":
    main(Config())
