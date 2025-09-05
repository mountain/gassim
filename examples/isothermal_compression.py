#!/usr/bin/env python3
"""
Phase 2 Example: Isothermal Compression

- x-min wall (wall_id=0) is a thermal bath at temperature T_bath
- x-max wall (wall_id=1) acts as a slowly moving piston compressing the gas inward
- Demonstrates:
  * set_thermal_wall(wall_id, T, accommodation)
  * set_piston(wall_id, velocity)
  * advance_to(t)
  * get_work_heat() -> (W_total, Q_total)
  * get_heat_history(window=..., wall_id=...) for recent heat events
  * First Law consistency: ΔU ≈ W + Q

Run:
  python examples/isothermal_compression.py
"""

import numpy as np

from gassim import GasSim


def kinetic_energy(vel: np.ndarray, mass: float) -> float:
    # vel shape (N, 3), mass scalar
    vsq = (vel * vel).sum(axis=1).sum()
    return 0.5 * mass * vsq


def main():
    # Configuration
    N = 128
    box = [20.0, 20.0, 20.0]
    radius = 0.2
    mass = 1.0
    seed = 2025

    # Thermal bath and piston
    thermal_wall_id = 0  # x-min
    T_bath = 1.0
    accommodation = 1.0  # full Maxwell accommodation
    piston_wall_id = 1  # x-max
    piston_speed = -0.02  # slow inward motion (negative for x-max moving inward)
    t_mix = 0.5  # pre-mix time with static walls
    t_final = 8.0

    # Initialize simulation
    sim = GasSim(
        num_particles=N,
        box_size=box,
        radius=radius,
        mass=mass,
        dim=3,
        seed=seed,
    )

    # Let system mix briefly
    sim.advance_to(t_mix)

    # Configure thermal wall and piston
    sim.set_thermal_wall(thermal_wall_id, T_bath, accommodation)
    sim.set_piston(piston_wall_id, piston_speed)

    # Record initial internal energy from velocities
    v0 = np.asarray(sim.get_velocities())
    U0 = kinetic_energy(v0, mass)

    # Advance to final time
    sim.advance_to(t_final)

    # Collect results
    W_total, Q_total = sim.get_work_heat()
    v1 = np.asarray(sim.get_velocities())
    U1 = kinetic_energy(v1, mass)
    dU = U1 - U0

    # Heat history in the last window
    window = 2.0
    heat_hist = np.asarray(sim.get_heat_history(window=window))

    # First Law residual
    residual = abs(dU - (W_total + Q_total))
    scale = max(1.0, abs(U0), abs(U1))
    rel_residual = residual / scale

    # Output summary
    print("=== Isothermal Compression (Phase 2) ===")
    print(f"N = {N}, box = {box}, radius = {radius}, mass = {mass}, seed = {seed}")
    print(f"Thermal wall: id = {thermal_wall_id}, T = {T_bath}, accommodation = {accommodation}")
    print(f"Piston: wall_id = {piston_wall_id}, speed = {piston_speed}")
    print(f"Times: t_mix = {t_mix}, t_final = {t_final}")
    print()
    print(f"U0 = {U0:.6f}")
    print(f"U1 = {U1:.6f}")
    print(f"W_total = {W_total:.6f}")
    print(f"Q_total = {Q_total:.6f}")
    print(f"ΔU = {dU:.6f}")
    print(f"First Law residual |ΔU - (W+Q)| = {residual:.6e}  (relative: {rel_residual:.6e})")
    if heat_hist.size > 0:
        print(f"Recent heat events in last {window} time units: {heat_hist.shape[0]} rows [t, dQ, wall_id]")
        # Show a small tail
        tail = heat_hist[-min(5, heat_hist.shape[0]) :]
        for row in tail:
            t, dq, wid = row.tolist()
            print(f"  t={t:.3f}, dQ={dq:.6f}, wall_id={int(wid)}")
    else:
        print(f"No heat events recorded in the last {window} time units.")

    # Suggestion: for stricter validation runs, reduce |piston_speed| further toward the quasi-static limit.


if __name__ == "__main__":
    main()
