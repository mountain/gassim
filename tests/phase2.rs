use gassim::core::Simulation;
use gassim::error::Result;

/// Phase 2: Thermal wall should accumulate heat and satisfy First Law ΔU = W + Q.
/// In this test we do not move the piston, so W ≈ 0 and ΔU ≈ Q.
#[test]
fn thermal_wall_heat_accumulates_and_first_law() -> Result<()> {
    let mut sim = Simulation::new(64, [20.0, 20.0, 20.0], 0.2, 1.0, Some(424242))?;
    // Set x-min wall (wall_id = 0) as thermal at T=1.0, full accommodation.
    sim.set_thermal_wall(0, 1.0, 1.0)?;

    // Initial internal energy
    let u0 = sim.kinetic_energy();

    // Advance to allow many wall interactions
    sim.advance_to(5.0)?;

    // Work and heat
    let (w, q) = sim.work_heat();
    let u1 = sim.kinetic_energy();
    let delta_u = u1 - u0;

    // First Law residual
    let residual = (delta_u - (w + q)).abs();

    // Tolerance scaled to energy magnitude (loose to account for stochastic path / ordering)
    let scale = u0.abs().max(u1.abs()).max(1.0);
    let tol = 1e-6 * scale + 1e-9;

    assert!(
        residual <= tol,
        "First Law residual too large: |ΔU - (W+Q)| = {residual}, tol = {tol}, ΔU = {delta_u}, W = {w}, Q = {q}"
    );

    // We expect at least one heat event recorded
    let heat_events = sim.heat_events();
    assert!(
        !heat_events.is_empty(),
        "expected at least one heat event when a thermal wall is configured"
    );

    Ok(())
}

/// Phase 2: Heat history should be populated when multiple thermal walls are active.
#[test]
fn thermal_walls_heat_history_present() -> Result<()> {
    let mut sim = Simulation::new(48, [15.0, 15.0, 15.0], 0.2, 1.0, Some(777))?;
    // Set both x-min (0) and x-max (1) walls as thermal with the same T
    sim.set_thermal_wall(0, 0.8, 1.0)?;
    sim.set_thermal_wall(1, 0.8, 1.0)?;

    sim.advance_to(3.0)?;

    let events = sim.heat_events();
    assert!(
        !events.is_empty(),
        "expected some heat events with dual thermal walls"
    );

    // Basic sanity checks on event values
    for (t, dq, wall_id) in events {
        assert!(
            t.is_finite() && t >= 0.0,
            "event time must be finite and non-negative"
        );
        assert!(dq.is_finite(), "heat increment must be finite");
        assert!(
            wall_id == 0
                || wall_id == 1
                || wall_id == 2
                || wall_id == 3
                || wall_id == 4
                || wall_id == 5,
            "wall_id out of expected range"
        );
    }

    Ok(())
}
