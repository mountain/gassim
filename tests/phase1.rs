use gassim::core::Simulation;

/// Phase 1 smoke test: moving piston performs mechanical work and records impulses.
/// This is a probabilistic test; it asserts that after enabling a moving wall and
/// advancing sufficiently in time, either the cumulative work changes or we have
/// recorded at least one piston impulse event.
#[test]
fn moving_piston_work_and_pressure_smoke() -> gassim::error::Result<()> {
    // Setup a moderately dense gas to encourage wall collisions
    let mut sim = Simulation::new(64, [20.0, 20.0, 20.0], 0.3, 1.0, Some(24680))?;
    // Let system run a bit with static walls
    sim.advance_to(0.5)?;

    // Choose x-max wall (axis 0 max => wall_id = 1) as piston; move inward slowly
    sim.set_piston_velocity(1, -0.05)?;
    let w0 = sim.work_total();

    // Advance long enough to allow several potential collisions
    sim.advance_to(10.0)?;
    let w1 = sim.work_total();

    // Check measurements are sane
    assert!(w1.is_finite(), "work_total must remain finite");

    // Either work changed or pressure events were recorded
    let events = sim.pressure_events();
    assert!(
        w1 != w0 || !events.is_empty(),
        "expected nonzero work change or at least one piston impulse event (events: {})",
        events.len()
    );

    Ok(())
}
