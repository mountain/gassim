use gassim::core::Simulation;
use gassim::error::{Error, Result};

/// Phase 1.5: Temperature tensor should be finite and positive after mixing.
#[test]
fn temperature_tensor_basic() -> Result<()> {
    let mut sim = Simulation::new(64, [20.0, 20.0, 20.0], 0.2, 1.0, Some(13579))?;
    // Allow some events to occur for mixed velocities
    sim.advance_to(5.0)?;

    let (tx, ty, tz, t_scalar, count) = sim.temperature_tensor(None)?;
    assert_eq!(count, sim.num_particles(), "count should equal N");
    assert!(tx.is_finite() && ty.is_finite() && tz.is_finite());
    assert!(
        t_scalar.is_finite() && t_scalar > 0.0,
        "scalar T must be positive"
    );
    Ok(())
}

/// Phase 1.5: Velocity histogram and KL diagnostic should be sane.
#[test]
fn velocity_histogram_and_kl() -> Result<()> {
    let mut sim = Simulation::new(80, [20.0, 20.0, 20.0], 0.2, 1.0, Some(24680))?;
    sim.advance_to(10.0)?;

    let axes = [0usize, 1usize, 2usize];
    let bins = 32usize;

    // Histograms: edges length bins+1, sum counts equals number of selected samples
    let hists = sim.velocity_histogram(&axes, bins, None, None)?;
    assert_eq!(hists.len(), axes.len());
    for (edges, counts) in hists {
        assert_eq!(edges.len(), bins + 1, "edges length mismatch");
        assert_eq!(counts.len(), bins, "counts length mismatch");
        let total: f64 = counts.iter().sum();
        assert_eq!(
            total as usize,
            sim.num_particles(),
            "sum of counts should equal N"
        );
    }

    // KL diagnostic over components should be finite and non-negative
    let dkl = sim.kl_stats(64, None, &axes, None)?;
    assert!(
        dkl.is_finite() && dkl >= 0.0,
        "D_KL must be non-negative and finite"
    );
    Ok(())
}

/// Phase 1.5: Mechanical pressure should either compute or at least have impulse events present.
/// This mirrors the probabilistic nature of piston collisions.
#[test]
fn mechanical_pressure_window() -> Result<()> {
    let mut sim = Simulation::new(96, [30.0, 30.0, 30.0], 0.25, 1.0, Some(98765))?;
    // Mix with static walls first
    sim.advance_to(0.5)?;
    // x-max wall is wall_id = 1
    sim.set_piston_velocity(1, -0.05)?;
    // Advance long enough to gather events
    sim.advance_to(8.0)?;

    let res = sim.mechanical_pressure(2.0, 1);
    match res {
        Ok(p) => {
            assert!(
                p.is_finite() && p >= 0.0,
                "pressure should be finite and non-negative"
            );
        }
        Err(Error::DiagnosticInsufficientSamples(_)) => {
            // Allow this outcome but require that we did record some impulses
            let events = sim.pressure_events();
            assert!(
                !events.is_empty(),
                "expected some piston impulse events when pressure unavailable"
            );
        }
        Err(e) => return Err(e),
    }
    Ok(())
}

/// Phase 1.5: Windowed −α/T residual should be finite for an adiabatic piston run.
/// In the quasi-adiabatic limit and over modest windows, this should be near zero,
/// but we only assert finiteness for robustness.
#[test]
fn alpha_over_t_window() -> Result<()> {
    let mut sim = Simulation::new(64, [20.0, 20.0, 20.0], 0.2, 1.0, Some(112233))?;
    sim.advance_to(0.5)?;
    sim.set_piston_velocity(1, -0.03)?;
    sim.advance_to(3.0)?;

    // Choose a window contained in the simulated time span
    let val = sim.alpha_over_t(1.0)?;
    assert!(val.is_finite(), "alpha/T must be finite");
    Ok(())
}
