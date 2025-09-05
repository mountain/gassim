use gassim::core::Simulation;

/// Phase 0 validation: Energy conservation in an NVE hard-sphere gas.
/// Advance for many collisions and assert total kinetic energy remains constant
/// within a tight floating-point tolerance.
#[test]
fn energy_conservation_nve() -> gassim::error::Result<()> {
    let mut sim = Simulation::new(64, [20.0, 20.0, 20.0], 0.2, 1.0, Some(12345))?;
    let e0 = sim.kinetic_energy();

    let t_end = 50.0;
    let steps = 200;
    for s in 1..=steps {
        let t = t_end * (s as f64) / (steps as f64);
        sim.advance_to(t)?;
    }

    let e1 = sim.kinetic_energy();
    let rel = ((e1 - e0) / e0).abs();
    assert!(
        rel < 1e-8,
        "relative energy drift {} too large (E0={}, E1={})",
        rel,
        e0,
        e1
    );
    Ok(())
}

/// Phase 0 validation: Approximate isotropy at equilibrium.
/// After evolving for a while, the per-axis mean squared velocities should be comparable.
/// We require max-min spread < 20% of the mean as a loose statistical check.
#[test]
fn isotropy_approx_after_collisions() -> gassim::error::Result<()> {
    let mut sim = Simulation::new(128, [30.0, 30.0, 30.0], 0.2, 1.0, Some(7777))?;
    sim.advance_to(100.0)?;

    let mut sum_sq = [0.0_f64; 3];
    let n = sim.num_particles() as f64;
    for p in &sim.particles {
        for (k, vk) in p.v.iter().enumerate() {
            sum_sq[k] += vk * vk;
        }
    }
    for x in &mut sum_sq {
        *x /= n;
    }
    let mean = (sum_sq[0] + sum_sq[1] + sum_sq[2]) / 3.0;
    let maxv = sum_sq.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let minv = sum_sq.iter().cloned().fold(f64::INFINITY, f64::min);
    let spread = (maxv - minv) / mean;
    assert!(
        spread < 0.2,
        "anisotropy too high: <vx^2>={}, <vy^2>={}, <vz^2>={}, spread={}",
        sum_sq[0],
        sum_sq[1],
        sum_sq[2],
        spread
    );
    Ok(())
}
