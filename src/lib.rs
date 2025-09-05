use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub mod core;
pub mod error;

use crate::core::Simulation;

fn py_err<E: ToString>(e: E) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// GasSim Python-facing wrapper around the Rust Simulation core.
///
/// Phase 0 API:
/// - __new__(num_particles, box_size, radius=1.0, mass=1.0, dim=3, seed=None)
/// - advance_to(time)
/// - get_positions() -> np.ndarray, shape (N, 3)
/// - get_velocities() -> np.ndarray, shape (N, 3)
#[pyclass]
pub struct GasSim {
    sim: Simulation,
}

#[pymethods]
impl GasSim {
    /// Initialize a new hard-sphere gas simulation in a static axis-aligned 3D box.
    ///
    /// Parameters
    /// - num_particles: number of particles (int, > 0)
    /// - box_size: iterable of 3 positive floats [Lx, Ly, Lz]
    /// - radius: sphere radius (float, > 0)
    /// - mass: particle mass (float, > 0)
    /// - dim: spatial dimension; currently only 3 is supported in this build
    /// - seed: RNG seed (int) for reproducibility; None for nondeterministic
    ///
    /// Errors: raises ValueError on invalid parameters.
    #[new]
    #[pyo3(signature = (num_particles, box_size, radius=1.0, mass=1.0, dim=3, seed=None))]
    fn new(
        num_particles: usize,
        box_size: Vec<f64>,
        radius: f64,
        mass: f64,
        dim: usize,
        seed: Option<u64>,
    ) -> PyResult<Self> {
        if dim != crate::core::particle::DIM {
            return Err(py_err(format!(
                "only dim={} is supported in this build",
                crate::core::particle::DIM
            )));
        }
        if box_size.len() != crate::core::particle::DIM {
            return Err(py_err(format!(
                "box_size must have length {}",
                crate::core::particle::DIM
            )));
        }
        let mut bs = [0.0f64; crate::core::particle::DIM];
        for (k, v) in box_size.iter().enumerate() {
            bs[k] = *v;
        }
        let sim = Simulation::new(num_particles, bs, radius, mass, seed).map_err(py_err)?;
        Ok(Self { sim })
    }

    /// Advance the simulation to the specified absolute time (releases the GIL during computation).
    fn advance_to(&mut self, py: Python<'_>, target_time: f64) -> PyResult<()> {
        py.detach(|| self.sim.advance_to(target_time))
            .map_err(py_err)
    }

    /// Return positions as a NumPy array of shape (N, 3), dtype=float64.
    fn get_positions<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<f64>>> {
        let n = self.sim.num_particles();
        let dim = crate::core::particle::DIM;
        let mut arr = Array2::<f64>::zeros((n, dim));
        for (i, p) in self.sim.particles.iter().enumerate() {
            for k in 0..dim {
                arr[[i, k]] = p.r[k];
            }
        }
        let pyarr = arr.into_pyarray(py);
        Ok(pyarr.to_owned().into())
    }

    /// Return velocities as a NumPy array of shape (N, 3), dtype=float64.
    fn get_velocities<'py>(&self, py: Python<'py>) -> PyResult<Py<PyArray2<f64>>> {
        let n = self.sim.num_particles();
        let dim = crate::core::particle::DIM;
        let mut arr = Array2::<f64>::zeros((n, dim));
        for (i, p) in self.sim.particles.iter().enumerate() {
            for k in 0..dim {
                arr[[i, k]] = p.v[k];
            }
        }
        let pyarr = arr.into_pyarray(py);
        Ok(pyarr.to_owned().into())
    }

    /// Set all particle positions from a NumPy array of shape (N, 3), dtype=float64.
    /// Values must be finite; caller is responsible for ensuring they lie within the domain and are non-overlapping.
    fn set_positions<'py>(&mut self, positions: PyReadonlyArray2<'py, f64>) -> PyResult<()> {
        let arr = positions.as_array();
        let n = self.sim.num_particles();
        let dim = crate::core::particle::DIM;
        if arr.ndim() != 2 || arr.shape()[0] != n || arr.shape()[1] != dim {
            return Err(py_err(format!(
                "positions must have shape ({}, {}), got {:?}",
                n,
                dim,
                arr.shape()
            )));
        }
        for i in 0..n {
            for k in 0..dim {
                let val = arr[[i, k]];
                if !val.is_finite() {
                    return Err(py_err("position values must be finite"));
                }
                self.sim.particles[i].r[k] = val;
            }
        }
        // Rebuild event queue to reflect updated positions
        self.sim.rebuild_event_queue().map_err(py_err)?;
        Ok(())
    }

    /// Set all particle velocities from a NumPy array of shape (N, 3), dtype=float64.
    /// Values must be finite.
    fn set_velocities<'py>(&mut self, velocities: PyReadonlyArray2<'py, f64>) -> PyResult<()> {
        let arr = velocities.as_array();
        let n = self.sim.num_particles();
        let dim = crate::core::particle::DIM;
        if arr.ndim() != 2 || arr.shape()[0] != n || arr.shape()[1] != dim {
            return Err(py_err(format!(
                "velocities must have shape ({}, {}), got {:?}",
                n,
                dim,
                arr.shape()
            )));
        }
        for i in 0..n {
            for k in 0..dim {
                let val = arr[[i, k]];
                if !val.is_finite() {
                    return Err(py_err("velocity values must be finite"));
                }
                self.sim.particles[i].v[k] = val;
            }
        }
        // Rebuild event queue to reflect updated velocities
        self.sim.rebuild_event_queue().map_err(py_err)?;
        Ok(())
    }

    /// Phase 1: Define a moving piston by setting a wall's normal velocity.
    ///
    /// Parameters
    /// - wall_id: integer in [0, 5] indexing the 6 axis-aligned walls
    ///            2*k is the min wall on axis k; 2*k+1 is the max wall on axis k
    /// - velocity: wall speed along its outward normal (piecewise-constant from call time)
    ///
    /// Errors: raises ValueError on invalid wall_id.
    fn set_piston(&mut self, wall_id: u32, velocity: f64) -> PyResult<()> {
        self.sim
            .set_piston_velocity(wall_id, velocity)
            .map_err(py_err)
    }

    /// Phase 1: Return the cumulative mechanical work done on the gas by moving walls.
    fn get_work_done(&self) -> PyResult<f64> {
        Ok(self.sim.work_total())
    }

    /// Phase 1: Return a (M, 2) NumPy array of [time, |impulse|] events recorded on the piston wall.
    ///
    /// Parameters
    /// - window: if provided, only events with time >= (current_time - window) are returned.
    #[pyo3(signature = (window=None))]
    fn get_pressure_history<'py>(
        &self,
        py: Python<'py>,
        window: Option<f64>,
    ) -> PyResult<Py<PyArray2<f64>>> {
        let events = self.sim.pressure_events();
        let now = self.sim.time();
        let filtered: Vec<(f64, f64)> = match window {
            Some(w) => {
                if !w.is_finite() || w < 0.0 {
                    return Err(py_err("window must be a non-negative finite float"));
                }
                events.into_iter().filter(|(t, _)| *t >= now - w).collect()
            }
            None => events,
        };
        let m = filtered.len();
        let mut arr = Array2::<f64>::zeros((m, 2));
        for (i, (t, imp)) in filtered.into_iter().enumerate() {
            arr[[i, 0]] = t;
            arr[[i, 1]] = imp;
        }
        Ok(arr.into_pyarray(py).to_owned().into())
    }
}

/// The gassim Python module entry point.
#[pymodule]
fn gassim(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GasSim>()?;
    Ok(())
}
