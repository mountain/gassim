# 001 semantic and api changes requests

Title: Split Work and Heat at Moving Thermal Walls; Per‑Wall W/Q Accounting & Audited Energy Closure (Phase‑2 Core Update)

---

## Summary

Add event‑level decomposition of particle–wall energy exchange at moving thermal walls so that:

ΔE = u_w·Δp (ΔW: mechanical work) + ΔQ (pure heat)

and expose per‑wall totals and histories for both W and Q. Make “heat” mean pure heat in public APIs by default (mechanical contribution removed), add energy‑audit assertions, and provide by‑wall impulse histories for pressure estimation. This change eliminates double‑counting at moving heat baths, makes the First Law hold per event and per run, and aligns the simulator’s bookkeeping with the contact‑geometric form of thermodynamics used in our theory.

---

## Motivation

1) Physics & correctness. In the current core, heat recorded at a moving thermal wall includes the mechanical work done by that wall; when W is also accumulated, work is counted twice, producing ΔU ≠ W + Q in analyses unless compensated externally. Separating W and Q at the source fixes this.

2) Theory alignment. Our thermodynamics/contact dictionary identifies the thermodynamic contact form α_TD = dU + p dV − T dS with the arithmetic‑expression contact form α = da − μ du − λ a dv via (S,V;U,T,p) = (v, −u; a, λ a, μ). The three channels (energy / work / heat) must be individually accounted to preserve the contactomorphism and to make −α/T‑type δ‑circuit diagnostics meaningful. The same structure (δ calculus, non‑commutativity and δ²a = μλ du∧dv) is developed in the longer manuscript, reinforcing why we must keep the channels disambiguated numerically.

---

## Problem Statement

- A moving wall configured as a thermal bath currently updates heat as
  Q += ½ m (|v_out|² − |v_in|²), which equals heat + mechanical work in the lab frame.
- Analyses then need to subtract W from Q post hoc to satisfy ΔU = W + Q.
- There is no per‑wall breakdown for W/Q, and histories mix concepts.
- We lack core‑level energy audits and unit tests that enforce ΔE = ΔW + ΔQ per event for all wall types.

---

## Proposed Solution

### A) Data model (Rust Simulation measurements)

- work_total: f64
- work_by_wall: [f64; 2*DIM] or Vec<f64>
- heat_total: f64  // pure heat
- heat_by_wall: [f64; 2*DIM]
- impulse_by_wall: [f64; 2*DIM]  // for pressure

Histories (ring buffers; Vec for initial patch is acceptable; ring buffer later):
- work_hist: Vec<(t: f64, dW: f64, wall_id: u32)>
- heat_hist: Vec<(t: f64, dQ: f64, wall_id: u32)>  // pure heat
- pressure_hist: Vec<(t: f64, |Δp_n|: f64, wall_id: u32)>

### B) Wall representation

Make kinematics and boundary condition orthogonal:

```rust
enum Kinematics { Static, Moving { vel: VecD } }
enum Boundary   { Adiabatic, Thermal { T: f64, accommodation: f64 } }

struct Wall {
    id: u32,
    n: VecD,
    pos: f64,
    kin: Kinematics,
    bc: Boundary,
}
```

This allows pure mechanical piston (Moving+Adiabatic), static bath (Static+Thermal), and moving bath (Moving+Thermal).

Note: The current engine maintains axis‑aligned planar walls via two WallKinematics arrays (min/max per axis) and a thermal_walls array. The above conceptual type guides future generalization; for this patch we keep the existing representation and add per‑wall accounting.

### C) Event‑level algorithm (particle–wall)

For any wall:

```rust
// inputs: v_in, normal axis, uw (wall normal speed), bc (Adiabatic/Thermal), m
let v_rel_in  = v_in - uw * n;               // along normal axis in this codebase
let v_rel_out = match bc {
  Adiabatic => specular_reflect(v_rel_in, n),
  Thermal { T, accommodation } => {
      if rng() > accommodation { specular_reflect(v_rel_in, n) }
      else { sample_half_maxwell_outgoing(T, m, n) }  // flux-weighted half-space
  }
};
let v_out = v_rel_out + uw * n;

let dp   = m * (v_out - v_in);                 // particle momentum change (lab)
let dE   = 0.5*m*(|v_out|^2 - |v_in|^2);       // particle energy change (lab)
let dW   = (uw * n) · dp;                       // mechanical work (wall → particle) = uw * Δp_n
let dQ   = dE - dW;                             // pure heat

meas.work_total      += dW;
meas.heat_total      += dQ;
meas.work_by_wall[id] += dW;
meas.heat_by_wall[id] += dQ;

let dp_n_abs = |dp · n|;
meas.impulse_by_wall[id] += dp_n_abs;

work_hist.push((t, dW, id));
heat_hist.push((t, dQ, id));
pressure_hist.push((t, dp_n_abs, id));

debug_assert!((dE - (dW + dQ)).abs() <= 1e-12);
```

Properties:
- Static thermal wall: u_w = 0 ⇒ dW = 0, dQ = dE.
- Moving adiabatic wall: thermal sampling disabled ⇒ dQ = 0, dE = dW.
- Moving thermal wall: both non‑zero but ΔE = ΔW + ΔQ by construction.

### D) Python API (PyO3)

- get_work_heat(*, true_heat: bool = True) -> tuple[float, float]
  - default returns (W_total, Q_total) with Q = pure heat.
  - if true_heat=False, return legacy raw Q (for backward comparison only).
- get_work_by_wall() -> np.ndarray[float64]  // shape (2*DIM,)
- get_heat_by_wall(*, true_heat: bool = True) -> np.ndarray[float64]
- get_heat_history(window: float | None = None, wall_id: int | None = None, *, true_heat: bool = True) -> np.ndarray
  - returns (M,3) rows [t, dQ, wall_id] with heat as pure heat by default.
- get_pressure_history(window: float | None = None, wall_id: int | None = None) -> np.ndarray
  - returns (M,2) rows [t, |impulse|]; if wall_id is None returns concatenated events for all walls, else filter.

Update docs to state unequivocally: “Heat” means pure heat (mechanical contribution removed). The true_heat flag exists solely for migration and reproducibility studies.

### E) Pressure and work estimation (windowed)

- Keep the impulse‑based pressure proxy per wall: P ≈ (Σ |Δp_n|) / (A Δt).
- Keep mechanical work from event sum dW = u_w · Δp; avoid back‑outs like dW/dV.

### F) Energy audits & diagnostics

- Per event (debug builds): assert abs(dE − (dW + dQ)) ≤ eps.
- Expose get_energy_residuals() -> dict with totals since start (optional).
- Metadata: mark whether any wall is moving+thermal for analysis clarity.

---

## Current Implementation Snapshot (baseline before patch)

Files reviewed:
- src/core/sim.rs
- src/lib.rs
- tests/phase2.rs

Observed behavior:
- Mechanical work: work_total accumulates u_w * Δp for piston collisions.
- Heat: heat_total accumulates ½ m (|v⁺|² − |v⁻|²) for all thermal wall interactions, irrespective of wall motion. This means “heat” currently includes the mechanical contribution at moving thermal walls.
- Events:
  - pressure_events: Vec<(t, |impulse|)> recorded only for the active piston wall (piston_wall_id).
  - heat_events: Vec<(t, dQ, wall_id)>, where dQ is the current “raw” value (includes mechanics at moving baths).
- Python API:
  - get_work_heat() -> (W_total, Q_total) without any true_heat flag.
  - get_pressure_history(window=None) without wall_id.
  - get_heat_history(window=None, wall_id=None) exposed, but always returns current core’s “raw” Q.

Conclusion: The proposed change is not yet implemented in code and requires:
- splitting dQ into pure heat by subtracting u_w · Δp at moving thermal walls,
- adding per‑wall accumulators and histories,
- extending APIs with by‑wall accessors and the true_heat switch for migration.

---

## Implementation Plan (file‑by‑file)

This plan reflects minimal, safe edits consistent with AGENT.md (fmt/clippy/tests/docs) and the current axis‑aligned wall model.

1) src/core/sim.rs
- Struct fields:
  - Add work_by_wall: [f64; N_WALLS]
  - Add heat_by_wall: [f64; N_WALLS]
  - Add impulse_by_wall: [f64; N_WALLS]
  - Add work_events: Vec<(f64, f64, u32)>  // optional; can derive from Δp and u_w; useful for debugging
  - Keep heat_events as is (now store pure heat increments)
  - Optionally generalize pressure_events to per‑wall or keep a flat event log with wall_id tagged: Vec<(f64, f64, u32)>
- Initialization:
  - Zero‑init new arrays; preserve existing behavior.
- resolve_p2w():
  - Compute v_before, v_after; compute Δp along axis (normal).
  - Compute dE = ½ m (|v_after|² − |v_before|²).
  - Compute dW = u_w * Δp (signed; “wall → gas”).
  - Compute dQ = dE − dW.
  - Accumulate totals and per‑wall arrays; push events with wall_id; push pressure event with wall_id and |Δp|.
  - Add debug_assert!((dE - (dW + dQ)).abs() <= 1e-12);
- Accessors:
  - Add pub fn work_by_wall(&self) -> Vec<f64> and heat_by_wall(&self) -> Vec<f64>.
  - Change pressure event accessor to optionally return all walls; or add pressure_events_all() returning Vec<(t, |impulse|, wall_id)>.
  - Add work_events() if implemented.
- mechanical_pressure(window, wall_id):
  - Use per‑wall impulse sums; support any wall_id, not just piston.

2) src/lib.rs (PyO3)
- get_work_heat(true_heat: bool = True) -> (f64, f64)
  - If true_heat, return (work_total, heat_total_pure).
  - If false, return (work_total, heat_total_raw), where “raw” sums dE when thermal, i.e., legacy behavior. For migration, you may compute “raw” by adding back u_w·Δp for thermal walls only, or maintain a separate accumulator while the deprecation window is open.
- get_work_by_wall() -> np.ndarray[float64]
- get_heat_by_wall(true_heat: bool = True) -> np.ndarray[float64]
- get_heat_history(window=None, wall_id=None, true_heat: bool = True) -> np.ndarray
- get_pressure_history(window=None, wall_id=None) -> np.ndarray
  - Return rows [t, |impulse|] filtered by wall or concatenated across walls.

3) tests/phase2.rs (and new tests)
- Add unit tests:
  - Static thermal wall: set_thermal_wall on a static wall; advance; assert dW=0 and ΔU ≈ Q (pure).
  - Moving adiabatic wall: set piston velocity; walls adiabatic; assert ΔU ≈ W and Q ≈ 0.
  - Moving thermal wall: set piston velocity and thermal bath on the same wall; assert per‑event |dE − (dW + dQ)| small (via internal debug or by integrating over short steps); over the run, residual ≤ tolerance.
  - Pressure by wall: activate different piston walls; assert per‑wall impulse accumulation indexes correctly.
- Extend integration test:
  - Moving thermal piston regression: with new split, First Law closes without analysis‑side corrections.

4) docs/design.md and docs/supplement.md
- Clarify that “heat” exposed by the API is pure heat by default.
- Note per‑wall totals and histories.
- Cross‑reference δ‑circuit helpers.

5) examples/
- Update isothermal_compression.py (if needed) to reflect the new API defaults and to show per‑wall heat/work printing and pressure history filtering.

6) CI
- Ensure all gates pass:
  - cargo fmt --check
  - cargo clippy -- -D warnings
  - cargo test
  - cargo doc --no-deps

---

## API Changes (Python)

New/updated methods on GasSim:

- get_work_heat(true_heat: bool = True) -> tuple[float, float]
  - Returns (W_total, Q_total). Q_total is pure heat by default.
  - Setting true_heat=False returns legacy raw heat (for one release window).
- get_work_by_wall() -> np.ndarray
  - Returns shape (2*DIM,) per wall id (0..2*DIM−1).
- get_heat_by_wall(true_heat: bool = True) -> np.ndarray
  - Same shape; returns per‑wall Q totals.
- get_heat_history(window: float | None = None, wall_id: int | None = None, true_heat: bool = True) -> np.ndarray
  - Rows [time, dQ, wall_id]; by default dQ is pure heat; filters by window and wall.
- get_pressure_history(window: float | None = None, wall_id: int | None = None) -> np.ndarray
  - Rows [time, |impulse|]; per‑wall filtering supported.

Notes:
- set_piston(wall_id, velocity) and set_thermal_wall(wall_id, T, accommodation) remain unchanged.
- “Heat” in public APIs henceforth means pure heat (mechanical removed).

---

## Migration Guide

- Semantics change: API “heat” now means pure heat. Previously, at moving thermal walls, heat counted mechanical work as well.
- To reproduce legacy analysis during the deprecation window:
  - Use get_work_heat(true_heat=False) and get_heat_by_wall(true_heat=False).
  - Or reconstruct raw heat as Q_raw ≈ Q_pure + Σ_events (u_w · Δp) for thermal walls if raw totals are not persisted.
- Analyses that previously subtracted W from Q should now use Q directly (pure), and W separately.

---

## Acceptance Criteria

Unit tests (Rust):
1) Static thermal wall: many events, assert dW=0, sum(dQ)=ΔU within tolerance.
2) Moving adiabatic wall: assert dQ=0, sum(dW)=ΔU.
3) Moving thermal wall: eventwise |dE − (dW + dQ)| < 1e-12 (double), cumulative closure across run.
4) Impulse/pressure: sign, scaling, and per‑wall accumulation correct.

Integration tests:
1) Isothermal compression baseline (single bath at x‑min, mechanical piston at x‑max, slow):
   - |T_end − T_bath|/T_bath ≤ 2%
   - |ΔU − (W + Q)|/max(U0,U1) ≤ 1e−3
   - W → W_rev and Σ_tot → 0⁺ as |u_pis| ↓
2) Moving thermal piston regression: with new split, First Law closes without analysis‑side corrections.

Performance:
- No measurable regression at typical N (per‑event ops are O(1); histories stored in Vec initially, optional ring buffer later).

---

## Risks & Mitigations

- API semantic change: mitigate with explicit flags, clear docs, and a one‑release deprecation period for “raw” heat.
- Floating‑point audits: use fused ops where available; tolerate ~1e-12 per event audits in double precision.
- Memory growth of histories: plan ring‑buffer compaction or windowed down‑sampling; expose window filters in Python.

---

## Documentation Updates

- design.md §10–§11–§13: redefine heat, add per‑wall interfaces, state event‑level identities and audit.
- Examples: provide two canonical Phase‑2 scripts:
  - Single bath + mechanical piston (mainline)
  - Moving thermal piston (demonstrates split)

---

## Rationale vs. Theory (traceability)

- The 1‑form α_TD = dU + p dV − T dS must pull back to the arithmetic contact form α = da − μ du − λ a dv; keeping dU, p dV (work), and T dS (heat) disentangled in the numerics enforces this correspondence at the level of recorded observables. The δ‑calculus and commutator relations provide the geometric backbone for our δ‑circuit entropy diagnostics; clean W/Q channels are required to compare −α/T window integrals with Clausius‑type statements.

---

## Minimal Patch (guided by current code)

In resolve_p2w (src/core/sim.rs), after computing v_before and deciding specular vs thermal branch:

```rust
// existing
let m = self.particles[i].mass;
let v_after_n = self.particles[i].v[axis];
let delta_p = m * (v_after_n - v_before_n); // signed normal impulse on particle

// NEW: compute pure heat split
let vsq_before: f64 = v_before_vec.iter().map(|&c| c * c).sum();
let vsq_after:  f64 = self.particles[i].v.iter().map(|&c| c * c).sum();
let dE = 0.5 * m * (vsq_after - vsq_before);
let dW = u_w * delta_p;           // wall normal speed * impulse (lab)
let dQ = dE - dW;                 // pure heat
self.work_total                 += dW;
self.heat_total                 += dQ;
self.work_by_wall[wall_id as usize] += dW;
self.heat_by_wall[wall_id as usize] += dQ;
self.impulse_by_wall[wall_id as usize] += delta_p.abs();
self.heat_events.push((self.time_now, dQ, wall_id));
// optionally: self.work_events.push((self.time_now, dW, wall_id));
// and record pressure event with wall_id tag if moving

debug_assert!((dE - (dW + dQ)).abs() <= 1e-12);
```

Adjust pressure event storage to include wall_id, or maintain a per‑wall vector.

---

## Release Plan

- Version: v0.0.5 (Phase‑2 Core Update)
- Scope:
  - Implement event‑level split and per‑wall accounting.
  - Add audits & tests.
  - Ship updated Python API + two example scripts.
- Deprecation:
  - Legacy “raw heat” available behind true_heat=False for one minor release window.
  - Document migration in CHANGELOG and README.

---

## CHANGELOG Snippet

feat(core): split work/heat at moving thermal walls; per‑wall W/Q/impulse totals and histories; event‑level First‑Law audits

- Heat semantics now default to pure heat (mechanical contribution removed).
- New per‑wall accumulators and histories for W, Q, and |Δp|.
- New Python API: get_work_by_wall(), get_heat_by_wall(), get_heat_history(..., true_heat), get_pressure_history(..., wall_id).
- Add per‑event debug assertions enforcing ΔE = ΔW + ΔQ.
- Migration helpers: get_work_heat(true_heat=False) to retrieve legacy “raw heat”.

---

## Open Questions

- Time‑varying accommodation per wall?
- Cellwise field sampling (ρ, u, T) in core vs. Python (initially keep in Python).
- Ring buffer implementation for histories (capacity and compaction policy).

---

## Developer Checklist (for this change)

- [ ] Update src/core/sim.rs fields for per‑wall and histories.
- [ ] Implement dE/dW/dQ split and audits in resolve_p2w().
- [ ] Add accessors for per‑wall work/heat/impulse and pressure events by wall.
- [ ] Extend mechanical_pressure() to all walls.
- [ ] Update src/lib.rs functions and signatures to include new methods and flags.
- [ ] Update/extend tests: unit (per‑event closure) and integration (moving thermal piston).
- [ ] Update docs: design.md, examples, README snippets; add CHANGELOG.
- [ ] Ensure all gates pass: fmt, clippy (no warnings), tests, docs.
