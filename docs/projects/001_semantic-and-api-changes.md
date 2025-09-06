
# 001 semantic and api changes requests

**Title**: Split Work and Heat at Moving Thermal Walls; Per‑Wall W/Q Accounting & Audited Energy Closure (Phase‑2 Core Update)

---

## Summary

Add **event‑level decomposition** of particle–wall energy exchange at **moving thermal walls** so that:

$$
\Delta E \;=\; \underbrace{u_w\!\cdot\!\Delta p}_{\Delta W\ \text{(mechanical work)}} \;+\; \underbrace{\Delta Q}_{\text{(pure heat)}}
$$

and expose **per‑wall totals and histories** for both $W$ and $Q$. Make **“heat” mean *pure heat*** in public APIs by default (mechanical contribution removed), add **energy‑audit assertions**, and provide **by‑wall impulse histories** for pressure estimation. This change eliminates double‑counting at moving heat baths, makes the First Law hold **per event and per run**, and aligns the simulator’s bookkeeping with the contact‑geometric form of thermodynamics used in our theory.

---

## Motivation

1. **Physics & correctness.** In the current core, heat recorded at a **moving thermal wall** includes the mechanical work done by that wall; when `W` is also accumulated, work is counted twice, producing $\Delta U \neq W + Q$ in analyses unless compensated externally. Separating $W$ and $Q$ at the source fixes this.

2. **Theory alignment.** Our thermodynamics/contact dictionary identifies the thermodynamic contact form

$$
\alpha_{\rm TD}=dU+p\,dV-T\,dS
$$

with the arithmetic‑expression contact form

$$
\alpha=da-\mu\,du-\lambda a\,dv
$$

via $(S,V;U,T,p)=(v,-u; a,\lambda a,\mu)$. The three channels (energy / work / heat) must be **individually** accounted to preserve the contactomorphism and to make $-\alpha/T$‑type δ‑circuit diagnostics meaningful.   The same structure (δ calculus, non‑commutativity and $\delta^2 a=\mu\lambda\,du\wedge dv$) is developed in the longer manuscript, reinforcing why we must keep the channels disambiguated numerically.&#x20;

---

## Problem Statement

* A **moving** wall configured as a **thermal bath** currently updates heat as

  $$
  Q \mathrel{+}= \tfrac12 m(|v_{\rm out}|^2-|v_{\rm in}|^2),
  $$

  which equals **heat + mechanical work** in lab frame.
* Analyses then need to subtract $W$ from $Q$ post hoc to satisfy $\Delta U = W + Q$.
* There is no **per‑wall** breakdown for $W$/$Q$, and histories mix concepts.
* We lack **core‑level energy audits** and **unit tests** that enforce $\Delta E=\Delta W+\Delta Q$ per event for all wall types.

---

## Proposed Solution

### A) Data model (Rust `Simulation.meas`)

* `work_total: f64`
* `work_by_wall: Vec<f64>`
* `heat_total: f64`  **// pure heat**
* `heat_by_wall: Vec<f64>` **// pure heat**
* `impulse_by_wall: Vec<f64>` **// for pressure**
* Histories (ring buffers):

    * `work_hist: RingBuf<(t: f64, dW: f64, wall_id: u32)>`
    * `heat_hist: RingBuf<(t: f64, dQ: f64, wall_id: u32)>` **// pure heat**
    * `pressure_hist: RingBuf<(t: f64, |Δp_n|: f64, wall_id: u32)>`

### B) Wall representation

Make kinematics and boundary condition **orthogonal**:

```rust
enum Kinematics { Static, Moving { vel: VecD } }
enum Boundary { Adiabatic, Thermal { T: f64, accommodation: f64 } }

struct Wall {
    id: u32,
    n: VecD,
    pos: f64,
    kin: Kinematics,
    bc: Boundary,
}
```

This allows *pure mechanical piston* (Moving+Adiabatic), *static bath* (Static+Thermal), and *moving bath* (Moving+Thermal).

### C) Event‑level algorithm (particle–wall)

For any wall:

```rust
// inputs: v_in, wall.n, uw (from Kinematics), bc (Adiabatic/Thermal), m
let v_rel_in  = v_in - uw;
let v_rel_out = match bc {
  Adiabatic => specular_reflect(v_rel_in, n),
  Thermal { T, accommodation } => {
      if rng() > accommodation { specular_reflect(v_rel_in, n) }
      else { sample_half_maxwell_outgoing(T, m, n) }  // half-space, flux weighted normal
  }
};
let v_out = v_rel_out + uw;

let dp   = m * (v_out - v_in);                 // particle momentum change (lab)
let dE   = 0.5*m*(v_out.norm2() - v_in.norm2());// particle energy change (lab)
let dW   = uw.dot(dp);                          // **mechanical work (wall → particle)**
let dQ   = dE - dW;                             // **pure heat**

meas.work_total += dW;
meas.heat_total += dQ;
meas.work_by_wall[id] += dW;
meas.heat_by_wall[id] += dQ;

let dp_n_abs = dp.dot(n).abs();
meas.impulse_by_wall[id] += dp_n_abs;

work_hist.push((t, dW, id));
heat_hist.push((t, dQ, id));
pressure_hist.push((t, dp_n_abs, id));
```

Properties:

* **Static thermal wall:** $u_w=0 \Rightarrow dW=0,\ dQ=dE$.
* **Moving adiabatic wall:** Thermal sampling disabled ⇒ $dQ=0,\ dE=dW$.
* **Moving thermal wall:** Both non‑zero but **$\Delta E=\Delta W+\Delta Q$** by construction.

### D) Python API (PyO3)

* `get_work_heat(*, true_heat: bool = True) -> tuple[float, float]`

    * default returns **(W\_total, Q\_total)** with $Q$ = *pure heat*.
    * if `true_heat=False`, return legacy *raw* Q (for backward comparison only).
* `get_work_by_wall() -> np.ndarray[float64]`
* `get_heat_by_wall(*, true_heat: bool = True) -> np.ndarray[float64]`
* `get_heat_history(window: float | None = None, wall_id: int | None = None, *, true_heat: bool = True) -> np.ndarray`
* `get_pressure_history(window: float | None = None, wall_id: int | None = None) -> np.ndarray`

Update docs to state unequivocally: **“Heat” means pure heat (mechanical contribution removed).**

### E) Pressure and work estimation (windowed)

* Keep the **impulse‑based** pressure proxy per wall:
  $P \approx \big(\sum |\Delta p_n|\big) / (A\,\Delta t)$.
* Keep **mechanical work** from event sum $dW=u_w\!\cdot\!\Delta p$; avoid theoretical back‑outs like $dW/dV$.

### F) Energy audits & diagnostics

* Per event (debug builds): assert `abs(dE - (dW + dQ)) <= eps`.
* Expose `get_energy_residuals() -> dict` with totals since start.
* Metadata: mark whether any wall is moving+thermal for analysis clarity.

---

## Acceptance Criteria

**Unit tests (Rust):**

1. **Static thermal wall:** many events, assert $dW=0$, `sum(dQ)=ΔU` within tolerance.
2. **Moving adiabatic wall:** assert $dQ=0$, `sum(dW)=ΔU`.
3. **Moving thermal wall:** eventwise `|dE-(dW+dQ)|<1e-12` (double), cumulative closure.
4. **Impulse pressure:** sign, scaling, and by‑wall accumulation correct.

**Integration tests:**

1. **Isothermal compression baseline** (single bath at x‑min, **mechanical** piston at x‑max, slow):

    * $|T_{\rm end}-T_{\rm bath}|/T_{\rm bath}\le 2\%$
    * $|\Delta U-(W+Q)|/\max(U_0,U_1)\le 10^{-3}$ (double)
    * $W\to W_{\rm rev}$ and $\Sigma_{\rm tot}\to 0^+$ as $|u_{\rm pis}|\downarrow$.
2. **Moving thermal piston** (for regression): with new split, First Law closes **without** analysis‑side corrections.

**Performance:**

* No measurable regression at typical N (histories are ring buffers; work/heat math is $O(1)$ per event).

---

## Backward Compatibility

* Default heat semantics change to **pure heat**.

    * Add release note and migration guide.
    * Provide `true_heat=False` switches to retrieve **legacy** raw $Q_{\rm core}$ for reproducibility studies.

---

## Risks & Mitigations

* **API semantic change:** mitigate with explicit flags, clear docs, and a one‑release deprecation period for “raw” heat.
* **Floating‑point audits:** use fused ops where available; tolerate `~1e-12` per event audits in double precision.

---

## Documentation Updates

* *Design Doc §10–§11–§13:* redefine heat, add by‑wall interfaces, state event‑level identities and audit.
* *Examples:* provide two canonical Phase‑2 scripts:
  (i) **Single bath + mechanical piston** (mainline);
  (ii) **Moving thermal piston** (demonstrates split).

---

## Open Questions

* Do we want to expose **accommodation** per wall as time‑varying profiles?
* Add optional **cellwise field sampling** (ρ, u, T) in core, or keep it in Python first?

---

## Rationale vs. Theory (traceability)

* The 1‑form **$\alpha_{\rm TD}=dU+p\,dV-T\,dS$** must pull back to the arithmetic contact form **$\alpha=da-\mu\,du-\lambda a\,dv$**; keeping $dU$, $p\,dV$ (work), and $T\,dS$ (heat) **disentangled** in the numerics enforces this correspondence at the level of recorded observables. **Proposition 3.1** establishes this contactomorphism precisely.&#x20;
* The δ‑calculus and commutator relations ($[D_u,D_v]=\mu\lambda\,\partial_a$, $\delta^2 a=\mu\lambda\,du\wedge dv$) provide the geometric backbone for our δ‑circuit entropy diagnostics; clean W/Q channels are required to compare $-\alpha/T$ window integrals with Clausius‑type statements.&#x20;

---

## Minimal Patch (pseudocode)

```rust
// in wall collision resolver
let v_rel_in  = v_in - uw;
let v_rel_out = match bc { /* adiabatic or sampled thermal outflux */ };
let v_out     = v_rel_out + uw;

let dp   = m * (v_out - v_in);
let dE   = 0.5*m*(v_out.norm_sqr() - v_in.norm_sqr());
let dW   = uw.dot(dp);
let dQ   = dE - dW;

meas.work_total += dW;
meas.heat_total += dQ;
meas.work_by_wall[id] += dW;
meas.heat_by_wall[id] += dQ;
meas.impulse_by_wall[id] += dp.dot(n).abs();

work_hist.push((t, dW, id));
heat_hist.push((t, dQ, id));
pressure_hist.push((t, dp.dot(n).abs(), id));

debug_assert!((dE - (dW + dQ)).abs() <= 1e-12);
```

---

## Requested Milestone

* **v0.0.5** (Phase‑2 Core Update):

    * Implement event‑level split and per‑wall accounting.
    * Add audits & tests.
    * Ship updated Python API + two example scripts.

This feature will make GasSim’s mechanical/thermal bookkeeping physically correct by construction, remove analysis‑side hacks, and bring the engine’s observables into exact correspondence with the contact‑geometric framework we use for δ‑circuit and entropy analysis.
