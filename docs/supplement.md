# GasSim Design Document — **Supplement**

> This section augments the original document. Existing sections are preserved; insert the following **⟨ADD⟩** items under the indicated headers.

---

## 2. High‑Level Architecture

**⟨ADD⟩** *Thermal contact and local thermometry*

* Core engine supports two thermal boundary primitives:

    1. **Thermal Wall** (Maxwell accommodation at a box face; piston may simultaneously be moving & thermal).
    2. **ThermoPlane** (an infinitesimally thin internal plane for *in‑bulk contact thermometry*; transmissive Maxwell “patch” with probability α).
* Observability plane/channel for **heat events** $dQ$ parallels the work/impulse channel; enables *contact‑temperature* readout via *zero‑heat‑flux condition* and closes the **δ‑circuit** with $dU = \delta Q - p\,dV$. The contact 1‑form $dU + p\,dV - T\,dS$ maps to AEG’s $\alpha = da - \mu\,du - \lambda a\,dv$ under the dictionary $(S,V;U,T,p)=(v,-u; a, \lambda a, \mu)$, guaranteeing structural consistency.&#x20;

---

## 5. Core Data Structures

**⟨ADD⟩** new/extended types

* **ThermalBcKind**

  ```rust
  enum ThermalBcKind {
      // Classical wall (boundary of domain)
      Wall { T: f64, accommodation: f64 },             // accommodation ∈ [0,1], 1=Maxwell, 0=specular
      // Internal transmissive patch for contact thermometry
      Plane { T: f64, accommodation: f64, normal: VecD, point: VecD, area: f64 },
  }
  ```
* **Boundary**

  ```rust
  struct Boundary {
      id: u32,
      kind: ThermalBcKind,
      // geometry: plane x·n = pos for faces; for Plane, (normal, point, finite area or infinite flag)
      n: VecD,
      pos: f64,      // for axis-aligned faces
  }
  ```
* **EventKind** **⟨MOD⟩**
  Add:

    * `P2Thermal { i: u32, bnd_id: u32 }` (covers both Wall & Plane)
* **Measurements** **⟨ADD⟩**

  ```rust
  struct HeatEvent { t: f64, dQ: f64, bnd_id: u32 }
  struct PressureEvent { t: f64, impulse_abs: f64, wall_id: u32 }
  struct WindowCfg { dt_window: f64, max_len: usize }
  struct Moments { sum_c: VecD, sum_cc: [[f64; D]; D], count: usize } // for temperature tensor
  ```

  Accumulators:

    * `heat_total: f64; heat_history: RingBuffer<HeatEvent>`
    * `work_total: f64; pressure_history: RingBuffer<PressureEvent>`
    * optional sliding-window **velocity moments** for anisotropy metrics and KL diagnostics (see §10).

---

## 6. Event Prediction and Resolution

### 6.4 Thermal wall (heat bath) — **detailed spec** **⟨ADD⟩**

* **Collision rule (Maxwell accommodation α):** with probability `α` resample in wall frame; with `1-α` do specular reflection.

    * Let wall normal $\mathbf n$ (into gas), wall velocity $\mathbf U$ (piston allowed). Pre‑collision: $\mathbf c^-=\mathbf v^- - \mathbf U$.
    * **Resampling** at bath temperature $T_b$ (m=mass, $k_B=1$ in code units):

        * Tangential: $c^+_{t1},c^+_{t2}\sim \mathcal N(0,\sigma^2), \ \sigma^2=T_b/m$.
        * Normal (into gas): $c^+_n \sim \text{Rayleigh}(\sigma)$ (flux‑weighted half‑space Maxwell).
    * Transform back: $\mathbf v^+ = \mathbf c^+ + \mathbf U$.
* **Heat accounting:** $\Delta Q = \tfrac{m}{2}\big(|\mathbf v^+|^2 - |\mathbf v^-|^2\big)$ (positive when wall → gas). Push `(t, ΔQ, bnd_id)` into `heat_history`; accumulate `heat_total`.
* **Notes:** In 2D, one tangential component; in 3D, two. The Rayleigh sampler may use inverse‑transform with $r=\sigma\sqrt{-2\ln(1-u)}$.

### 6.5 ThermoPlane (internal contact thermometer) **⟨ADD⟩**

* Geometric predicate: free‑flight line–plane intersection (like P2W), but plane is **transmissive**.
* On crossing, with probability `α` perform **transmissive Maxwellization** (same sampler as thermal wall but keep the sign of the normal so that particle continues across the plane); otherwise pass through unchanged.
* Log each modified crossing as a **heat event** with the same $\Delta Q$ definition.
* This enables **in‑bulk contact thermometry** by tuning $T_b$ so that the window‑averaged $\dot Q \to 0$ on that plane; see Python control loop in §11.

---

## 7. Event Loop and Time Advancement

**⟨ADD⟩** measurement hooks

* After every resolved event:

    * If `P2W` on a piston wall: accumulate impulse & mechanical work.
    * If `P2Thermal`: accumulate `ΔQ`.
    * Update sliding‑window **velocity moments** for anisotropy/KL if enabled (O(1) amortized per event with Welford‑style updates).
* Provide a low‑overhead **window compactor** (down‑sample) to keep histories bounded.

---

## 8. Numerical Robustness and Determinism

**⟨ADD⟩** thermal randomness & reproducibility

* Thermal sampling uses a dedicated RNG stream per boundary id; seeds recorded in metadata for replay.
* Tie‑break precedence: `P2P` > `P2W` > `P2Thermal` on equal time within ε, then by (i,j,wall\_id) index for determinism.

---

## 9. Performance Design

**⟨ADD⟩** selective diagnostics

* Temperature‑tensor & KL diagnostics are **opt‑in** and can be computed on a **spatial mask** (bulk vs near‑wall), or on a coarse slab grid to cap cost.
* For large N, export periodic **moments** instead of full velocities; Python builds tensors & histograms.

---

## 10. Measurements and Observables

**⟨ADD⟩ 10.1 Temperature definitions (operational)**

* **Drift‑removed kinetic temperature (scalar):**
  $\mathbf u = \langle \mathbf v\rangle,\ \mathbf c_i=\mathbf v_i-\mathbf u,\  T = \tfrac{m}{D N}\sum_i \lVert \mathbf c_i\rVert^2.$
* **Temperature tensor:**
  $T_{jk}=\tfrac{m}{N}\sum_i c_{i,j}c_{i,k},\quad T=\tfrac{1}{D}\mathrm{tr}\,T.$
* **Anisotropy indicators:**
  $\eta_\infty=\max_j|T_{jj}/T-1|,\ \ \eta_F=\|T - T I\|_F /(T\sqrt{D}).$

**⟨ADD⟩ 10.2 Pressure and work (clarification)**

* **Mechanical pressure at piston:** use **work–volume slope** $p_{\rm mech} = -\Delta W/\Delta V$ or windowed **impulse/area/time** (already present). For EOS checks, compare bulk $p_{\rm th}=\tfrac{1}{D}\text{tr}\,P\approx n k_BT$ vs wall‑normal measure; differences quantify near‑wall anisotropy.

**⟨ADD⟩ 10.3 Heat and contact temperature**

* **Heat flux per boundary:**
  $\dot Q_b \approx \sum_{k\in \text{window}} \Delta Q_k / \Delta t.$
* **Contact thermometry (operational):** a Python controller adjusts boundary $T_b$ until $\dot Q_b\to 0$ → returns **contact temperature** $T_{\rm contact}$ for that site/plane.

**⟨ADD⟩ 10.4 Distribution diagnostics (KL/KS)**

* 1D (per component) KL: histogram $P$ of $c_x$ vs $Q=\mathcal N(0,\sigma^2)$ with $\sigma^2=T/m$:
  $D_{\rm KL}(P\|Q)=\sum_j p_j \ln(p_j/q_j)$.
* 3D speed version: compare radial histogram to Maxwell speed PDF.
* Report time series of KL (bulk and near‑wall).

**⟨ADD⟩ 10.5 δ‑Circuit measurement (non‑equilibrium entropy accounting)**

* **Baseline** (ideal monatomic gas, $k_B=1$):
  $dS_{\rm IG} = \tfrac{1}{T}dU + \tfrac{p_{\rm th}}{T}dV,\quad \Delta S_{\rm IG}= N\!\left[\tfrac{3}{2}\ln\frac{T_2}{T_1}+\ln\frac{V_2}{V_1}\right].$
* **δ‑circuit 1‑form:**
  $\alpha_{\rm IG} = dU - \underbrace{(-p_{\rm mech}\,dV + T\,dS_{\rm IG})}_{\delta a_{\rm IG}}.$
  The **vertical residual** $-\alpha_{\rm IG}/T$ is the running *entropy‑production measure*; integrate over time/windows.
* **AEG consistency:** under the contactomorphism $(S,V;U,T,p)=(v,-u; a,\lambda a,\mu)$, $\alpha_{\rm TD}=dU+p\,dV-T\,dS$ pulls back to $\alpha = da - \mu\,du - \lambda a\,dv$ (Eq. (33)–(35)), so the δ‑circuit is geometrically faithful. &#x20;

---

## 11. Python API (PyO3)

**⟨ADD⟩ Phase 1.5 diagnostics**

* `get_temperature_tensor(mask: Optional[np.ndarray]) -> (Txx,Tyy,Tzz,T_scalar,count)`
* `get_velocity_histogram(axes=("x","y","z"), bins=80, range=None, mask=None) -> dict`
* `get_kl_stats(mask=None) -> float` (returns current window KL for chosen mode)

**⟨ADD⟩ Thermal contact / thermometry (Phase 2)**

* `set_thermal_wall(wall_id: int, T: float, accommodation: float = 1.0) -> None`
* `set_thermo_plane(plane_id: int, point: Sequence[float], normal: Sequence[float], T: float, accommodation: float = 0.02, area: Optional[float]=None) -> None`
* `get_heat_history(boundary_id: int, window: float) -> np.ndarray  # [t, dQ]`
* **Controller pattern (Python‑side):** `find_contact_temperature(boundary_id, tol=..., max_iter=...)` adjusts `T` to drive `⟨\dot Q⟩ → 0`.

**⟨ADD⟩ δ‑circuit helpers**

* `get_mechanical_pressure(window: float, wall_id: int) -> float`
* `get_entropy_IG(time_series: bool = False) -> (S_vals or ΔS)`
* `get_alpha_over_T(window: float) -> float  # windowed −α/T`
* `get_work_heat() -> (W_total, Q_total)`

Implementation notes: all data are **windowed**; masks allow bulk vs near‑wall selection.

---

## 12. Error Handling and Logging

**⟨ADD⟩** diagnostics guardrails

* KL/anisotropy modules validate binning and sample counts; return `DiagnosticInsufficientSamples` if below threshold.

---

## 13. Testing and Validation

**⟨ADD⟩ Unit/Integration (new)**

1. **Thermal sampler tests**: moments of sampled $\mathbf c^+$ match Maxwell half‑space within statistical error; Rayleigh law for $c_n$.
2. **First Law closure** (Phase 2): $\Delta U = \sum_b Q_b + W$ over runs (tolerance relative to $U$).
3. **Contact temperature**: single thermal wall at $T_b$ relaxes to $T_b$; thermo‑plane controller recovers $T_{\rm contact}$ with $|T_{\rm contact}-T_{\rm kin}|$ bulk median ≤ 5%.
4. **δ‑circuit “must‑pass”** (Phase 1/2):

    * Same‑end‑volume compression ensemble: all runs share $V_{\rm end}$.
    * Quasi‑static limit: $\langle -\alpha/T\rangle$ increases monotonically with |piston speed|, →0 as speed→0.
    * Work–energy residual small (`|ΔU - W|` in Phase 1 adiabatic; `|ΔU - W - Q|` in Phase 2).
    * Baseline correlation: deviations in $T$ & wall pressure vs adiabatic/EOS baselines positively correlate with $\langle -\alpha/T\rangle$.
      These checks mirror the contact‑geometric **horizontal/vertical** split in AEG (§2.1–§2.3).&#x20;

**⟨ADD⟩ Diagnostic thresholds (defaults)**

* Bulk anisotropy: $\eta_F < 0.1$ typical; warn if exceeded persistently.
* KL (bulk): $D_{\rm KL}<10^{-2}$ typical for near‑equilibrium; warn on sustained spikes.

---

## 14. Reproducibility and Randomness

**⟨ADD⟩** record per‑boundary RNG seeds and controller states (`find_contact_temperature` iterations) in metadata for exact reruns.

---

## 16. Milestones and Alignment with docs/plan.md

**⟨ADD⟩ Phase 1.5 (Diagnostics & δ‑circuit)**

* Deliver anisotropy/KL diagnostics, δ‑circuit accounting, and “must‑pass” CI gate.
* Scripts: `examples/phase1_adiabatic_alpha.py`, `examples/diagnostics_qc.py`.

**⟨ADD⟩ Phase 2 (Thermal contact & thermometry)**

* Implement `Thermal Wall` and `ThermoPlane`; expose heat histories; First‑Law regression.
* Example: piston compression with one thermal wall; conduction test between `T_L ≠ T_R`.

---

## 18. Future Enhancements

**⟨ADD⟩ Local fields and geometry**

* Slab/cell binning to export $\rho(\mathbf r), \mathbf u(\mathbf r), T(\mathbf r)$ and **contact‑temperature field** from sparse thermo‑planes (kriging/IDW blend with kinetic‑temperature map).
* Variable accommodation, moving thermo‑planes (Lagrangian probes), and specular/roughness models.

---

## 19. Appendix (new): Maxwell sampling cookbook **⟨ADD⟩**

* **Tangential components:** Box–Muller or Ziggurat for $\mathcal N(0,\sigma^2)$.
* **Positive normal (flux‑weighted):** sample $R\sim \text{Rayleigh}(\sigma)$, set $c_n^+=R$.
* **Specular mixing:** with prob. $1-\alpha$, reflect $c_n^+ = -c_n^-$, keep tangentials.
* **Heat increment:** $\Delta Q = \tfrac{m}{2}(\|\mathbf v^+\|^2-\|\mathbf v^-\|^2)$.
* **Contact controller (pseudo‑code):**

  ```python
  # every Δt window
  qdot = sum(dQ)/Δt
  if abs(qdot) < tol: done
  T_b -= kappa * qdot  # simple integral control; or bracket root by bisection
  ```

---

### Notes on Theory Alignment (for the Methods section of the paper)

* The simulator’s δ‑circuit is the numerical instantiation of the contact 1‑form $\alpha_{\rm TD}=dU+p\,dV-T\,dS$; under the AEG dictionary $(S,V;U,T,p)=(v,-u;a,\lambda a,\mu)$, it pulls back to $\alpha=da-\mu\,du-\lambda a\,dv$. The **horizontal** part $\delta a=\mu\,du+\lambda a\,dv$ records reversible bookkeeping; the **vertical** residual $\alpha$ gives the non‑commutative “area/entropy‑production” measure $-\alpha/T$. This correspondence justifies using macroscopic $U,V,T,p$ alone to assess irreversibility, independent of microstate enumerations. &#x20;
