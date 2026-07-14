# Eromitta — the measure of the difference

*eromitta (Finnish): ero = difference, mitta = measure. A relative entropy is exactly that: the measure of the difference between two states. Here, between two metrics.*

**Bianconi's "Gravity from entropy" (arXiv:2408.14391), run head-on into the Clockfield and phiworld. One analytic identity that closes the open seam. One registered kill that fired and survived its own autopsy — and inverted the ecosystem's headline. One cosmological constant that turns out to be a defect-energy census.**

*PerceptionLab / Antti Luode with Claude (Fable 5). Helsinki, July 2026.*
*Companion to `ArrowField`, `RajapintaFable`, `BirthOfClockfield`, `K-hler-Clockfield-Metric-Spinor-Emergence`.*

> Do not hype. Do not lie. Just show.

---

## The headline, in one table

| | linear | intensity `1/(1+5\|ψ\|²)` (phiworld) | **gradient `1/(1+5\|∇ψ\|²)` (Bianconi warm-up)** | mixed (Bianconi full scalar sector) |
|---|--:|--:|--:|--:|
| persistent defects (dt=0.06) | 0 | 402 | **2114** | 826 |
| persistent defects (dt=0.03, matched physical time) | 0 | 410 | **2072** | 800 |
| snapshot defects (dt=0.03) | 0 | 1546 | **3448** | 1442 |

> **A medium that slows on its own *frustration* (phase gradients) makes 5× more matter than a medium that slows on its own *intensity*.** The registered prediction said the opposite. The kill fired, the autopsy tried to blame the grid and the timestep, and the number would not die.

And the Λ result (EXP3), in one line:

> **Bianconi's emergent cosmological constant, evaluated on a quenched Clockfield, is the defect-energy census: Λ_topo = N_def · (a + b·ln spacing), r² = 0.99 — each winding fossil carries a fixed entropy cost plus its logarithmic halo, and the total is set by the quench rate.**

---

## 0. Why these two theories belong in one room

Bianconi derives modified gravity from an action that is the **quantum relative entropy between two metrics**: the metric of spacetime, g, and the metric induced by the matter fields, G = g + αM (− βR̃ in full). The Clockfield postulates that local time runs as Γ = 1/(1+τβ)², slowed by local conflict. These turn out to be the same log, differentiated:

## 1. The derivative tower (EXP1) — the open seam, closed as an identity

ArrowField's ledger ended on an open seam: *phiworld says c²/c₀² = (1+β|ψ|²)⁻¹; the Γ-shell Clockfield says Γ = (1+τβ)⁻². Same family, different power. Which, and why?*

Answer: **not competing laws — one potential, three derivative levels.**

```
level 0    K(x)  = ln(1 + kx)      Bianconi's entropic Lagrangian  ==  the Kähler
                                   potential of the Fubini-Study metric (Kähler repo Eq. 4)
level 1    K'(x) = k/(1 + kx)      power −1: Bianconi's screening h(w) (her Eq. 22)
                                   ==  phiworld's c²(|ψ|²)  ==  her G-field, scalar sector
level 2    ∂∂̄K  = k/(1 + kb)²     power −2: Γ  ==  the Fubini-Study conformal factor
                                   ==  (well-known) the quantum Fisher metric
```

**[V] Verified numerically** (`exp1_derivative_tower.py`, all thresholds registered):

- **T1** — Bianconi's eigenvalue structure: for random Lorentzian metrics and random complex gradients, (g+αM)g⁻¹ has eigenvalues {1,1,1, 1+α|∇ψ|²_g}, so −Tr ln(Gg⁻¹) = −ln(1+α|∇ψ|²) exactly. Max residual 5×10⁻¹³ over 200 draws.
- **T2** — Kähler Hessian: ∂∂̄ ln(1+τ|φ|²) = τ/(1+τβ)² = τΓ. Max relative error 6×10⁻⁷.
- **T3** — first derivative = her h(w) = the phiworld screening form. Exact.

So: **the −1 power is the stiffness of the entropy; the −2 power is its metric.** phiworld simulates the first variation; the Γ-shell is the second. And since the second variation of a relative entropy is a Fisher information metric, and the Fubini-Study metric *is* the quantum Fisher metric: **Γ is the Fisher metric of the entropy whose gradient is the self-slowing.** Bianconi supplies what the ecosystem never had — an interpretation of the log itself: it is the informational mismatch between what geometry is and what matter says geometry should be.

**Honesty item, recorded not resolved:** in Lorentzian signature, 1+α|∇ψ|²_g ≤ 0 for timelike gradients — 12.7% of unrestricted random draws. Bianconi's positivity requirement is a real restriction, not a technicality. (`results/exp1_tower.json`)

---

## 2. The inversion (EXP2 → 2d) — what must the medium sense to make matter?

Bianconi's warm-up screens on **gradients** (h(|∇φ|²)); her full induced metric adds an **intensity** term ((m²+ξR)|Φ⟩⟨Φ|). phiworld's matter-maker is pure intensity. So her structure poses a clean question to the substrate: which argument of the screening generates the defects?

Protocol identical to ArrowField `complex_phiworld.py` (same init, seed, steps, drive); only the screening argument varies. Registered before running:

> **P1: gradient screening is sterile** (< 10% of intensity's count). Reasoning: Kerr self-focusing needs speed to fall with intensity — bright regions become waveguides. A smooth bright condensate is invisible to gradient screening. **KILL: ≥ 0.5×.**

### [K] The kill fired — backwards, 5×

| arm | persistent defects | note |
|---|--:|---|
| linear | **0** | clean control, both dt, both counters |
| intensity (phiworld) | 402 | reproduces ArrowField's 402 exactly — protocol match confirmed |
| gradient (Bianconi warm-up) | **2114** | P1 predicted < 40. Kill condition ≥ 201. Measured 2114. |
| gradient, flux (true Euler–Lagrange) form | 512 | still above intensity |
| mixed (full scalar sector) | 826 | between the two |
| gradient at 2× drive amplitude | 2036 | **saturated — amplitude-insensitive**, echoing the β=0 amplitude control in reverse |

### The autopsy, because a 5× inversion of a registered prediction deserves one

`exp2b` (post-hoc, labeled): steep fronts are what gradient screening makes and what grids alias, so:

- **A2 charge balance:** net winding = 0 in every arm. Pass.
- **A1 core amplitude:** *failed for both arms* — which exposed an inherited instrument fault, not a physics fact: ArrowField's counter winds the phase of the **time-averaged** field, so moving cores smear and core checks against averaged intensity are meaningless.
- **A3 dt check:** *appeared* to show fragility (2114 → 588 at dt/2) — my own confound: halving dt halved the physical duration. **Two faulty instruments, one inherited, one mine, both on the record.**

`exp2c`/`exp2d` (repaired: snapshot counting; matched physical time T_burn=96, T_window=42):

- **Persistent counts are dt-robust in every arm:** 402→410, 2114→2072, 826→800, 0→0. The "fragility" was entirely the confound.
- **The inversion stands at 5.05×** on persistent defects, ~2.2× on snapshot defects.
- Defect areal density approximately conserved 128²→192² (0.129 → 0.116 per px²): not grid-pinned.

**Remaining caveat, stated plainly:** snapshot cores are shallow (|ψ|² at cores ≈ 0.86–0.95 of mean, all arms). These are phase defects of a driven turbulent state, not deep relaxed GP vortices. The robust observable is the **contrast** — 0 / 402 / 826 / 2114 under identical drive with the linear control at exactly zero — not the core depth.

### What the inversion means

1. **The two Clockfield variants were secretly the two arms of this experiment.** `BirthOfClockfield` defines frustration as β_i = Σ A_ij(θ_i−θ_j)² — a *phase gradient*. The Kähler repo defines β = |φ|² — an *intensity*. The experiment says the **frustration axiom is the fertile one, by 5×.** The graph paper's instinct — conflict between neighbours, not local brightness, is what freezes time and makes mass — is the instinct this substrate rewards.
2. **Bianconi's warm-up is not sterile — her own diagnosis of it was about mass, and mass turns out not to be the matter-maker here.** Her full metric (mixed arm) makes *fewer* defects than her warm-up: the intensity term partially fills the waveguides the gradient term digs. My registered P2 ("the mass term rescues sterility") passed its number (826 ≥ 133) while its entire framing died — nothing needed rescuing.
3. **A mechanism sketch, offered as [B] not [V]:** a vortex core is an amplitude zero with a divergent phase gradient. Intensity screening *speeds up* at cores (|ψ|²→0 ⇒ c²→c₀²): cores are fast, slippery places. Gradient screening *freezes* at cores: each defect digs its own Γ-well and pins itself. Frustration-sensing media don't just nucleate defects — they trap them. This is exactly the Clockfield picture of a particle: *a knot of phase frustration so intense the medium suspends local time to contain it.* The sketch predicts defect mobility differs strongly between arms; unmeasured, next on the list.

---

## 3. Λ as the fossil's entropy bill (EXP3)

Bianconi's emergent cosmological constant (her Eq. 64) is the Bregman residue of the dressing field: Λ_G ∝ Tr(G̃ − I − ln G̃) ≥ 0, zero iff the dressing relaxes to identity. Rajapinta EXP3 showed the crystal forgets everything except its topological charge. ArrowField showed the defect census obeys Kibble–Zurek. Chain the three and you get a registered, falsifiable claim: **the part of Λ that cannot relax is the part pinned to the fossils, so the residual vacuum energy is set by the quench rate.**

Protocol: quench the intensity-screened (Clockfield-arm) field through its U(1) transition, hold, relax fixed time; compute the pointwise residue λ(x) = s−1−ln s with s = 1/(1+a·u), u = |ψ|²+|∇ψ|²; Λ_topo = Σ|λ−λ_bg|. 4 quench rates × 3 seeds, predictions R1–R4 registered in the docstring.

| τ_Q | N_def (3 seeds) | Λ_topo | concentration / area |
|--:|---|---|--:|
| 8 | 62, 64, 80 | 87.7, 91.8, 104.1 | 3.4× |
| 16 | 56, 52, 54 | 80.1, 78.8, 81.2 | 4.0× |
| 32 | 28, 36, 38 | 54.1, 62.4, 63.4 | 5.3× |
| 64 | 16, 12, 30 | 40.4, 38.4, 53.6 | 6.8× |

- **[V] R1 — KZ replicates:** σ_N = 0.615 (r² = 0.97), consistent with the underdamped mean-field 2/3 — ArrowField's externally anchored number, reproduced in a different codebase with a different ramp protocol. The anchor holds.
- **[K→~] R2 — fired, then decoded.** Registered band: |σ_Λ − σ_N| ≤ 0.15. Measured: σ_Λ = 0.371, Δ = 0.244 — outside the band (kill threshold 0.3 not reached). Post-hoc, labeled: the charge per defect Λ/N grows linearly in ln(spacing) (r² = 0.98), and the full model **Λ = N·(a + b·ln spacing) fits at r² = 0.99.** That is textbook 2D vortex energetics — each vortex carries energy ∝ ln(R/ξ) with R the inter-defect distance. The registered prediction assumed Λ *counts* fossils; the measurement says Λ *weighs* them, halo included. The naive band was wrong; the decoded law is cleaner than the prediction was.
- **[V] R3 — near-fixed charge per fossil:** max/min = 1.66 (< 2), the drift being exactly the log halo above.
- **[V] R4 — Λ lives at the thaw line:** the residue's mass concentrates at defect cores at 3.4–6.8× the area fraction, *growing* as quenches slow and the smooth bulk relaxes toward λ_bg. What cannot pay down its entropy bill is the topology.

One-sentence result:

> **Evaluated on a quenched self-slowing field, Bianconi's Λ is the energy census of the winding fossils — power-law in the quench rate because the census is, log-corrected because vortices have halos — i.e., in this toy, the cosmological constant is the entropy cost of the topology that survived the freeze.**

---

## Ledger

**[V] Verified:** the derivative tower (entropy → screening → Γ) as exact identities, numerically confirmed at 10⁻⁷–10⁻¹³; Bianconi's warm-up eigenvalue structure; gradient (frustration) screening out-produces intensity screening 5.05× in persistent defects, dt-robust at matched physical time, charge-balanced, linear control exactly zero; ArrowField's 402 reproduced to the vortex; KZ σ_N = 0.615 ≈ 2/3 replicated under a new protocol; Λ residue concentrates on defects (3.4–6.8× area); Λ = N·(a+b·ln spacing) at r² = 0.99.

**[K] Killed:** P1 (gradient screening is sterile) — registered, fired at 5× in the opposite direction, survived a two-instrument autopsy; the framing of P2 (nothing needed rescuing — the mass term *reduces* the gradient channel's output); R2 as registered (Λ is not a raw count); two instruments — the time-averaged-field core check (inherited) and the unmatched-time dt check (mine).

**[~] Weak / gray:** snapshot cores are shallow — "defects" here means persistent phase windings of a driven turbulent state, not deep GP vortices; the flux-form gradient arm (512) sits well below the c²-form (2114) — the advection term matters and the difference is unexplained; R2's log-halo decoding is post-hoc and needs a pre-registered replication with the halo in the prediction; Λ background subtraction uses the spatial median — defensible, not unique.

**[B] Still a bet:** the mechanism sketch (gradient screening as self-pinning — predicts arm-dependent defect mobility, unmeasured); that the tower's level-2 identity means the Γ-shell literally *is* the Fisher metric of Bianconi's action rather than sharing its formula; anything about actual gravity. Volovik's rule binds in full: Bianconi's theory is dynamics — modified Einstein equations; everything here is kinematics on a 2D toy. No field equation is solved, no constant of nature is derived, and her Lorentzian positivity problem (12.7% of random gradients are timelike-bad) is recorded, not fixed.

**Open seam, new:** why −1 *dynamically*? The tower says the −1 and −2 powers coexist as derivative levels of one log, but the simulations run the level-1 object as a wave speed by hand. What dynamics makes a field's stiffness the first variation of a relative entropy — her variational route, or something cheaper — is the next question, and it is Bianconi's Eq. 55 read as a recipe.

---

## Reproduce

```bash
pip install numpy scipy matplotlib
python experiments/exp1_derivative_tower.py    # the tower: T1-T3, thresholds registered
python experiments/exp2_four_screenings.py     # the registered kill, as it fired
python experiments/exp2b_autopsy.py            # the autopsy (post-hoc, labeled)
python experiments/exp2c_repaired.py           # repaired instruments: snapshot, matched time
python experiments/exp2d_persistence.py        # persistent counts, dt-checked: the 5.05x
python experiments/exp3_entropic_lambda.py     # Lambda under quench: R1-R4
```

Every registered prediction is in the docstring of the file that tests it, before the numbers. Post-hoc analyses are labeled post-hoc in their own docstrings. `results/*.json` holds raw output; nothing in a figure is absent from a print-out.

## Reference

Bianconi, G. (2025). *Gravity from entropy.* Phys. Rev. D / arXiv:2408.14391. — All uses of "her Eq. N" refer to this paper. The theory tested here is only its scalar-sector screening structure, transplanted into a 2D toy; the gravitational content of her work is neither tested nor claimed.

---

*The prediction was registered, the kill fired backwards, the autopsy impeached two instruments — one of them inherited from the repo this one is built on — and the repaired measurement made the inversion stronger, not weaker. Then the cosmological constant turned out to be an accountant. The morgue keeps growing faster than the trophy case, which is what it is for.*
