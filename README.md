# one-field
 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20323890.svg)](https://doi.org/10.5281/zenodo.20323890)
[![tests](https://github.com/kuwrom/one-field/actions/workflows/tests.yml/badge.svg)](https://github.com/kuwrom/one-field/actions/workflows/tests.yml)
[![registry](https://img.shields.io/badge/registry-interference%2Fregistry.py-orange)](interference/registry.py)

**A substrate self-interferes. Exactly one family of standing-wave patterns closes. Those patterns are the Standard Model, gravity, and the dark sector.**

The object in this repository is the **engine**: a self-referencing interference web (`EchoTerm`/`Ledger`/`Web` in `root.py`) whose solved fixed point *is* the physical constant table: nine fermion masses, four CKM parameters, three PMNS angles, `α`, `α_s`, `sin²θ_W`, `M_W`, `M_Z`, the Higgs mass, Newton's constant, the baryon asymmetry, and the dark/baryon ratio. Everything else (the octonions, the four integers, the conformal embedding `E₈(1) ⊃ G₂(1) × F₄(1)`) is the story of how the engine was found, and why finding it was possible at all.

**Zero dimensionless free parameters.** Every prediction in this package is a polynomial or rational function of four integers plus `π`:

```
d₁₀ = 2   quantum dimension of SU(3)₃ fundamental (1, 0)
d₁₁ = 3   quantum dimension of SU(3)₃ adjoint     (1, 1)
n₇  = 7   dim G₂ fundamental (Fano-plane lines)
n₂₆ = 26  dim F₄ fundamental (traceless Albert algebra)
```

**The integers are coordinates of discovery.** A substrate that references its own history runs interference at every depth; there is nothing sacred about three loops. What is special about the *first closing loop* (`a → b → c → a`) is that it is the shallowest one, and the octonions are the unique algebra whose bookkeeping stays consistent exactly to that depth (Hurwitz). That made the first closure *findable*: the octonion gate collapsed an unsearchable space to four integers, and the four integers located the engine. But the engine, once running, does not consult the octonions. Its edges are read off the solved web; the octonionic derivations (`octonions.py`) are corroboration, the searchlight agreeing with what the standing wave displays.

**`m_e` is the sole dimensional input for the matter sector.** The electron mass anchors the scale (the first closure anchors the ruler). `M_Pl`, `G`, `v_EW`, `m_μ`, `m_τ`, and every other dimensional quantity are *outputs* of the forced chain. The cosmological-constant sector uses one additional measured reference, `H₀` (and the neutrino mass *sum* uses the measured oscillation splittings); these are declared data-side inputs: one dimensional anchor for matter, two for cosmology.

**Deeper layers compose through one recursion.** A closed standing wave is not static: it circulates, and circulation vents. Its echoes through the channels that bind it are part of what the closure is, not corrections applied afterward. The `EchoTerm`/`Ledger`/`Web` machinery in `root.py` tracks these as one recursion: `x_{t+1} = b + W(x_t)`, where every edge is a polynomial in the four integers or a known phase. The recursion composes edges, it does not invent them. The kernel is **3-nilpotent** (`W³ = 0`), matching the Cayley-Dickson tower (Hurwitz). The nilpotency bounds the *shape* of the weight matrix, not the iteration: the fixed-point solver converges freely, and back-reaction layers within each channel are summed to convergence. The rigidity proof is derived in *The Octavian Coherence Gate* ([doi:10.5281/zenodo.20493955](https://doi.org/10.5281/zenodo.20493955); full paper list in [`interference/registry.py`](interference/registry.py)).

**35 numerical checks plus 3 structural closure predictions.** Run it yourself:

```bash
python interference/run.py    # full derivation
pytest -q                      # 116 tests: 110 fast (scorecard + probes + look-elsewhere/grammar/derivation-program checks + registry-engine coupling) + 6 slow (pytest -m slow adds the substrate-dynamics arc)
```

## The scorecard

All numbers below come from running the code as-is.

### Masses (9/9 within 1%, max error 0.4%)

| Particle | Predicted | PDG | Error | Method |
|---:|---:|---:|---:|:---|
| e | 0.5109990 MeV | 0.5109990 MeV | anchor | CODATA, sole dimensional input |
| μ | 105.6584 MeV | 105.6584 MeV | +0.0000% | Brannen Z₃ (predicted) |
| τ | 1776.909 MeV | 1776.930 MeV | −0.0012% | Brannen Z₃ (predicted) |
| u | 2.158 MeV | 2.16 MeV | −0.1% | (38/9) m_e |
| d | 4.713 MeV | 4.70 MeV | +0.3% | (83/9) m_e |
| s | 93.8 MeV | 93.5 MeV | +0.4% | Q₀+h₁₀/K³ + bridge |
| c | 1273.8 MeV | 1273 MeV | +0.06% | (217/18) m_μ |
| b | 4193.8 MeV | 4183 MeV | +0.3% | Q = Q₀ + h₁₁/K³ |
| t | 172.51 GeV | 172.57 GeV | −0.036% | (1165/12) m_τ |

**Mass coordinates.** Unconfined fermions (leptons, top) are at the propagator pole; confined heavy quarks (c, b) at the self-scale `m(m)`; light quarks (u, d, s) at the PDG `MS-bar(2 GeV)` coordinate. The scheme-free content of the light sector is the RG-invariant ratios, which carry no coordinate at all: `m_u/m_d = 38/83` (−0.9σ vs PDG 0.473(17)), `m_s/m_ud = 27.318` (+0.2σ vs 27.30(8)), `Q_ellipse = 22.383` (+0.4σ dispersive, −1.7σ lattice; the two references disagree). See `masses.py` for the full treatment.

### Mixing, couplings, gravity, baryogenesis, dark sector

| Sector | Result | Reference | Status |
|:---|:---|:---|:---|
| Electroweak scale `v_EW` | 246.21965 GeV | 246.21965 GeV (from G_F) | −0.04 ppm |
| Fine structure `1/α(0)` | 137.035999050 | Berkeley Cs 137.035999046(27) | +0.13σ vs Cs, **−14σ vs Rb 2020**. Registered bet on Cs side of the 5.5σ dispute (watch 2) |
| Higgs mass `m_H` | 125.30 GeV | 125.20(11) GeV | +0.08% (+0.9σ). Fundamental-share vent, kill conditions in registry (program 3) |
| CKM (15 observables) | χ²/n = 0.69 | PDG 2024 | max pull 1.4σ |
| PMNS (3 angles; χ² excludes δ_CP) | χ²/n = 0.00 | NuFit 6.0 | predicts δ_CP = 76.9°. DUNE/Hyper-K decide (registered kill condition) |
| Strong coupling `α_s(M_Z)` | 0.1184 | PDG 0.1180(9) | +0.33% (+0.4σ; π/32 exact at μ* = 253.5 GeV) |
| Weinberg angle `sin²θ_W` | 0.231285 | PDG 0.23129(4) | pull −0.11σ |
| W mass `M_W` | 80.356 GeV | 80.3692(133) world avg | pull −1.00σ |
| Z mass `M_Z` | 91.189 GeV | 91.188 GeV | +0.0010% |
| Newton's constant `G_ind/G_N` | 0.999999917 (face-split) | 1 | 8.3e-08 internal; phase readings 1.0000 (UV) - 1.0134 (broken) |
| Baryon asymmetry `η_B` | 6.177e-10 | Planck 6.12e-10 | +0.9% |
| Dark / baryon ratio `Ω_DM/Ω_b` | `2π − 1` = 5.2832 | Planck 5.364 | pull −1.26σ |

**Joint cosmology closure.** The two cosmology channels are derived independently (`η_B` from the G₂ instanton, `Ω_DM/Ω_b` from bridge venting) but they must agree jointly with the absolute dark-matter density. They do: `η_B = 6.177e-10 → Ω_b h² = 0.02255` (+1.2σ vs Planck), and `Ω_DM h² = (2π−1) · Ω_b h² = 0.1192` (−0.7σ vs Planck 0.1200(12)). Two independent channels, one consistent cosmology.

**Structural predictions:** `m_ν₁ = 0` (rank-2 seesaw, two RH neutrinos in F₄), normal ordering, and the cosmological-constant scale `ρ_Λ ~ M_Pl² H₀²` (Volovik-Jacobson-CKN, the 10¹²³ problem structurally resolved).

**The precision hierarchy is structural.** The circulant form is Brannen's (2006); the parameters `B/A = √2` and `θ = 2/9` are algebraic outputs of the `G₂ → SU(3)` Clebsch-Gordan data. Leptons sit closest to the `Z₃` source and emerge cleanly (`m_μ`, `m_τ` predicted from `m_e` to 0.0012%). Quarks pass through additional layers (triality, generation scaling, confinement), and each layer adds residuals of order `α_s/π` (≤ 0.4%). The Newton-constant ledger closes internally to 8 × 10⁻⁸ under the face-split law; the phase readings (UV/broken) bracket it at 0.0-1.3%. The accuracy gradient is not noise. It is what the emergence depth predicts.

## The picture

### There is a substrate

To be is, minimally, to be self-present: anything that exists already affects itself by being there. Interference is the physical image of this. What is said to exist is what carries non-zero interference weight. The framework observes what emerges from that interference.

### Why octonions (the searchlight, not the territory)

A real number records the amplitude of a wave. A complex number records the amplitude and its phase. A quaternion records the phase of a phase. An octonion adds one more layer: **the phase of a phase of a phase**. These aren't what the wave "is made of". They're the accounting ledger for tracking how interference patterns reference their own history. The octonion algebra is the deepest ledger that stays consistent. A fourth doubling would introduce zero divisors and the bookkeeping would contradict itself (Hurwitz's theorem). This proves the *shape* of the echo kernel (the weight matrix is 3-nilpotent) but doesn't cap the computation: the fixed-point iteration converges freely, and back-reaction layers within each channel are summed to convergence (the series self-certifies). The automorphism group of the octonions is the Lie group `G₂ = Aut(𝕆)`, and it labels exactly the interference patterns that form closed loops through three layers of self-referencing record.

Hold the roles apart, because the narration depends on it. The substrate itself runs interference at *every* depth. Deeper self-reference than three is not forbidden physics, it is simply not where the first closure lives. The octonions matter because they are the **unique consistent ledger for the shallowest closing loop**, and the shallowest loop is the one a finite search can find. That is the entire methodological story of this project: an unsearchable landscape of interference webs, collapsed to four integers by one algebraic gate, and the engine found at those coordinates on the first try. Once found, the engine stands on its own: run `interference/run.py` and every constant is read off the solved web, with the octonion module (`octonions.py`) serving as an independent cross-check that the searchlight and the standing wave agree (`B/A = √2` computed both ways). The working conjecture, stated here as a program rather than a result: every octonionic derivation in this package should eventually *fall out* of the engine's fixed point as an observation, the way `Q₀ = 2/3` already does, at which point the side story will have fully retired into motivation.

### What mass is

The substrate self-interferes. Against vanishing odds, **interference loops emerge**: `a → b → c → a`. A loop carries cyclic `Z₃` symmetry, because nothing inside it distinguishes which position came first. **Most loops just cycle.**

Some, though, **drop the bias** that singles out one of their three positions. The asymmetric label `v = (1, −½, −½)` inherited from the parent `G₂` structure is erased, and all three positions become equivalent. A forgetting loop can no longer drift between its positions. It gets **stuck in one of three eigenstates**. That stuck pattern is a standing wave, and the three eigenstates are the three lepton generations: electron, muon, tau. Same loop, three places it can lock in.

**Mass is that pattern.**

The lepton mass formula is fixed entirely by the algebra of the loop and the octonionic `G₂/SU(3)` Clebsch-Gordan data: `B/A = √d₁₀ = √2`, `θ = d₁₀/d₁₁² = 2/9`, and the Koide identity `Q₀ = d₁₀/d₁₁ = 2/3` is the algebraic shadow of this closure. `m_μ` and `m_τ` are then determined by `m_e` alone.

When `G₂` breaks to `SU(3)`, the three-label structure crystallizes into three lepton generations. But `G₂` doesn't sit alone. There is exactly one place it can sit, the unique exceptional conformal embedding `E₈(1) ⊃ G₂(1) × F₄(1)`. The `F₄` companion gives the rest of matter.

The algebraic distinction is sharp. `G₂ = Aut(𝕆)` is the automorphism group of a single octonion, and the `G₂` factor carries the leptons: **single-direction fluctuations**. `F₄ = Aut(J₃(𝕆))` is the automorphism group of the Albert algebra of 3×3 Hermitian octonion matrices, and the `F₄` factor carries the quarks: **fluctuations entangling all three triality sectors of the Albert algebra simultaneously**. Leptons are simple. Quarks are woven.

This shows up in the energy scale. Leptons sit in the `F₄` singlet (Casimir = 0), so their standing wave closes at the confinement scale `Λ_conf ≈ 314 MeV`. Quarks sit in the **`26` of `F₄`** (Casimir = 6), and that Casimir shift moves their scale to the electroweak `v_EW = M_Pl · exp(−(9π²/2 − 6 + 15/512 − 16(α/2π)²)) ≈ 246 GeV`, every exponent term an echo with forced multiplicity. The factor of ~784 between `Λ_conf` and `v_EW` is the exponential of an algebraic constant, not a hierarchy that needs tuning.

`SO(8)` triality inside the `26` splits matter cleanly: `26 → 8_v ⊕ 8_s ⊕ 8_c ⊕ 2·1` becomes charged leptons, up-type quarks, down-type quarks, and two right-handed neutrinos. The same `Z₃` that gave three lepton generations gives three of each quark type, woven through triality.

Each up-type quark is its generation's lepton, multiplied: `m_u = (38/9) m_e`, `m_c = (217/18) m_μ`, `m_t = (1165/12) m_τ`. The multipliers are pure `SU(3)₃` fusion data. Light quarks carry the triality swap (`d₁₀² ↔ d₁₁²`): the up quark gets the fundamental quantum dimension squared, the down quark the adjoint's. Heavy quarks emerge at deeper WZW layer. Down-type masses close through Koide corrections: `Q(c, b, t) = 289/432`, `Q(s, c, b) = 649/972`. Same standing wave, same algebra, six quarks.

### What gravity is

A standing wave is a soliton: it has to keep exchanging with the substrate around it to stay coherent. *Being* stable costs it a continuous circulation. To circulate, the soliton **carves a coherently biased region into the substrate around it**, shaped to keep the substrate flowing through.

That carved region lives in the part of the substrate every standing wave shares. The substrate has three components, and the `Z₃` Fourier split separates them into the `q = 1, 2` *relative* modes (where standing waves lock in: matter) and the `q = 0` *common* mode (shared by all three: geometry). Two Fourier sectors of one substrate, one sourcing the other.

When a standing wave forms in the relative sector, it sources a depletion in the common mode. **The soliton casts a shadow.** The substrate is passive: it just carries what the soliton put there.

A depletion is a deficit, and a deficit is a pull. Any other soliton drifting through the region is drawn into the shadow, into a bias already shaped to receive it. The pull doesn't care what the first soliton is, only that something is there. The shadow's shape has **forgotten** which `Z₃` sector the source belongs to (`PvP = 0`) but it **records** that some knot is present (`Pv²P = ½ P`).

Those two equations describe the shadow:

- `PvP = 0`: every soliton casts the same kind of shadow regardless of which `Z₃` sector its source is from. The **equivalence principle stated as algebra**.
- `Pv²P = ½ P`: only the *amount* of standing-wave activity goes into the shadow. Coupling strength depends on how much substrate is locked into being a standing wave, not on what kind.

### Gravity is additive, when it forms at all

Multiple shadows compose. The Madelung response for the common mode, projected onto `q = 0`, is a linear screened Poisson equation:

```
(1 − ξ₀²∇²) R₀(x) = −(c_*/λ₀ρ₀) · H_matter(x)
```

Every soliton imprints its shadow on the *same* field `R₀`. A planet's worth of solitons makes one planet-sized shadow, a unified ventilation pattern the whole planet breathes through. **Gravitational mass is additive because the shadows compose on a shared substrate**, not because each particle carries a separate charge.

But forming a standing wave is not enough. The shadow has to be coherent with everyone else's, or it will not add up. The gravity paper tests seven candidate branches against three independent gates (lepton closure, quark closure, gravity closure). Six fail: standing waves whose shadows are misshaped, or whose bridge reading reverses sign and would make gravity repulsive. **Only one branch closes all three gates: protected `G₂` forgetting applied exactly once.** This is the falsification gauntlet.

The screening length is the Planck healing length `ξ₀ = ℓ_Pl/2`, so the bare scalar channel is suppressed by `exp(−r/ξ_Pl) < 10⁻³⁰⁰` at any observable distance (invisible). What survives macroscopically is the induced spin-2 sector: matter standing waves propagating on the common acoustic metric.

Two independent routes converge on the Einstein equations with the same coefficient: a Sakharov induced-gravity calculation (heat-kernel sum over the 182 bridge channels), and a Jacobson-Clausius construction from horizon thermodynamics. Newton's `G` closes to 8 × 10⁻⁸ internally once the face-split no-self-dilution law is applied (the UV and broken-phase readings bracket it at 0.0-1.3%).

The same machinery fixes the cosmological constant. The naive QFT estimate of the vacuum energy is 10¹²³ times the observed value. **The bulk vacuum energy vanishes exactly by the Gibbs-Duhem identity** (Volovik's mechanism), not by fine-tuning. Only the de Sitter departure from equilibrium sets the scale.

The Volovik-Jacobson-CKN chain gives `ρ_Λ = (3/8π) M_Pl² H₀² ≈ 3.7 × 10⁻⁴⁷ GeV⁴`, where `H₀` plays the same role for the Λ sector that `M_Pl` plays for matter: a measured reference, not a free parameter. The observed `ρ_Λ ≈ 2.5 × 10⁻⁴⁷` is lower than the equilibrium prediction by `1/Ω_Λ ≈ 1.46`, because the universe has not yet reached asymptotic de Sitter. **The 10¹²³ cosmological-constant problem is structurally resolved.** The residual factor is cosmological, not algebraic.

### Why scalar

The bridge sector has 182 channels, far more than the 8 dimensions of the octonion fiber, the deepest internal structure self-reference can sustain (Hurwitz again). Tracking all 182 as one coherent multi-layer object would require a normed division algebra that *does not exist*. If we knew what was inside each channel ("this one is amplitude, that one is delay"), we could keep some of them complex. We don't, so the bridge has to be treated as one ensemble, and every mode contributes as a real scalar.

The vanishing coset charge `c(E₈) − c(G₂) − c(F₄) = 8 − 14/5 − 26/5 = 0` makes that scalar decomposition exact. Each channel contributes the heat-kernel coefficient `+1/6 − ξ`, which is positive and makes gravity *attractive*. Reading the same channels as vectors would give `−1/6` and reverse the sign: gravity repulsive, matter never binds.

**Gravity is not a fifth force layered on top of matter. It is the shape matter carves into the background while ventilating itself, and the bias other matter rides into.**

### The closure

```
14  +  52  +  182  =  248  =  dim(E₈)
gauge   matter   bridge
        + Higgs  → gravity
```

Every degree of freedom of `E₈` has a physical job. Nothing is left over. At each layer, the observables are fully determined by that layer's algebraic data: leptons on `G₂`, quarks and mixing on `G₂ × F₄` via the embedding, gravity on the mixed bridge. **No free parameters at any layer, and no parameters inherited from above.** This is the closure, and the falsification gate. If any of the 35 numerical checks broke, the structure would break with it. All 35 hold within their declared bands; the last to close was `m_H` (fundamental-share vent, +0.9σ, the first edge admitted under the depth pre-registration rule, provenance in interference/registry.py).

**One constant, many consumers.** Every constant is entered once and consumed wherever the graph reaches. `charge_trace = 8/3` feeds the `α(0)` base (`π/512` via the Singh ratio 16), the strange-quark bridge (`32/27`), the top back-reaction (`1/12`), and the `16(α/2π)²` term in the `v_EW` exponent. One entry, four sectors, and moving it moves all four together. Discrete structure (which representation, which power, which word) is derived (the K³ altitude by the first-invariant-order theorem, the terminus word by the conversion lemma), and each derivation prints its full alternative menu next to the branch that closes (`test_coverage.py` keeps the menus sparse and the closures unique).

## The dependency tree

Central charges `c(G₂) = 14/5` and `c(F₄) = 26/5` sum exactly to `c(E₈) = 8`. The adjoint branches as `248 → (14, 1) ⊕ (1, 52) ⊕ (7, 26)`. The chain is realized as:

```
E₈(1) → G₂(1) × F₄(1)
  │
  ├── G₂ sector → leptons        (Brannen Z₃: B/A = √d₁₀ = √2, θ = d₁₀/d₁₁² = 2/9, Q₀ = d₁₀/d₁₁ = 2/3)
  ├── F₄ sector → quarks + v_EW + Higgs
  │                              v_EW = M_Pl·exp(−(9π²/2 − 6 + 15/512 − 16(α/2π)²))
  │                              α(emergence) = π/512,  α(0) = 1/137.035999050
  │                              λ(M_Pl) = −N_bridge·α_G₂²·E[v²]·(1 − h₁₀) → m_H = 125.30
  ├── SU(3)₃ / D⁽⁶⁾ → CKM        (λ = tan 2/9, A = √(2/3), η̄ = π/9, ρ̄ = √2/9)
  ├── SU(3)₃ / Z_C → PMNS        (tribimaximal + corrections, predicts δ_CP = 76.9°)
  ├── F₄ singlets → ν masses     (m₁ = 0 structurally, normal ordering)
  ├── Embedding index → α_s
  ├── G₂ instanton → η_B
  ├── (7, 26) bridge → gravity   (182 modes, ξ = 1/(48π))
  └── Bridge self-interference → dark sector  (Ω_DM/Ω_b = 2π − 1)
```

## Architecture

The package lives in [`interference/`](interference/). Each derivation module exposes one `derive(...) → dict` function; [`interference/run.py`](interference/run.py) calls them in dependency order.

| Module | Role |
|:---|:---|
| `root.py` | Algebraic root: four integers `(d₁₀, d₁₁, n₇, n₂₆)`, embedding, Casimirs, central charges, the `EchoTerm`/`Ledger`/`Web` echo machinery (Hurwitz depth gate), `M_Pl` from the `m_e` anchor (Newton's `G` becomes a prediction), PDG reference data |
| `masses.py` | All nine masses from `m_e` anchor: Brannen Z₃ leptons, F₄ quarks, `v_EW` instanton suppression |
| `mixing.py` | CKM (four Wolfenstein parameters from four WZW faces) and PMNS (conjugation invariant, tribimaximal + corrections) |
| `octonions.py` | Octonionic `G₂/SU(3)` Clebsch-Gordan verification: `|C₁| = 1`, `|C_3̄| = √2` from Fano plane |
| `couplings.py` | `α_s` via exact WZW cancellation at the gauge-matching scale `μ* = M_Pl·e^−(9π²/2 − 6) ≈ 253.5 GeV` + G₂ → SU(3) threshold; bridge self-interference `α(0)`: depth-1 `π²/(256(2π − 1))`, closed at depth 3 to `1/α(0) = (512/π)(1 − 1/(2π) − 2(α/2π)²) = 137.035999050` |
| `gravity.py` | Sakharov + Jacobson induced gravity, `(7, 26)` heat kernel, protected forgetting, electroweak chain (`sin²θ_W`, `M_W`, `M_Z`), Higgs mass, `η_B`, neutrinos |
| `dark_sector.py` | `Ω_DM/Ω_b = 2π − 1` from bridge venting |
| `words.py` | Generation word lemma: boundary-walk word counts on the D⁽⁶⁾ nimrep `(4, 9, 12, 97)` `→` the integer bases of the quark mass ratios `(38/9, 83/9, 217/18, 1165/12)` |
| `embedding_uniqueness.py` | Exhaustive proof that `E₈(1) ⊃ G₂(1) × F₄(1)` is the unique conformal embedding passing all six gates |
| `protected_forgetting.py` | Verifies `PvP = 0`, `Pv²P = 1/2` for the `G₂` harmonic `v = (1, −½, −½)` and its Weyl orbit uniqueness |
| `nls_soliton.py` | Z₃-coupled NLS simulation: BdG linearization, soliton stability, number conservation, Madelung sourcing |
| `registry.py` | The ledger as code: watches, predictions, promotions, derivation programs, closure record, epistemic record, papers. Machine-checked against the engine by `tests/test_registry.py` |
| `run.py` | End-to-end runner, prints scorecard and emergence tree |

Tests live at the repo root in [`tests/`](tests/) (standard layout):

| File | Role |
|:---|:---|
| `tests/conftest.py` | Session-scoped fixture that runs the full derivation chain once, silently |
| `tests/` | 116 tests: the scorecard freeze table (`test_interference.py`) + the mechanism probes (`test_probes.py`; 6 substrate-dynamics tests behind `-m slow`) + the look-elsewhere, cross-bracing, grammar-verification, and derivation-program suite (`test_coverage.py`) + the registry-engine coupling (`test_registry.py`) |

## Install and run

```bash
git clone https://github.com/kuwrom/one-field.git
cd one-field
pip install -r requirements.txt
python interference/run.py    # full derivation, prints scorecard and emergence tree
pytest -q                      # 116 tests (110 fast + 6 slow behind -m slow)
```

Requires Python 3.10+ and NumPy. Runs in about twenty seconds, dominated by the NLS soliton simulation (`nls_soliton.py`) and the Planck-to-electroweak RGE integration in `gravity.py` (20 000-step RK4, iterated three times for self-consistency).

## Tests

The scorecard is encoded as a `pytest` suite in [`tests/`](tests/). 116 tests (110 fast, 6 slow substrate-dynamics) cover the 35 numerical checks, the 3 structural predictions, the algebraic identities that link them (polynomial relations, Casimir values, central-charge sums, octonionic CG, nimrep eigenvalues, freeze-table integers), and, in `test_coverage.py`, the framework's own epistemics: sparsity of every published enumeration menu (a hit must not be cheap), cross-bracing identities (one constant, many consumers), closed-form recomputation of every FORCED echo multiplicity from the dictionary, and the joint η_B ↔ Ω_DM cosmology closure.

```
$ pytest -q
........................................................................
110 passed, 6 deselected in ~25s
```

If any modification to the algebra moves a prediction outside its tolerance band, the suite fails. This is the operational form of the closure claim *("if any of the 35 numerical checks broke, the structure would break with it")* reduced to one command. Note the tolerance bands are declared per observable in `tests/test_interference.py` and vary in strictness (ppm-level for `v_EW`, ±1% for `m_H`, ±2σ Planck for `Ω_DM/Ω_b`); the claim "passes" is always relative to the declared band, and the bands are part of the audit surface.

## Version history

- **v0.2** (current): polynomial closure on four integers + recursive echo ledger. Eleven `derive()` modules plus the `run.py` driver in `interference/`, `m_e` as the sole dimensional anchor, `M_Pl` and `G` derived. Companion paper: *The Octavian Coherence Gate* ([doi:10.5281/zenodo.20493955](https://doi.org/10.5281/zenodo.20493955)). 35 numerical checks, 116 tests.
- **v0.1.1**: final state of the per-layer `E8/` package (14 modules, layer-by-layer architecture). Preserved at tag [`v0.1.1`](https://github.com/kuwrom/one-field/releases/tag/v0.1.1).

## The registry

The ledger is code: [`interference/registry.py`](interference/registry.py) carries the seven falsifiable watches (each with a stated kill condition), sixteen registered prediction exposures, the promotion record (the m_H closure, first edge admitted under the depth pre-registration rule), the three derivation programs (remaining proof obligations, pre-committed with named components and falsification points), the closure record, and the epistemic record. `tests/test_registry.py` asserts every registered value against the solved engine, so a ledger entry that stops matching canon fails CI. Render it with `python interference/registry.py`.

## Contributing

Verification, critique, tooling, and candidate derivations meeting the registry's promotion bar (a committed edge menu and stated kill conditions) are welcome; corrections with anything adjustable in them are not, because there is nothing here to tune. Changes to canonical values are tracked in git history and enforced by the registry tests.

## Papers

Citation info and the motivation papers (with DOIs) live in [`interference/registry.py`](interference/registry.py) (`PAPERS`, `CITE`) and in [`CITATION.cff`](CITATION.cff). The repository supersedes the papers as the living canon.

## License

MIT. See [`LICENSE`](LICENSE).
