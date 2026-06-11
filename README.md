# one-field

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20323890.svg)](https://doi.org/10.5281/zenodo.20323890)
[![tests](https://github.com/kuwrom/one-field/actions/workflows/tests.yml/badge.svg)](https://github.com/kuwrom/one-field/actions/workflows/tests.yml)

**End-to-end derivation of the Standard Model, Newton's constant, and the dark/baryon ratio through the conformal embedding `E₈(1) ⊃ G₂(1) × F₄(1)`.**

**Zero dimensionless free parameters.** Every prediction in this package is a polynomial or rational function of four irreducible integers plus `π`:

```
d₁₀ = 2   quantum dimension of SU(3)₃ fundamental (1, 0)
d₁₁ = 3   quantum dimension of SU(3)₃ adjoint     (1, 1)
n₇  = 7   dim G₂ fundamental (Fano-plane lines)
n₂₆ = 26  dim F₄ fundamental (traceless Albert algebra)
```

**`m_e` is the sole dimensional input.** The electron mass anchors the scale (the first closure anchors the ruler). `M_Pl`, `G`, `v_EW`, `m_μ`, `m_τ`, and every other dimensional quantity are *outputs* of the forced chain.

**Deeper layers compose through one recursion.** Every quantity is a `Ledger` node (a base amplitude plus a stack of `EchoTerm` ripples), composed by the kernel

```
x_{t+1} = b + W(x_t)
```

where `W` is the derived-weight edge map. Every edge is a theorem (a polynomial in the four integers, or a known phase), so no numerology can enter through the recursion itself. The **Hurwitz gate** caps echo depth at 3, matching the Cayley–Dickson tower `R → C → H → O`. The cap is on the bookkeeping, not on the substrate. If an observable closes within three echo layers (as `α(0)` does at ~10⁻⁹), that is *evidence* the standing wave folds no deeper. The four-number rigidity proof and the depth cap are derived in *The Octavian Coherence Gate* (companion paper, below).

**35 numerical checks plus 3 structural closure predictions.** Run it yourself:

```bash
python interference/run.py    # full derivation
pytest -q                      # 58 tests including the canonical freeze table
```

## The scorecard

All numbers below come from running the code as-is.

### Masses (9/9 within 1%, max error 0.3%)

| Particle | Predicted | PDG | Error | Method |
|---:|---:|---:|---:|:---|
| e | 0.5109990 MeV | 0.5109989 MeV | anchor | CODATA, sole dimensional input |
| μ | 105.6584 MeV | 105.6584 MeV | +0.0000% | Brannen Z₃ (predicted) |
| τ | 1776.909 MeV | 1776.930 MeV | −0.0012% | Brannen Z₃ (predicted) |
| u | 2.158 MeV | 2.16 MeV | −0.1% | (38/9) m_e |
| d | 4.712 MeV | 4.70 MeV | +0.3% | (83/9) m_e |
| s | 93.2 MeV | 93.5 MeV | −0.3% | Rivero + bridge |
| c | 1273.8 MeV | 1273 MeV | +0.06% | (217/18) m_μ |
| b | 4193.8 MeV | 4183 MeV | +0.3% | Q = Q₀ + h₁₁/K³ |
| t | 172.51 GeV | 172.57 GeV | −0.036% | (1165/12) m_τ |

### Mixing, couplings, gravity, baryogenesis, dark sector

| Sector | Result | Reference | Status |
|:---|:---|:---|:---|
| Electroweak scale `v_EW` | 246.21 GeV | 246.22 GeV | −0.00% |
| Higgs mass `m_H` | 124.85 GeV | 125.20 GeV | −0.28% |
| CKM (15 observables) | χ²/n = 0.69 | PDG 2025 | max pull 1.4σ |
| PMNS (3 observables) | χ²/n = 0.00 | NuFit 6.0 | predicts δ_CP = 76.9° |
| Strong coupling `α_s(M_Z)` | 0.1177 | PDG 0.1180 | −0.27% |
| Weinberg angle `sin²θ_W` | 0.231234 | PDG 0.23129(4) | pull −1.41σ |
| W mass `M_W` | 80.365 GeV | 80.369 GeV | pull −0.33σ |
| Z mass `M_Z` | 91.196 GeV | 91.188 GeV | +0.0088% |
| Newton's constant `G_ind/G_N` | 1.0000 (UV) – 1.0134 (broken) | 1 | 0.0–1.3% |
| Baryon asymmetry `η_B` | 6.177 × 10⁻¹⁰ | Planck 6.12 × 10⁻¹⁰ | +0.9% |
| Dark / baryon ratio `Ω_DM/Ω_b` | `2π − 1` = 5.2832 | Planck 5.364 | pull −1.26σ |

**Structural predictions:** `m_ν₁ = 0` (rank-2 seesaw, two RH neutrinos in F₄), normal ordering, `θ_QCD = 0` (forced by `π₃(G₂) → π₃(SU(3))` with the unique E₈(1) vacuum module), and the cosmological-constant scale `ρ_Λ ~ M_Pl² H₀²` (Volovik–Jacobson–CKN, the 10¹²³ problem structurally resolved).

**The precision hierarchy is structural.** Leptons sit closest to the `Z₃` source and emerge cleanly (`m_μ`, `m_τ` predicted from `m_e` to 0.0012%). Quarks pass through additional layers (triality, generation scaling, confinement), and each layer adds residuals of order `α_s/π` (≤ 0.3%). Gravity sits two layers further out, at 0.0–1.3%. The accuracy gradient is not noise. It is what the emergence depth predicts.

## The picture

### There is a substrate

To be is, minimally, to be self-present: anything that exists already affects itself by being there. Interference is the physical image of this: a single photon traversing both slits is a substrate whose fluctuations reference their own history. This framework doesn't model what the substrate is made of. It tracks what the substrate *does* when it references itself.

### Why octonions

A real number records the amplitude of a wave. A complex number records the amplitude and its phase. A quaternion records the phase of a phase. An octonion adds one more layer: **the phase of a phase of a phase**. These aren't what the wave "is made of". They're the accounting ledger for tracking how interference patterns reference their own history. The octonion algebra is the deepest ledger that stays consistent. A fourth doubling would introduce zero divisors and the bookkeeping would contradict itself (Hurwitz's theorem). The automorphism group of the octonions is the Lie group `G₂ = Aut(𝕆)`, and it labels exactly the interference patterns that form closed loops through all three layers of this self-referencing record.

### What mass is

The substrate self-interferes. Against vanishing odds, **interference loops emerge**: `a → b → c → a`. A loop carries cyclic `Z₃` symmetry, because nothing inside it distinguishes which position came first. **Most loops just cycle.**

Some, though, **drop the bias** that singles out one of their three positions. The asymmetric label `v = (1, −½, −½)` inherited from the parent `G₂` structure is erased, and all three positions become equivalent. A forgetting loop can no longer drift between its positions. It gets **stuck in one of three eigenstates**. That stuck pattern is a standing wave, and the three eigenstates are the three lepton generations: electron, muon, tau. Same loop, three places it can lock in.

**Mass is that pattern.**

The lepton mass formula is fixed entirely by the algebra of the loop and the octonionic `G₂/SU(3)` Clebsch–Gordan data: `B/A = √d₁₀ = √2`, `θ = d₁₀/d₁₁² = 2/9`, and the Koide identity `Q₀ = d₁₀/d₁₁ = 2/3` is the algebraic shadow of this closure. `m_μ` and `m_τ` are then determined by `m_e` alone.

When `G₂` breaks to `SU(3)`, the three-label structure crystallizes into three lepton generations. But `G₂` doesn't sit alone. There is exactly one place it can sit, the unique exceptional conformal embedding `E₈(1) ⊃ G₂(1) × F₄(1)`. The `F₄` companion gives the rest of matter.

The algebraic distinction is sharp. `G₂ = Aut(𝕆)` is the automorphism group of a single octonion, and the `G₂` factor carries the leptons: **single-direction fluctuations**. `F₄ = Aut(J₃(𝕆))` is the automorphism group of the Albert algebra of 3×3 Hermitian octonion matrices, and the `F₄` factor carries the quarks: **fluctuations entangling all three triality sectors of the Albert algebra simultaneously**. Leptons are simple. Quarks are woven.

This shows up in the energy scale. Leptons sit in the `F₄` singlet (Casimir = 0), so their standing wave closes at the confinement scale `Λ_conf ≈ 314 MeV`. Quarks sit in the **`26` of `F₄`** (Casimir = 6), and that Casimir shift moves their scale to the electroweak `v_EW = M_Pl · exp(−(9π²/2 − 6 + 15/512)) ≈ 246 GeV`. The factor of ~784 between `Λ_conf` and `v_EW` is the exponential of an algebraic constant, not a hierarchy that needs tuning.

`Spin(8)` triality inside the `26` splits matter cleanly: `26 → 8_v ⊕ 8_s ⊕ 8_c ⊕ 2·1` becomes charged leptons, up-type quarks, down-type quarks, and two right-handed neutrinos. The same `Z₃` that gave three lepton generations gives three of each quark type, woven through triality.

Light quark masses follow a clean `1:2:3` ratio with the electron, set by the `Z₃` centre charges `{0, 1/3, 2/3, 1}` of `{ν, d, u, e}`. Heavy quarks follow `SU(3)₃` WZW emergence: `m_t = (1165/12) m_τ`, `m_c = (217/18) m_μ`, with `Q(c, b, t) = 289/432`. Same closure logic, deeper layer.

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

But forming a standing wave is not enough. The shadow has to be coherent with everyone else's, or it will not add up. The gravity paper tests seven candidate branches against three independent gates (lepton closure, quark closure, gravity closure). Five form standing waves whose shadows are misshaped, or whose bridge reading reverses sign and makes gravity repulsive. **Only one branch closes all three gates: protected `G₂` forgetting applied exactly once.** This is the falsification gauntlet.

The screening length is the Planck healing length `ξ₀ = ℓ_Pl/2`, so the bare scalar channel is suppressed by `exp(−r/ξ_Pl) < 10⁻³⁰⁰` at any observable distance (invisible). What survives macroscopically is the induced spin-2 sector: matter standing waves propagating on the common acoustic metric.

Two independent routes converge on the Einstein equations with the same coefficient: a one-loop heat-kernel (Sakharov) calculation of the induced action, and a Jacobson–Clausius construction from horizon thermodynamics. Newton's `G` falls out at 0.0–1.3%.

The same machinery fixes the cosmological constant. The naive QFT estimate of the vacuum energy is 10¹²³ times the observed value. **The bulk vacuum energy vanishes exactly by the Gibbs–Duhem identity** (Volovik's mechanism), not by fine-tuning. Only the de Sitter departure from equilibrium sets the scale.

The Volovik–Jacobson–CKN chain gives `ρ_Λ = (3/8π) M_Pl² H₀² ≈ 3.7 × 10⁻⁴⁷ GeV⁴`, where `H₀` plays the same role for the Λ sector that `M_Pl` plays for matter: a measured reference, not a free parameter. The observed `ρ_Λ ≈ 2.5 × 10⁻⁴⁷` is lower than the equilibrium prediction by `1/Ω_Λ ≈ 1.46`, because the universe has not yet reached asymptotic de Sitter. **The 10¹²³ cosmological-constant problem is structurally resolved.** The residual factor is cosmological, not algebraic.

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

Every degree of freedom of `E₈` has a physical job. Nothing is left over. At each layer, the observables are fully determined by that layer's algebraic data: leptons on `G₂`, quarks and mixing on `G₂ × F₄` via the embedding, gravity on the mixed bridge. **No free parameters at any layer, and no parameters inherited from above.** This is the closure, and the falsification gate. If any of the 35 numerical checks broke, the structure would break with it. They don't.

## The dependency tree

Central charges `c(G₂) = 14/5` and `c(F₄) = 26/5` sum exactly to `c(E₈) = 8`. The adjoint branches as `248 → (14, 1) ⊕ (1, 52) ⊕ (7, 26)`. The chain is realized as:

```
E₈(1) → G₂(1) × F₄(1)
  │
  ├── G₂ sector → leptons        (Brannen Z₃: B/A = √d₁₀ = √2, θ = d₁₀/d₁₁² = 2/9, Q₀ = d₁₀/d₁₁ = 2/3)
  ├── F₄ sector → quarks + v_EW + Higgs
  │                              v_EW = M_Pl · exp(−(9π²/2 − 6 + 15/512))
  │                              α(emergence) = π/512,  α(0) = 1/137.04
  ├── SU(3)₃ / D⁽⁶⁾ → CKM        (λ = tan 2/9, A = √(2/3), η̄ = π/9, ρ̄ = √2/9)
  ├── SU(3)₃ / Z_C → PMNS        (tribimaximal + corrections, predicts δ_CP = 76.9°)
  ├── F₄ singlets → ν masses     (m₁ = 0 structurally, normal ordering)
  ├── Embedding index → α_s
  ├── G₂ instanton → η_B
  ├── (7, 26) bridge → gravity   (182 modes, ξ = 1/(48π))
  └── Bridge self-interference → dark sector  (Ω_DM/Ω_b = 2π − 1)
```

## Architecture

The package lives in [`interference/`](interference/). Each module exposes one `derive(deps) → dict` function. The order is explicit in [`interference/run.py`](interference/run.py).

| Module | Role |
|:---|:---|
| `root.py` | Algebraic root: four integers `(d₁₀, d₁₁, n₇, n₂₆)`, embedding, Casimirs, central charges, PDG reference data |
| `masses.py` | All nine masses from `m_e` anchor: Brannen Z₃ leptons, F₄ quarks, `v_EW` instanton suppression |
| `mixing.py` | CKM (four WZW structures) and PMNS (conjugation invariant, tribimaximal + corrections) |
| `octonions.py` | Octonionic `G₂/SU(3)` Clebsch–Gordan verification: `|C₁| = 1`, `|C_3̄| = √2` from Fano plane |
| `couplings.py` | `α_s`, `α(0) = π²/(256(2π − 1))`, bridge self-interference, electroweak imports, Higgs |
| `gravity.py` | Sakharov + Jacobson induced gravity, `(7, 26)` heat kernel, protected forgetting, `η_B`, neutrinos |
| `dark_sector.py` | `Ω_DM/Ω_b = 2π − 1` from bridge venting |
| `words.py` | Generation word lemma: D⁽⁶⁾ nimrep eigenvalues `→` integer mass ratios `(38/9, 83/9, 217/18, 1165/12)` |
| `run.py` | End-to-end runner, prints scorecard and emergence tree |

Tests live at the repo root in [`tests/`](tests/) (standard layout):

| File | Role |
|:---|:---|
| `tests/conftest.py` | Session-scoped fixture that runs the full derivation chain once, silently |
| `tests/test_interference.py` | 58 tests including the canonical freeze table |

## Install and run

```bash
git clone https://github.com/kuwrom/one-field.git
cd one-field
pip install -r requirements.txt
python interference/run.py    # full derivation, prints scorecard and emergence tree
pytest -q                      # 58 tests
```

Requires Python 3.10+ and NumPy. Runs in well under a second.

## Tests

The scorecard is encoded as a `pytest` suite in [`tests/test_interference.py`](tests/test_interference.py). 58 tests cover the 35 numerical checks, the 3 structural predictions, and the algebraic identities that link them (polynomial relations, Casimir values, central-charge sums, octonionic CG, nimrep eigenvalues, freeze-table integers).

```
$ pytest -q
..........................................................
58 passed in 4.56s
```

If any modification to the algebra moves a prediction outside its tolerance band, the suite fails. This is the operational form of the closure claim *("if any of the 35 numerical checks broke, the structure would break with it")* reduced to one command.

## Version history

- **v0.2** (current): polynomial closure on four integers + recursive echo ledger. Eight `derive()` modules plus the `run.py` driver in `interference/`, `m_e` as the sole dimensional anchor, `M_Pl` and `G` derived. Companion paper: *The Octavian Coherence Gate* ([doi:10.5281/zenodo.20493955](https://doi.org/10.5281/zenodo.20493955)). 35 numerical checks, 58 tests.
- **v0.1.1**: final state of the per-layer `E8/` package (14 modules, layer-by-layer architecture). Preserved at tag [`v0.1.1`](https://github.com/kuwrom/one-field/releases/tag/v0.1.1).

## Citing

Companion papers, in dependency order:

1. Kahsay, Kibrom Kidane (2026). *The Innocent Lepton*. Zenodo. [doi:10.5281/zenodo.19899091](https://doi.org/10.5281/zenodo.19899091)
2. Kahsay, Kibrom Kidane (2026). *One Substrate, Three Generations*. Zenodo. [doi:10.5281/zenodo.20069456](https://doi.org/10.5281/zenodo.20069456)
3. Kahsay, Kibrom Kidane (2026). *The Echo of Standing Waves*. Zenodo. [doi:10.5281/zenodo.20144381](https://doi.org/10.5281/zenodo.20144381)

v0.2 algebraic foundation (rigidity proof for the four integers):

4. Kahsay, Kibrom Kidane (2026). *The Octavian Coherence Gate: Rigidity of the Four Irreducible Integers of the E₈(1) ⊃ G₂(1) × F₄(1) Conformal Embedding*. Zenodo. [doi:10.5281/zenodo.20493955](https://doi.org/10.5281/zenodo.20493955)

## License

MIT. See [`LICENSE`](LICENSE).
