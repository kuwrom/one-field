# one-field

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20323890.svg)](https://doi.org/10.5281/zenodo.20323890)
[![tests](https://github.com/kuwrom/one-field/actions/workflows/tests.yml/badge.svg)](https://github.com/kuwrom/one-field/actions/workflows/tests.yml)

To be is, minimally, to be self-present. Anything that exists already affects itself by being there. Interference is the physical image of this. What is said to exist is what carries non-zero interference weight. The framework observes what emerges from interference.

**A substrate self-interferes. Exactly one family of standing wave patterns closes to give our universe: the Standard Model, gravity, and the dark sector.**

Concretely, the code starts from one dimensional anchor (the electron mass) and the closure demand, whose unique solution is the algebra. From these it produces the other eight fermion masses, the CKM and PMNS scorecards, the electroweak scale, the fine structure and strong couplings, Newton's constant, the baryon asymmetry, and the dark matter fraction. It has zero tunable dimensionless parameters. Each result carries an explicit epistemic tag: theorem, enumeration, identification hypothesis, coordinate convention, arithmetic observation, or cited import with a stated band. The runner prints the scorecard and the derivation graph. The substrate simulations are the physical realization program. They test whether the structure emerges, and they supply no input to the scorecard.

## Run it

```bash
git clone https://github.com/kuwrom/one-field.git
cd one-field
pip install -r requirements.txt
python3 interference/run.py   # full derivation, prints scorecard and derivation graph
pytest -q                     # 311 tests (295 fast + 16 slow behind -m slow)
```

Python 3.10+, NumPy, SciPy. Twenty seconds to a minute on one core.

## Masses (one anchor in, eight out, max error 0.4%)

| Particle | Predicted | PDG | Error | Method | Status |
|---:|---:|---:|---:|:---|:---|
| e | 0.5109990 MeV | 0.5109990 MeV | anchor | CODATA, sole dimensional input | input |
| μ | 105.6584 MeV | 105.6584 MeV | +0.0000% | Brannen Z₃ (predicted) | Thm + IH |
| τ | 1776.909 MeV | 1776.930 MeV | −0.0012% | Brannen Z₃ (predicted) | Thm + IH |
| u | 2.158 MeV | 2.16 MeV | −0.1% | (38/9) m_e | Thm + IH + conv |
| d | 4.713 MeV | 4.70 MeV | +0.3% | (83/9) m_e | Thm + IH + conv |
| s | 93.8 MeV | 93.5 MeV | +0.4% | Q₀+h₁₀/K³ + bridge | Thm + IH + conv |
| c | 1273.8 MeV | 1273 MeV | +0.06% | (217/18) m_μ | Thm + IH + conv |
| b | 4193.8 MeV | 4183 MeV | +0.3% | Q = Q₀ + h₁₁/K³ | Thm + IH + conv |
| t | 172.51 GeV | 172.57 GeV | −0.036% | (1165/12) m_τ | Thm + IH + conv |

<details>
<summary><b>Fine print: mass coordinates and status labels</b></summary>

**Status labels.** Thm = executable theorem in this repo. IH = identified (a mathematical object matched with a physical one, where no other reading closes, separately falsifiable). conv = stated scheme convention.

The framework produces one completed reading per fermion, with every echo included. It does not run to a scale, so each comparison names the coordinate of the external reference. Leptons are compared at the propagator pole, and the top at the PDG direct mass reference. Confined heavy quarks (c, b) are compared at the self-scale `m(m)`. Light quarks (u, d, s) are compared at the PDG `MS-bar(2 GeV)` coordinate. With the alternative charm reference near 1.67 GeV, the charm residual moves from +0.06% to roughly −24%, with no change to the framework. Likewise, the Koide relation among the measured lepton masses, transported to `M_Z`, holds only at the 0.17% level. The claim is always the completed reading at its stated coordinate, never a relation at every scale. The coordinate free content of the light sector is its RG invariant ratios: `m_u/m_d = 38/83` (−0.9σ vs PDG 0.473(17)), `m_s/m_ud = 27.318` (+0.2σ vs 27.30(8)), `Q_ellipse = 22.383` (+0.4σ dispersive, −0.2σ FLAG lattice). See `masses.py` for the full treatment.

</details>

**Lepton mass chain.** The calculation uses Brannen's 2006 circulant form. The `G₂ → SU(3)` Fano plane calculation gives the Clebsch-Gordan ratio `√2`. The closing lane's `SU(3)₃` modular data have canonical positive spherical dimension `d(1,0) = 2`. The channel reading map combines these data with two cycle orientations and gives `B/A = √2`. The circulant identity then gives `Q₀ = 2/3`. The lightest triality changing primary has Sugawara weight `h = 2/9`, and the mass phase map sets `θ = h`. This is the only positive, nondegenerate value in the committed angle menu. The reading is one ledger that includes every interference layer. The forced depth-4 EM backreaction depends on flavour, so it survives in the mass ratios. It contributes −9.8 ppm to the μ/e ratio, and the measurement resolves it ([`masses.py`](interference/masses.py)). The framework computes the physical particle with its electromagnetic self-interference included. That is why the propagator pole entry is the comparison coordinate. With `m_e` as the dimensional anchor, the calculation reproduces `m_μ` to 0.0000003% and `m_τ` to 0.0012%.

## What mass is

A particle is a closed standing wave of the substrate, whose phase circulates at one fixed frequency. In the empty field it is realized as a self-trapped lump ([`lump_in_empty_field.py`](tests/probes/lump_in_empty_field.py)). That frequency is the mass, `E = ħω`. A heavier particle is "heavy" only in that its internal clock is faster.

The mass construction has three clock settings because it uses three child labels with `Z₃` cyclic symmetry (`a → b → c → a`). A graph can of course also have self-loops or two cycles (the Web uses both), so three is a framework commitment, not a graph theoretic minimum. When the bias that singles out one child label is erased (protected forgetting), the resulting `Z₃`-circulant has three eigenstates identified with the electron, muon, and tau generations. The circulant fixes the mass ratios. The electron fixes their dimensional unit. `m_μ` and `m_τ` are then determined by `m_e` alone.

Any configuration that closes to give the leptons also produces quarks and gravity, from the same structure. The universe does not have to get lucky twice. The lepton standing wave closes without involving the embedding and sits at the confinement scale. The quark standing wave closes through the embedding and sits at the electroweak scale. The three generations of each come from the same `Z₃`, not from triality. The mass sector is one generator plus a grammar: the lepton circulant is the generator, and the quark masses are arithmetic on it. The full mass chain is in [`interference/__init__.py`](interference/__init__.py).

## Why octonions

The searchlight is the octonionic algebra used to locate the structure. It is organized around the three label closure. Its starting premise is that a closed record is a finite dimensional unital real composition algebra: records compose bilinearly, and their positive definite norm, read as probability weight, is multiplicative. Hurwitz's classification then leaves exactly `R`, `C`, `H`, and `O`. Relative to the real amplitude baseline, the three nonreal stages `C`, `H`, and `O` give searchlight depths 1–3. The next Cayley–Dickson algebra fails the criterion because it has zero divisors. This is a complete derivation given the framework's ledger premise. It is not a claim that Hurwitz's theorem, on its own, says anything about physical self-reference. [`searchlight/octonions.py`](interference/searchlight/octonions.py) verifies bilinearity, norm multiplicativity, the explicit sedenion zero divisor, and the Fano/MTC channel map. The canonical Web then carries an explicit depth-3 `e↔q` loop ledger. Every depth-4 backreaction names that ledger as its parent and reads it.

## Mixing, couplings, gravity, baryogenesis, dark sector

| Sector | Result | Reference | Status |
|:---|:---|:---|:---|
| Electroweak scale `v_EW` | 246.219645 GeV | 246.219640 GeV (from G_F) | +0.02 ppm |
| Fine structure `1/α(0)` | 137.035999050 | Berkeley Cs 137.035999046(27) | +0.13σ vs Cs, **−14σ vs Rb 2020**. The framework commits to the Cs side of the 5.5σ dispute. A resolution toward Rb kills it |
| Higgs transport | `lambda(m_t) = 0.12651` | Internal transport | Candidate bridge boundary in the imported beta coordinate. The depth-4 Web reading is unassigned |
| CKM (15 comparison coordinates, correlations unmodeled) | diagonal pull score/coordinate = 0.21 | PDG 2024 | max marginal pull 0.7σ. The reference block mixes correlated global fit outputs with a separately quoted phase constraint. Without a joint covariance model this is a descriptive score, not χ². Where mixing lives is itself an executable theorem. Two `Z₃`-equivariant transports commute and give `V = 1` exactly, so mixing is not a sector transport. It enters as four readings of the solved `D⁽⁶⁾*` boundary, one per face. These are fixed engine constants with no continuous freedom ([`mixing_z3_dichotomy.py`](tests/probes/mixing_z3_dichotomy.py)). Each reading has its own grade. `η̄ = π/9`, the braiding monodromy, is computed. `λ` is enumerated: the overlap's Hermitian pair closes it on its 12-entry menu without data ([`lambda_admissibility.py`](tests/probes/lambda_admissibility.py)), and it is the unique 1σ hit. The attachment to the physical Wolfenstein coordinates is identified. The per face grading is in [`mixing.py`](interference/mixing.py) |
| PMNS (3 correlated angles) | diagonal pull score/angle = 0.00 | Frozen NuFit 6.0 NO snapshot, IC19 without SK atmospheric data | Identified conjugation orbit TBM map with `U_e = R₂₃(−φ)R₁₃(θ_C,δ)`. The construction uses no new engine numbers. The TBM block is computed from the conjugation orbit compressed fusion matrices. The backreaction terms are functions of the derived mass spectrum, `φ = √(m_e/m_μ)/√d₁₀`, and the echoed `λ`. Symmetric marginal errors and no covariance make this descriptive rather than χ². The phase `−arg(U_e3) = −77.0°` depends on convention and is not reported as the physical `δ_CP`. The rephasing invariant `δ_CP = 275.74°` lies inside all NuFIT 6.0 3σ ranges. Two data blind gates fix its sign: the lepton sector's own pairing and the charge dictionary. The NuFIT ranges confirm it ([`pmns_sign_bits.py`](tests/probes/pmns_sign_bits.py)) |
| Strong coupling `α_s(M_Z)` | 0.1184 | PDG 0.1180(9) | +0.4σ (π/32 exact at μ* = 253.5 GeV) |
| Weinberg angle `sin²θ_W` | 0.231285 | PDG 0.23129(4) | pull −0.11σ |
| W mass `M_W` | 80.356 GeV | 80.3625(77) world avg, without CDF II | pull −0.86σ |
| Z mass `M_Z` | 91.189 GeV | 91.188 GeV | +0.0010% |
| Newton's constant `G` | 6.673936×10⁻¹¹ | CODATA 6.67430(15)×10⁻¹¹ | −2.43σ. `M_Pl` comes from the anchor inversion `m_e → M_Pl` (+22 ppm). `G` reads it through the common mode (the field component that carries no label): `G M_b² = 3Q_l − 1`, the lepton circulant's Koide moment at depth four. One of nine readings of this count survives ([`common_mode_count.py`](tests/probes/common_mode_count.py)). IH |
| Baryon asymmetry `η_B` | 6.178e-10 | Planck 6.12e-10 | +0.9% |
| Dark / baryon ratio `Ω_DM/Ω_b` | `2π − 1` = 5.2832 | Planck 5.364 | pull −1.26σ |

**Joint cosmology.** The two cosmology channels, `η_B` from the G₂ instanton and `Ω_DM/Ω_b` from bridge venting, share one engine input: the normalized current response `h/2π`. `η_B = 6.178e-10 → Ω_b h² = 0.02256` (+1.2σ vs Planck), and `Ω_DM h² = (2π−1) · Ω_b h² = 0.1192` (−0.7σ vs Planck 0.1200(12)). The BBN conversion η₁₀ = 273.9 Ω_b h² is a cited change of coordinate applied to a printed number. No engine step uses it.

**Conditional structural constraint and cosmology bookkeeping.** Suppose the two F₄ singlets are right handed neutrinos, and the mass mechanism is a minimal type-I seesaw with no other light neutrino mass operators. Then the light mass matrix has rank at most two, so at least one state is massless. The ordering is not selected (`m₁ = 0` for normal or `m₃ = 0` for inverted). Separately, `(3/8π) M_Pl² H₀²` is the Friedmann critical density, not `ρ_Λ`. The current dark energy density is `Ω_Λ` times it. Jacobson's thermodynamic derivation leaves Λ as an integration constant, and CKN saturation is an identification hypothesis.

## What gravity is

Mass, above, is a clock. Gravity is the effect of that ticking on the shared background, seen from outside. In a standing wave the substrate circulates in place. That circulation takes up room in the shared background and depletes it. The knot and its depletion form together (`tests/probes/skyrmion_3d.py`). **The soliton casts a shadow**: a depleted region around it. A depletion is a deficit, and a deficit pulls. Every standing wave depletes the same shared background, and the depletions add. A planet's worth of circulations makes one large depletion, a single venting flow shared by every particle in the planet. Having mass and sourcing gravity are the same event.

**Gravity is not a fifth force added on top of matter. It is the depletion that matter makes in the background as it vents, and the bias that other matter follows.** The full treatment is in [`gravity.py`](interference/gravity.py) and [`interference/__init__.py`](interference/__init__.py).

The sign is checked from both directions. The minimal matter spectrum alone gives the wrong sign: `Σ_matter = −10/3` in the standard combination `(N₀ + N_½ − 4N₁)/6`, which is repulsive gravity. The 182 bridge channels make the stiffness attractive. The wiring winner forces their count, and each enters minimally coupled at `1/6`. The total is `Σ = 27`, and `25.83` on the one boundary list that the running `b₀` also reads ([`boundary.py`](interference/boundary.py)). The ledger fixes the sign and the Planck order of the stiffness, not its coefficient: `6π/Σ` moves with the regulator by a factor of two ([`sakharov_normalisation.py`](tests/probes/sakharov_normalisation.py)). The strength of gravity is the anchor inversion read through the common mode. The gravity lane has six gates: sign, universality, magnitude, one cone, tensor carrier, and no preferred location. Each names the tests that derive it and the case in which it fails ([`gravity.py`](interference/gravity.py), [`test_gravity_gates.py`](tests/common_mode/test_gravity_gates.py)).

## What closes

A substrate self-interferes, and many patterns form. Most do not close. Some lanes (candidate closure routes) make electron like objects whose ratios do not carry on to the next sector. Others make quark like objects that never bind. Six consistency gates eliminate them. Only one lane passes all six: protected `G₂` forgetting applied exactly once. The gates and their scan are the searchlight's bookkeeping, not a mechanism. The substrate either closes or does not. The scan is the computed record that no other lane passes. It runs on Lie table invariants alone, and it prints its one real contender, `E₈ ⊃ E₇ × A₁`, beside the winner ([`embedding_uniqueness.py`](interference/searchlight/embedding_uniqueness.py)). The seven branch gauntlet is reproduced from the gravity record (The Echo of Standing Waves, Table V). Its winning row is the executable one. The 504 hits that the closure search reports are mass grammar realizations within this fixed lane, not additional closure lanes.

That surviving lane sits inside one algebra:

```
14  +  52  +  182  =  248  =  dim(E₈)          E₈(1) ⊃ G₂(1) × F₄(1)
gauge   matter   bridge                         c: 14/5 + 26/5 = 8, exact
        + Higgs  → gravity                      248 → (14,1) ⊕ (1,52) ⊕ (7,26)
```

Every degree of freedom of `E₈` has a physical job. Nothing is left over. The audit is in [`interference/__init__.py`](interference/__init__.py).

The same rule works at every scale. Particles, atoms and orbits last because they are closed ([`spinning_lump.py`](tests/probes/spinning_lump.py), [`dimension_by_closure.py`](tests/probes/dimension_by_closure.py)), and solar systems and galaxies are the structures that closed and held together. What we see is what closed.

## Quantum numbers

The Standard Model quantum numbers are read off the four integers by one dictionary, an identification hypothesis stated in [`root.py`](interference/root.py):

```
Q_u  = d₁₀/d₁₁ = 2/3         (up-type electric charge)
Q_d  = 1/d₁₁   = 1/3         (down-type electric charge)
N_c  = d₁₁     = 3           (colour multiplicity)
sin²θ_W = d₁₁/(d₁₀²+d₁₁²) = 3/13   (Weinberg angle, zero edge base)
```

The charge trace `Σ Q² N_c = (d₁₀²+1+d₁₁)/d₁₁ = 8/3` gives the Singh ratio `(8/3)·C₂(26) = 16`. This ratio converts the `G₂` coupling at the end of its lever arm, `α_G₂(μ*) = π/32`, to the zero edge electromagnetic base (the value before any echo edge): `α_EM = α_G₂/16 = π/512`.

**CP violation** comes from the oriented Fano plane: its seven directed lines select one orientation sign. Under the declared map from Fano orientation to phase sign, reversing the orientation flips `η̄`, `J_CKM`, `J_lep`, and `η_B` together, and leaves CP-even quantities unchanged ([`tests/probes/orientation_bit.py`](tests/probes/orientation_bit.py)). There is one global CP orientation bit, not several independently inserted signs.

**Strong CP**: the index one `G₂ → SU(3)_c` embedding preserves winding. In the framework's ontology, which has no Lagrangian, the physical angle is the complete matter coupled character `θ̄`. The full record rule makes that character act on the vacuum record algebra `A = 1 + B`. Multiplication leaves only the trivial phase automorphism, so `θ̄ = 0 mod 2π`. This conclusion follows within that ontology, not from local MTC data alone. A measured neutron EDM attributable to `θ̄ ≠ 0` kills it ([derivation](interference/couplings.py), [supporting probe](tests/probes/strong_cp.py)).

## Derivation flow

The scorecard is an executable dependency pipeline, run in order by [`interference/run.py`](interference/run.py):

1. [`wiring.py`](interference/wiring.py) scans the finite candidate menus (lattices, conformal pairs, levels) and finds the unique closure. The depth-3 searchlight reads `d₁₀ = 2`, `d₁₁ = 3`, `n₇ = 7`, `n₂₆ = 26`, and `α_G₂(M_Pl) = 1/(24π)` off the winner. [`root.py`](interference/root.py) uses that vocabulary, and `run.py` asserts that both stages share the same wiring record.
2. [`root.py`](interference/root.py) builds the recursive web, its couplings and scales, and derives `M_Pl` by inverting the electron anchor.
3. [`words.py`](interference/words.py) derives the generation walk counts. [`masses.py`](interference/masses.py) derives the nine fermion masses and `v_EW`.
4. [`mixing.py`](interference/mixing.py) computes the WZW data, exact D⁽⁶⁾* connection, and mixing scorecards. [`couplings.py`](interference/couplings.py) computes `α_s` and `α(0)`.
5. [`gravity.py`](interference/gravity.py) computes the induced gravity, baryogenesis, and neutrino outputs. [`higgs_ew.py`](interference/higgs_ew.py) computes the Higgs quartic boundary and electroweak observables (SM RGE transport, M_W, M_Z). [`dark_sector.py`](interference/dark_sector.py) computes the cosmological ratios.

The substrate is the physical realization of this chain. Its probes test whether the structure emerges. They do not supply scorecard inputs.

## Architecture

The derivation lives in the [`interference/`](interference/) package as one dependency chain. Package relative imports give one canonical module graph and one Web singleton. The checks live in [`tests/`](tests/), one folder per branch (foundations, masses, couplings, mixing, common mode, substrate, audit). Freezes pin the numbers, probes cover mechanisms and falsifiers, and the substrate and audit instruments provide numerical checks. Each derivation module exposes one `derive(...) → dict`, and [`run.py`](interference/run.py) calls them in dependency order. The full narrative follows that order in [`interference/__init__.py`](interference/__init__.py).

| Module | Role |
|:---|:---|
| `run.py` | End to end runner: `python3 interference/run.py` |
| `root.py` | Embedding, Casimirs, central charges, the recursive echo web, `M_Pl` from the `m_e` anchor. The four integers are the depth-3 searchlight's vocabulary, read from `wiring.py` |
| `masses.py` | All nine masses: Brannen Z₃ leptons, F₄ quarks, `v_EW` |
| `mixing.py` | WZW/D⁽⁶⁾ boundary data, braiding phase, CKM parameter map, PMNS conjugation orbit TBM construction |
| `couplings.py` | `α_s` (WZW cancellation at μ*), `α(0)` from bridge self-interference |
| `gravity.py` | Induced gravity ledger, target factorizations, `η_B`, neutrinos |
| `higgs_ew.py` | SM RGE transport, Higgs quartic boundary, M_W, M_Z, import audit |
| `dark_sector.py` | `Ω_DM/Ω_b = 2π − 1` from bridge venting |
| `words.py` | Generation word lemma: nimrep walk counts → quark mass bases |
| `boundary.py` | The rows at the boundary: one list that the running (`b₀`) and the induced stiffness both read |
| `wiring.py` | Closure residual scan over the candidate menus. The four integers are never an input |
| `reference.py` | Comparison data, single home, quarantined by `tests/audit/test_blindness.py` |
| `searchlight/` | The octonions: `G₂/SU(3)` CG verification, embedding uniqueness through six coherence conditions on parent/factor pairs, protected forgetting. The one subdirectory, and the engine never imports it |
| `tests/substrate/` | Z₃-NLS substrate: BdG, stability, Madelung sourcing, and the same equation in 3D |
| `tests/audit/` | The committed grammar, and the same grammar run as a generator: blind gates first, then a single comparison with data |


## Physics FAQ

The derivation probes in [`tests/probes/`](tests/probes/) answer standard physics objections. Each states its own verdict and kill condition. The table below is the discoverability map.

| Question | Answer | Code | Status |
|:---|:---|:---|:---|
| Why 3+1 dimensions? | A bound orbit stays bound only below four space dimensions, and it closes on itself exactly in three. Of the 64 sign patterns, only the Lorentzian ones with one time direction keep the energy bounded below. Structures that both close and persist therefore fix 3+1 | [`dimension_by_closure.py`](tests/probes/dimension_by_closure.py) | Derived |
| Why are solitons stable? | On the committed coupling sign, a particle in the empty field is a self-trapped lump. Its closure is a circulation: one quantum of common phase winding around a ring persists, protected by the conserved charge and angular momentum. A `π₃` texture placed in the lump unwinds, as it does in the bare NLS on the gas sign ([`knot_charge.py`](tests/probes/knot_charge.py)) | [`lump_stability.py`](tests/probes/lump_stability.py), [`spinning_lump.py`](tests/probes/spinning_lump.py) | Computed |
| Is there confinement? | Vortex proliferation is measured, the area law fit beats the perimeter law fit decisively, and the Gaussian excitation spectrum is gapped | [`confined_phase.py`](tests/probes/confined_phase.py) | Classical confinement signature and Gaussian gap computed |
| How does gravity enter? | Written without a background metric, the carrier defines a metric with no action of its own: the gradient bilinear summed over the three components. Its tensor stiffness is the heat kernel weight per state, and the venting ledger sums those weights. The strength of gravity is the anchor inversion | [`gravity.py`](interference/gravity.py), [`composite_metric.py`](tests/probes/composite_metric.py), [`composite_pole.py`](tests/probes/composite_pole.py) | Six gates pass on three statements: no background metric, the self-bound sign, the empty vacuum |
| What about chirality/CP? | Closure leaves the mass protected chiral branch, dimensionality gives 3+1D, and the embedding branch supplies the gauge representations. The radial Skyrmion–Dirac operator is built on the repository's actual `π₃` hedgehog and exhibits the single localized crossing associated with unit winding. Its scanned coupling is a spectral family coordinate, not a fitted or additional framework parameter. Independently, reversing the Fano orientation flips every computed CP-odd sign and preserves CP-even magnitudes. These branches compose under the stated fermion/knot and representation identifications | [`zero_mode.py`](tests/probes/zero_mode.py), [`orientation_bit.py`](tests/probes/orientation_bit.py), [`dimensionality.py`](tests/probes/dimensionality.py) | Complete derivation within the stated identifications. Knot spectral flow computed |
| Are EW quantities secretly fitted? | Import audit: the top contribution uses the framework's own `G_F` and `m_t`. The hadronic VP is classed as measurement data. Higher order pieces are loop mathematics evaluated at the measured `m_H`, which is imported as data. No free theory parameters enter | [`ew_internal.py`](tests/probes/ew_internal.py) | Audited |
| Is the algebraic wiring arbitrary? | `(7,26)` bridge is the unique nontrivial boson of `G₂(1) × F₄(1)`. Face channels confined. Sector graph and cycle inventory derived from the modular tensor category | [`wiring_theorem_probe.py`](tests/probes/wiring_theorem_probe.py) | Derived. The 16 and 2 responses carry the unique grade |
| Can comparison data leak into predictions? | `test_blindness.py` perturbs reference values and requires predictions to remain identical. AST audit of `reference.py` access patterns | [`test_blindness.py`](tests/audit/test_blindness.py) | Enforced |
| Does the knot form in 3D? | The relative mode RMS grows ×15.6 under `G₂` coupling and stays flat in the `g₁ = 0` control. The common mode density dips at the core | [`skyrmion_3d.py`](tests/probes/skyrmion_3d.py) | Measured on the gas sign over a held uniform background |
| What does the classical substrate claim, exactly? | The split between layers: the classical field supplies the slots (a stable loop, three bound modes, channel lifting). The quantized bookkeeping supplies the ratios. A classical measurement of `Q = 2/3` falsifies this split | [`circulant_normalization.py`](tests/probes/circulant_normalization.py) | Committed with falsifier |

## Where the questions end

The engine takes one number, and even that number is a choice of unit, not physics. The ratio `m_e/M_Pl` is itself derived, so the dimensionless output takes no input at all. How strongly the substrate interferes is not a second input either. It is what is being derived, one channel at a time. What the axiom cannot explain is why the amplitude is nonzero at all. Self-presence is universal, not conditional: the empty world is present to itself too, and precisely because it is nothing the self-touching generates no amplitude. That is not an exception to the axiom but its trivial fixed point: the zero orbit of `x = b + W(x)`, zero in and zero out. Nothing does not escape the axiom. It simply never closes into anything further.


## Contributing

Verification, critique, tooling, and candidate derivations meeting the promotion bar (a committed edge menu and stated kill conditions) are welcome. Corrections with anything adjustable in them are not, because there is nothing here to tune.

## Paper

- The Octavian Coherence Gate: The Four Irreducible Integers of the E₈⁽¹⁾ ⊃ G₂⁽¹⁾ × F₄⁽¹⁾ Conformal Embedding ([10.5281/zenodo.20493955](https://doi.org/10.5281/zenodo.20493955)). The four integers and the uniqueness of the embedding. The chain from the integers to the scorecard is this repository.

Related records: The Innocent Lepton ([10.5281/zenodo.19899091](https://doi.org/10.5281/zenodo.19899091)), One Substrate, Three Generations ([10.5281/zenodo.20069456](https://doi.org/10.5281/zenodo.20069456)), The Echo of Standing Waves ([10.5281/zenodo.20144381](https://doi.org/10.5281/zenodo.20144381)).

Cite as: Kahsay, Kibrom Kidane (2026). One-field. https://github.com/kuwrom/one-field, with the Octavian Coherence Gate paper. See also [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).
