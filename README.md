# one-field

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20323890.svg)](https://doi.org/10.5281/zenodo.20323890)
[![tests](https://github.com/kuwrom/one-field/actions/workflows/tests.yml/badge.svg)](https://github.com/kuwrom/one-field/actions/workflows/tests.yml)

To be is, minimally, to be self-present. Anything that exists already affects itself by being there. Interference is the physical image of this. What is said to exist is what carries non-zero interference weight. The framework observes what emerges from interference.

**A substrate self-interferes. Exactly one family of standing-wave patterns closes to give our universe: the Standard Model, gravity, and the dark sector.**

Concretely, from one dimensional anchor (the electron mass) and the closure demand whose unique solution is the algebra, the code produces the other eight fermion masses, the CKM and PMNS scorecards, the electroweak scale, the fine-structure and strong couplings, Newton's constant, the baryon asymmetry, and the dark-matter fraction, with zero tunable dimensionless parameters. Each result carries an explicit epistemic tag: theorem, enumeration, identification hypothesis, coordinate convention, arithmetic observation, or cited import with a stated band. The runner prints the scorecard and derivation graph. The substrate simulations are the physical-realization program: they test emergence and supply no scorecard input.

## Run it

```bash
git clone https://github.com/kuwrom/one-field.git
cd one-field
pip install -r requirements.txt
python3 interference/run.py   # full derivation, prints scorecard and derivation graph
pytest -q                     # 170 tests (155 fast + 15 slow behind -m slow)
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
| t | 172.51 GeV | 172.57 GeV | −0.036% | (1165/12) m_τ | Thm + IH |

<details>
<summary><b>Fine print: mass coordinates and status labels</b></summary>

**Status labels.** Thm = executable theorem in this repo. IH = identification hypothesis (a mathematical object matched with a physical one, separately falsifiable). conv = stated scheme convention. Full taxonomy in [What is derived vs identified](#what-is-derived-vs-identified).

The framework produces one completed reading per fermion, every echo included. It does not run to a scale, so each comparison names the coordinate attached to the external reference. Unconfined fermions (leptons, top) sit at the propagator pole. Confined heavy quarks (c, b) sit at the self-scale `m(m)`. Light quarks (u, d, s) sit at the PDG `MS-bar(2 GeV)` coordinate. The alternative charm reference near 1.67 GeV moves the charm Koide residual from +0.06% to roughly −24% without any change to the framework, and the Koide relation among the measured lepton masses, transported to `M_Z`, holds only at the 0.17% level. The claim is always the completed reading at its stated coordinate, never a relation at every scale. The invariant content of the light sector is its RG-invariant ratios, which carry no coordinate at all: `m_u/m_d = 38/83` (−0.9σ vs PDG 0.473(17)), `m_s/m_ud = 27.318` (+0.2σ vs 27.30(8)), `Q_ellipse = 22.383` (+0.4σ dispersive, −1.7σ lattice, and the two references disagree). See `masses.py` for the full treatment.

</details>

**Lepton mass chain.** The calculation uses Brannen's 2006 circulant form. The `G₂ → SU(3)` Fano-plane calculation gives the Clebsch-Gordan ratio `√2`. The closing lane's `SU(3)₃` modular data have canonical positive spherical dimension `d(1,0) = 2`. The channel-reading map combines these data with two cycle orientations and gives `B/A = √2`. The circulant identity then gives `Q₀ = 2/3`. The lightest triality-changing primary has Sugawara weight `h = 2/9`, and the mass-phase map sets `θ = h`. This is the sole positive, nondegenerate value in the committed angle menu. The reading is one ledger, every interference layer included: the forced depth-4 EM backreaction, flavour-dependent so it survives in the mass ratios, contributes −9.8 ppm on the μ/e ratio, and the measurement resolves it ([`masses.py`](interference/masses.py)). The framework computes the physical particle with its electromagnetic self-interference included, which is why the propagator-pole entry is the comparison coordinate. The Koide relation and Brannen parameterization predate this implementation, so the lepton scorecard is postdictive. With `m_e` as the dimensional anchor, the calculation reproduces `m_μ` to 0.0000003% and `m_τ` to 0.0012%.

## What mass is

A particle is a self-trapped patch of the substrate, a lump whose internal phase circulates at one fixed frequency. That rate is the mass, `E = ħω`. A heavier particle is a faster internal clock, and nothing else about it is "heavy".

The mass construction has three clock settings because it uses three child labels with `Z₃` cyclic symmetry (`a → b → c → a`). A graph can of course also have self-loops or two-cycles (the Web uses both), so three is a framework commitment, not a graph-theoretic minimum. When the bias that singles out one child label is erased (protected forgetting), the resulting `Z₃`-circulant has three eigenstates identified with the electron, muon, and tau generations. The circulant fixes the mass ratios. The electron fixes their dimensional unit. `m_μ` and `m_τ` are then determined by `m_e` alone.

Any configuration that closes to give the leptons automatically produces quarks and gravity from the same structure. The universe does not have to get lucky twice. The lepton standing wave closes without engaging the embedding and sits at the confinement scale. The quark standing wave closes through the embedding and sits at the electroweak scale. Three generations of each are the same `Z₃`, not triality. The mass sector is one generator plus a grammar: the lepton circulant is the generator, and the quark masses are arithmetic on it. The full mass chain is in [`interference/__init__.py`](interference/__init__.py).

## Why octonions

The searchlight is organized around the three-label closure. Its starting premise is that a closed record is a finite-dimensional unital real composition algebra: records compose bilinearly and their positive-definite norm, read as probability weight, is multiplicative. Hurwitz's classification then leaves exactly `R`, `C`, `H`, and `O`. Relative to the real-amplitude baseline, the three non-real stages `C`, `H`, and `O` give searchlight depths 1–3. The next Cayley–Dickson algebra fails the criterion through zero divisors. This is a complete derivation from the framework's ledger premise, not a claim that Hurwitz, in isolation, is a theorem about physical self-reference. [`searchlight/octonions.py`](interference/searchlight/octonions.py) verifies bilinearity, norm multiplicativity, the explicit sedenion zero divisor, and the Fano/MTC channel map. The canonical Web then carries an explicit depth-3 `e↔q` loop ledger, and every depth-4 backreaction names and actually reads that parent.

## Mixing, couplings, gravity, baryogenesis, dark sector

| Sector | Result | Reference | Status |
|:---|:---|:---|:---|
| Electroweak scale `v_EW` | 246.219645 GeV | 246.219651 GeV (from G_F) | −0.02 ppm |
| Fine structure `1/α(0)` | 137.035999050 | Berkeley Cs 137.035999046(27) | +0.13σ vs Cs, **−14σ vs Rb 2020**. The commitment is to the Cs side of the 5.5σ dispute. The dispute resolving toward Rb kills it |
| Higgs transport | `lambda(m_t) = 0.12651` | Internal transport | Candidate bridge boundary in the imported beta coordinate. The depth-4 Web reading is unassigned |
| CKM (15 comparison coordinates, correlations unmodeled) | diagonal pull score/coordinate = 0.21 | PDG 2024 | max marginal pull 0.7σ. The reference block mixes correlated global-fit outputs with a separately quoted phase constraint. Without a joint covariance model this is a descriptive score, not χ². Where mixing lives is itself an executable theorem: two `Z₃`-equivariant transports commute, giving `V = 1` exactly, so mixing is not a sector transport. It enters as four readings of the solved `D⁽⁶⁾*` boundary, one per face, rigid engine constants with no continuous freedom ([`mixing_z3_dichotomy.py`](tests/probes/mixing_z3_dichotomy.py)). Each reading carries its own grade, from computed (`η̄ = π/9`, the braiding monodromy) through enumerated (`λ`, the unique 1σ hit over the engine's echo menu, frozen in CI) to identified (the attachment to the physical Wolfenstein coordinates). The per-face grading is in [`mixing.py`](interference/mixing.py) |
| PMNS (3 correlated angles) | diagonal pull score/angle = 0.00 | Frozen NuFit 6.0 NO snapshot, IC19 without SK atmospheric data | Identified conjugation-orbit TBM map with `U_e = R₂₃(−φ)R₁₃(θ_C,δ)`. The construction consumes no new engine numbers: the TBM block is computed from the conjugation-orbit-compressed fusion matrices, and the corrections are functions of the derived mass spectrum, `φ = √(m_e/m_μ)/√d₁₀`, and the echoed `λ`. Symmetric marginal errors and no covariance make this descriptive rather than χ². The current `−arg(U_e3) = 77.0°` is convention-dependent and is not reported as the physical `δ_CP` |
| Strong coupling `α_s(M_Z)` | 0.1184 | PDG 0.1180(9) | +0.4σ (π/32 exact at μ* = 253.5 GeV) |
| Weinberg angle `sin²θ_W` | 0.231285 | PDG 0.23129(4) | pull −0.11σ |
| W mass `M_W` | 80.356 GeV | 80.3692(133) world avg | pull −1.00σ |
| Z mass `M_Z` | 91.189 GeV | 91.188 GeV | +0.0010% |
| Newton's constant `G_ind/G_N` | 0.999999917 (apportioned) | 1 | Internal closure under the Planck/healing-length identification: `Σ_target = 2πd₁₁ = 6π`. The probe's two factorizations are the same cutoff identity. Phase readings: 1.0000 (UV), 1.0134 (broken) |
| Baryon asymmetry `η_B` | 6.178e-10 | Planck 6.12e-10 | +0.9% |
| Dark / baryon ratio `Ω_DM/Ω_b` | `2π − 1` = 5.2832 | Planck 5.364 | pull −1.26σ |

**Joint cosmology closure.** The two cosmology channels are derived independently (`η_B` from the G₂ instanton, `Ω_DM/Ω_b` from bridge venting) but they must agree jointly with the absolute dark-matter density. They do. `η_B = 6.178e-10 → Ω_b h² = 0.02256` (+1.2σ vs Planck), and `Ω_DM h² = (2π−1) · Ω_b h² = 0.1192` (−0.7σ vs Planck 0.1200(12)). Two independent channels, one consistent cosmology. The BBN conversion η₁₀ = 273.9 Ω_b h² is cited, not derived: it is a change of coordinate on a printed number, and no engine step consumes it.

**Conditional structural constraint and cosmology bookkeeping.** If the two F₄ singlets are right-handed neutrinos and the mass mechanism is a minimal type-I seesaw with no other light-neutrino mass operators, then the light mass matrix has rank at most two. Under that assumption at least one state is massless. The ordering remains unselected (`m₁ = 0` for normal or `m₃ = 0` for inverted). Separately, `(3/8π) M_Pl² H₀²` is the Friedmann critical density, not `ρ_Λ`. The current dark-energy density is `Ω_Λ` times it. Jacobson's thermodynamic derivation leaves Λ as an integration constant, and CKN saturation is an identification hypothesis.

## What gravity is

Mass, above, was a clock. Gravity is the price of the ticking, seen from outside. In a standing wave the substrate circulates in place. That circulation takes up room in the shared background and depletes it, and the knot and its depletion form together (`tests/probes/skyrmion_3d.py`). **The soliton casts a shadow.** A depletion is a deficit, and a deficit is a pull. Every standing wave pushes on the same shared background, and their depletions compose: a planet's worth of circulations makes one giant shadow, a unified ventilation that all the particles in the planet breathe through. Having mass and sourcing gravity are one event.

**Gravity is not a fifth force layered on top of matter. It is the shape matter carves into the background while ventilating itself, and the bias other matter rides into.** The full treatment is in [`gravity.py`](interference/gravity.py) and [`interference/__init__.py`](interference/__init__.py).

The normalisation is audited from both directions. The minimal matter spectrum alone gives the wrong sign, `Σ_min ≈ −10.2`, repulsive gravity. The 182 bridge channels, a count forced by the wiring winner and entering linearly with no adjustable multiplicity, supply the `+29.126` that closes `G_ind/G_N` to `8×10⁻⁸`. The unfactored cutoff schemes miss on their own: `t₀ = ξ₀²` gives `G_ind/G_N = 1/12`, `2ξ₀²` gives `1/6`, `ℓ_Pl²` gives `1/3`, so the target `Σ = 2πd₁₁` sits on no unfactored scheme. The closure is internal, and its one measured contact is the derived `M_Pl` at +22 ppm ([`sakharov_normalisation.py`](tests/probes/sakharov_normalisation.py)).

## What closes

A substrate self-interferes, and many patterns form. Most do not close: lanes that make electron-like things whose ratios fail to continue, quark-like things that never bind. Six consistency gates eliminate them. Only one lane closes all six: protected `G₂` forgetting applied exactly once. The gates and their scan are the searchlight's bookkeeping, not a mechanism. The substrate closes or fails to close, and the scan is the computed record that nothing else passes, run on Lie-table invariants alone with its one genuine contender, `E₈ ⊃ E₇ × A₁`, printed beside the winner ([`embedding_uniqueness.py`](interference/searchlight/embedding_uniqueness.py)). The seven-branch gauntlet is the computed record of the shorter trees. The 504 hits reported later are mass grammar realizations within this fixed lane, not additional closure lanes.

That surviving lane sits inside one algebra:

```
14  +  52  +  182  =  248  =  dim(E₈)          E₈(1) ⊃ G₂(1) × F₄(1)
gauge   matter   bridge                         c: 14/5 + 26/5 = 8, exact
        + Higgs  → gravity                      248 → (14,1) ⊕ (1,52) ⊕ (7,26)
```

Every degree of freedom of `E₈` has a physical job. Nothing is left over. The audit is in [`interference/__init__.py`](interference/__init__.py).

## Quantum numbers

The `E₈ → G₂ × F₄` representation theory derives the Standard Model quantum numbers directly ([`root.py`](interference/root.py)):

```
Q_u  = d₁₀/d₁₁ = 2/3         (up-type electric charge)
Q_d  = 1/d₁₁   = 1/3         (down-type electric charge)
N_c  = d₁₁     = 3           (colour multiplicity)
sin²θ_W = d₁₁/(d₁₀²+d₁₁²) = 3/13   (Weinberg angle, at the embedding scale)
```

The charge trace `Σ Q² N_c = (d₁₀²+1+d₁₁)/d₁₁ = 8/3` gives the Singh ratio `(8/3)·C₂(26) = 16` that converts the `G₂` coupling to the electromagnetic one: `α_EM = α_G₂/16 = π/512`.

**CP violation** traces to the oriented Fano plane: seven directed lines select one orientation sign. Reversing it flips `η̄`, `J_CKM`, `J_lep`, and `η_B` together while leaving CP-even quantities unchanged ([`tests/probes/orientation_bit.py`](tests/probes/orientation_bit.py)). One global CP-orientation bit, not several independently inserted signs.

**Strong CP**: the index-one `G₂ → SU(3)_c` embedding preserves winding. In the no-Lagrangian ontology, the physical angle is the complete matter-coupled character `θ̄`, and the full record rule makes that character act on the vacuum record algebra `A = 1 + B`. Multiplication leaves only the trivial phase automorphism, so `θ̄ = 0 mod 2π`. This conclusion follows within that ontology rather than from local MTC data alone. A measured neutron EDM attributable to `θ̄ ≠ 0` kills it ([derivation](interference/couplings.py), [supporting probe](tests/probes/strong_cp.py)).

## Derivation flow

The scorecard is an executable dependency pipeline, run in order by [`interference/run.py`](interference/run.py):

1. [`wiring.py`](interference/wiring.py) scans the finite candidate menus (lattices, conformal pairs, levels) and finds the unique closure. The depth-3 searchlight reads `d₁₀ = 2`, `d₁₁ = 3`, `n₇ = 7`, `n₂₆ = 26`, and `α_G₂(M_Pl) = 1/(24π)` off the winner. [`root.py`](interference/root.py) consumes that vocabulary, and `run.py` asserts that both stages share the same wiring record.
2. [`root.py`](interference/root.py) builds the recursive web, its couplings and scales, and derives `M_Pl` by inverting the electron anchor.
3. [`words.py`](interference/words.py) derives the generation walk counts. [`masses.py`](interference/masses.py) derives the nine fermion masses and `v_EW`.
4. [`mixing.py`](interference/mixing.py) computes the WZW data, exact D⁽⁶⁾* connection, and mixing scorecards. [`couplings.py`](interference/couplings.py) computes `α_s` and `α(0)`.
5. [`gravity.py`](interference/gravity.py) computes the induced-gravity, baryogenesis, and neutrino outputs. [`higgs_ew.py`](interference/higgs_ew.py) computes the Higgs quartic boundary and electroweak observables (SM RGE transport, M_W, M_Z). [`dark_sector.py`](interference/dark_sector.py) computes the cosmological ratios.

The substrate is the physical realization of this chain. Its probes test emergence. They do not supply scorecard inputs.

## Architecture

The derivation lives in the [`interference/`](interference/) package as one dependency chain. Package-relative imports produce one canonical module graph and one Web singleton. The scrutiny lives in [`tests/`](tests/): freezes pin the numbers, probes cover mechanisms and falsifiers, and the substrate and audit instruments provide numerical checks. Each derivation module exposes one `derive(...) → dict`, and [`run.py`](interference/run.py) calls them in dependency order. The full narrative follows that order in [`interference/__init__.py`](interference/__init__.py).

| Module | Role |
|:---|:---|
| `run.py` | End-to-end runner: `python3 interference/run.py` |
| `root.py` | Embedding, Casimirs, central charges, the recursive echo web, `M_Pl` from the `m_e` anchor. The four integers are the depth-3 searchlight's vocabulary, read from `wiring.py` |
| `masses.py` | All nine masses: Brannen Z₃ leptons, F₄ quarks, `v_EW` |
| `mixing.py` | WZW/D⁽⁶⁾ boundary data, braiding phase, CKM parameter map, PMNS conjugation-orbit TBM construction |
| `couplings.py` | `α_s` (WZW cancellation at μ*), `α(0)` from bridge self-interference |
| `gravity.py` | Induced-gravity ledger, target factorizations, `η_B`, neutrinos |
| `higgs_ew.py` | SM RGE transport, Higgs quartic boundary, M_W, M_Z, import audit |
| `dark_sector.py` | `Ω_DM/Ω_b = 2π − 1` from bridge venting |
| `words.py` | Generation word lemma: nimrep walk counts → quark mass bases |
| `wiring.py` | Closure-residual scan over the candidate menus. The four integers are never an input |
| `reference.py` | Comparison data, single home, quarantined by `tests/test_blindness.py` |
| `searchlight/` | The octonions: `G₂/SU(3)` CG verification, embedding uniqueness through all six gates, protected forgetting. The one subdirectory, and the engine never imports it |
| `tests/substrate/` | Z₃-NLS substrate: BdG, stability, Madelung sourcing, and the same equation in 3D |
| `tests/audit/` | The committed grammar, and the grammar run as a generator: blind gates first, data contact once |


## Physics FAQ

The derivation probes in [`tests/probes/`](tests/probes/) answer standard physics objections. Each states its own verdict and kill condition. The table below is the discoverability map.

| Question | Answer | Code | Status |
|:---|:---|:---|:---|
| Why 3+1 dimensions? | The `Z₃` relative sector is `C²`, so the order parameter space is `S³`. `π₃(S³) = Z` first permits winding in 3 spatial dimensions. Topology-changing zeros of a four-real-component field have codimension 4, making them point-events exactly in `D = 4` | [`dimensionality.py`](tests/probes/dimensionality.py) | Derived |
| Why are solitons stable? | The bare NLS unwinds (intentional negative result in [`knot_charge.py`](tests/probes/knot_charge.py)). The Faddeev-Skyrme quartic stabilising term is then constructively derived from the framework's own energy functional | [`effective_action.py`](tests/probes/effective_action.py) | Derived |
| Is there confinement? | Vortex proliferation is measured, the area-law fit beats the perimeter-law fit decisively, and the Gaussian excitation spectrum is gapped | [`confined_phase.py`](tests/probes/confined_phase.py) | Classical confinement signature and Gaussian gap computed |
| How does gravity enter? | The exact Z₃ split isolates the common mode. Its sourced response and long-wavelength dispersion are computed. Given the declared common-mode/acoustic-geometry and Planck/healing-length identifications, the standard acoustic-metric and DeWitt–Seeley heat-kernel theorems produce the induced Einstein–Hilbert term, whose coefficient is evaluated by the venting ledger | [`gravity.py`](interference/gravity.py), [`sakharov_normalisation.py`](tests/probes/sakharov_normalisation.py) | Complete derivation within the stated identifications |
| What about chirality/CP? | Closure leaves the mass-protected chiral branch, dimensionality gives 3+1D, and the embedding branch supplies the gauge representations. The radial Skyrmion–Dirac operator is built on the repository's actual `π₃` hedgehog and exhibits the single localized crossing associated with unit winding. Its scanned coupling is a spectral-family coordinate, not a fitted or additional framework parameter. Independently, reversing the Fano orientation flips every computed CP-odd sign and preserves CP-even magnitudes. These branches compose under the stated fermion/knot and representation identifications | [`zero_mode.py`](tests/probes/zero_mode.py), [`orientation_bit.py`](tests/probes/orientation_bit.py), [`dimensionality.py`](tests/probes/dimensionality.py) | Complete derivation within the stated identifications. Knot spectral flow computed |
| Are EW quantities secretly fitted? | Import audit: top contribution from the framework's own `G_F`, `m_t`. Hadronic VP classified as measurement data, higher-order pieces as loop mathematics. No free theory parameters enter | [`ew_internal.py`](tests/probes/ew_internal.py) | Audited |
| Is the algebraic wiring arbitrary? | `(7,26)` bridge is the unique nontrivial boson of `G₂(1) × F₄(1)`. Face channels confined. Sector graph and cycle inventory derived from the modular tensor category | [`wiring_theorem_probe.py`](tests/probes/wiring_theorem_probe.py) | Derived |
| Can comparison data leak into predictions? | `test_blindness.py` perturbs reference values and requires predictions to remain identical. AST audit of `reference.py` access patterns | [`test_blindness.py`](tests/test_blindness.py) | Enforced |
| Does the knot form in 3D? | Relative-mode RMS grows ×15.6 under `G₂` coupling, flat under `g₁ = 0` control. Common-mode shadow forms at the core | [`skyrmion_3d.py`](tests/probes/skyrmion_3d.py) | Derived |
| What does the classical substrate claim, exactly? | The layer scoping: the classical field supplies the slots (a stable loop, three bound modes, channel lifting). The quantized bookkeeping supplies the ratios. A classical measurement of `Q = 2/3` falsifies the scoping | [`circulant_normalization.py`](tests/probes/circulant_normalization.py) | Committed with falsifier |

## Where the questions end

The engine takes one number, and even that number is a unit name rather than physics: the ratio `m_e/M_Pl` is itself derived, so the dimensionless output takes no input at all. How strongly the substrate interferes is not a second input either. It is the thing being derived, one channel at a time. What the axiom cannot reach is why the amplitude is nonzero at all. Self-presence is universal, not conditional: the empty world is present to itself too, and precisely because it is nothing the self-touching generates no amplitude. That is not an exception to the axiom but its trivial fixed point, the zero orbit of `x = b + W(x)`, zero in and zero out. Nothing does not escape the axiom. It just does not close down the road.


## Contributing

Verification, critique, tooling, and candidate derivations meeting the promotion bar (a committed edge menu and stated kill conditions) are welcome. Corrections with anything adjustable in them are not, because there is nothing here to tune.

## Papers

- The Innocent Lepton ([10.5281/zenodo.19899091](https://doi.org/10.5281/zenodo.19899091))
- One Substrate, Three Generations ([10.5281/zenodo.20069456](https://doi.org/10.5281/zenodo.20069456))
- The Echo of Standing Waves ([10.5281/zenodo.20144381](https://doi.org/10.5281/zenodo.20144381))
- The Octavian Coherence Gate: The Four Irreducible Integers of the E₈⁽¹⁾ ⊃ G₂⁽¹⁾ × F₄⁽¹⁾ Conformal Embedding ([10.5281/zenodo.20493955](https://doi.org/10.5281/zenodo.20493955))

Cite as: Kahsay, Kibrom Kidane (2026). One-field. https://github.com/kuwrom/one-field (the repository supersedes the papers as the living canon). See also [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).
