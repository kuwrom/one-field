# one-field

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20323890.svg)](https://doi.org/10.5281/zenodo.20323890)
[![tests](https://github.com/kuwrom/one-field/actions/workflows/tests.yml/badge.svg)](https://github.com/kuwrom/one-field/actions/workflows/tests.yml)

To be is, minimally, to be self-present. Anything that exists already affects itself by being there. Interference is the physical image of this. What is said to exist is what carries non-zero interference weight. The framework observes what emerges from interference.

**A substrate self-interferes. Exactly one family of standing-wave patterns closes to give our universe: the Standard Model, gravity, and the dark sector.**

Concretely, from one dimensional anchor (the electron mass) and the closure condition that forces the algebra, the code in this repository derives the other eight fermion masses, the CKM and PMNS matrices, the electroweak scale, the fine-structure and strong couplings, Newton's constant, the baryon asymmetry, and the dark-matter fraction, with zero dimensionless free parameters. Every predicted number in this README is printed by the code, and every link states its own status (derived or identified) in the module that computes it. The fastest way to read it is to run it first.

## Run it

```bash
git clone https://github.com/kuwrom/one-field.git
cd one-field
pip install -r requirements.txt
python3 interference/run.py   # full derivation, prints scorecard and derivation graph
pytest -q                     # 169 tests (154 fast + 15 slow behind -m slow)
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

Unconfined fermions (leptons, top) sit at the propagator pole. Confined heavy quarks (c, b) sit at the self-scale `m(m)`. Light quarks (u, d, s) sit at the PDG `MS-bar(2 GeV)` coordinate. The scheme-free content of the light sector is the RG-invariant ratios, which carry no coordinate at all: `m_u/m_d = 38/83` (−0.9σ vs PDG 0.473(17)), `m_s/m_ud = 27.318` (+0.2σ vs 27.30(8)), `Q_ellipse = 22.383` (+0.4σ dispersive, −1.7σ lattice, and the two references disagree). See `masses.py` for the full treatment.

</details>

**The precision hierarchy is structural.** The circulant form is Brannen's (2006). The parameters `B/A = √2` and `θ = 2/9` are derived here as algebraic outputs of the `G₂ → SU(3)` Clebsch-Gordan data. Leptons sit closest to the `Z₃` source and emerge cleanly (`m_μ`, `m_τ` predicted from `m_e` to 0.0012%). Quarks pass through additional layers (triality, generation scaling, confinement), and each layer adds residuals of order `α_s/π` (≤ 0.4%). Gravity, two layers further out, lands at 0.0-1.3%. The accuracy gradient is not noise. It is what the emergence depth predicts.

## What mass is

A particle is a self-trapped patch of the substrate, a lump whose internal phase circulates at one fixed frequency. That rate is the mass, `E = ħω`. A heavier particle is a faster internal clock, and nothing else about it is "heavy".

The clock comes in exactly three settings because an interference loop needs at least three nodes, and closure was found at the minimum. The three components couple in a cycle (`a → b → c → a`), and when the bias that singles out one position is erased (protected forgetting), the resulting `Z₃`-circulant has three eigenstates: the electron, muon, and tau generations. The circulant fixes the mass ratios. The electron fixes their dimensional unit. `m_μ` and `m_τ` are then determined by `m_e` alone.

Any configuration that closes to give the leptons automatically produces quarks and gravity from the same structure. The universe does not have to get lucky twice. The lepton standing wave closes without engaging the embedding and sits at the confinement scale. The quark standing wave closes through the embedding and sits at the electroweak scale. Three generations of each are the same `Z₃`, not triality. The mass sector is one generator plus a grammar: the lepton circulant is the generator, and the quark masses are arithmetic on it. The full mass chain is in [`interference/__init__.py`](interference/__init__.py).

## Why octonions

The smallest self-interfering loop needs three nodes. The deepest consistent phase-tracking algebra (the octonions, by Hurwitz's theorem) resolves exactly three nestings. Because those two independent threes coincide, the octonions reduced the search from blind fluctuations to seven directions and a finite menu, and luckily the closing structure was found there. The engine's arithmetic never evaluates an octonion. Once found, it stands on its own. [`searchlight/octonions.py`](interference/searchlight/octonions.py) serves as an independent cross-check, not a dependency.

## Mixing, couplings, gravity, baryogenesis, dark sector

| Sector | Result | Reference | Status |
|:---|:---|:---|:---|
| Electroweak scale `v_EW` | 246.219645 GeV | 246.219651 GeV (from G_F) | −0.02 ppm |
| Fine structure `1/α(0)` | 137.035999050 | Berkeley Cs 137.035999046(27) | +0.13σ vs Cs, **−14σ vs Rb 2020**. The commitment is to the Cs side of the 5.5σ dispute; the dispute resolving toward Rb kills it |
| Higgs transport | `lambda(m_t) = 0.12651` | — | Candidate bridge boundary in the imported beta coordinate; the depth-4 Web reading is not identified here |
| CKM (15 observables) | χ²/n = 0.21 | PDG 2024 | max pull 0.7σ. The 9 \|V_ij\| have parameter rank 4, so the 15 entries are 4 independent measurements against 4 engine numbers, not 15. Of the four: η̄ = π/9 is derived (a braiding phase); λ, A, ρ̄ are identifications, and the mixing mechanism itself is a no-go, executed twice (the label-algebra dichotomy and the six-quark crossing lift): the web computes frame-blind data, and a finite menu is all such data can carry |
| PMNS (3 angles) | χ²/n = 0.00 | NuFit 6.0 | 3 parameters → 3 angles (zero dof: the χ² carries no information; the content is δ_CP). Predicts δ_CP = 77.0°, **in present tension with the NuFit 6.0 NO global fit, which centres near 177° and disfavours the first quadrant**. DUNE/Hyper-K decide (kill condition) |
| Strong coupling `α_s(M_Z)` | 0.1184 | PDG 0.1180(9) | +0.4σ (π/32 exact at μ* = 253.5 GeV) |
| Weinberg angle `sin²θ_W` | 0.231285 | PDG 0.23129(4) | pull −0.11σ |
| W mass `M_W` | 80.356 GeV | 80.3692(133) world avg | pull −1.00σ |
| Z mass `M_Z` | 91.189 GeV | 91.188 GeV | +0.0010% |
| Newton's constant `G_ind/G_N` | 0.999999917 (face-split) | 1 | 8.3e-08 internal, phase readings 1.0000 (UV) - 1.0134 (broken) |
| Baryon asymmetry `η_B` | 6.177e-10 | Planck 6.12e-10 | +0.9% |
| Dark / baryon ratio `Ω_DM/Ω_b` | `2π − 1` = 5.2832 | Planck 5.364 | pull −1.26σ |

**Joint cosmology closure.** The two cosmology channels are derived independently (`η_B` from the G₂ instanton, `Ω_DM/Ω_b` from bridge venting) but they must agree jointly with the absolute dark-matter density. They do. `η_B = 6.177e-10 → Ω_b h² = 0.02255` (+1.2σ vs Planck), and `Ω_DM h² = (2π−1) · Ω_b h² = 0.1192` (−0.7σ vs Planck 0.1200(12)). Two independent channels, one consistent cosmology. The BBN conversion η₁₀ = 273.9 Ω_b h² is cited, not derived: it is a change of coordinate on a printed number, and no engine step consumes it.

**Structural predictions.** `m_ν₁ = 0` (rank-2 seesaw, two RH neutrinos in F₄), normal ordering, and the cosmological-constant scale `ρ_Λ ~ M_Pl² H₀²` (the 10¹²³ problem structurally resolved).

## What gravity is

Mass, above, was a clock. Gravity is the price of the ticking, seen from outside. In a standing wave the substrate circulates in place. That circulation takes up room in the shared background and depletes it, and the knot and its depletion form together (`tests/probes/skyrmion_3d.py`). **The soliton casts a shadow.** A depletion is a deficit, and a deficit is a pull. Every standing wave pushes on the same shared background, and their depletions compose: a planet's worth of circulations makes one giant shadow, a unified ventilation that all the particles in the planet breathe through. Having mass and sourcing gravity are one event.

**Gravity is not a fifth force layered on top of matter. It is the shape matter carves into the background while ventilating itself, and the bias other matter rides into.** The full treatment is in [`gravity.py`](interference/gravity.py) and [`interference/__init__.py`](interference/__init__.py).

## What closes

A substrate self-interferes, and many patterns form. Most do not close: lanes that make electron-like things whose ratios fail to continue, quark-like things that never bind. Six consistency gates eliminate them. Only one lane closes all six: protected `G₂` forgetting applied exactly once. The seven-branch gauntlet is the computed record of the shorter trees.

That surviving lane sits inside one algebra:

```
14  +  52  +  182  =  248  =  dim(E₈)          E₈(1) ⊃ G₂(1) × F₄(1)
gauge   matter   bridge                         c: 14/5 + 26/5 = 8, exact
        + Higgs  → gravity                      248 → (14,1) ⊕ (1,52) ⊕ (7,26)
```

Every degree of freedom of `E₈` has a physical job. Nothing is left over. The full gate-by-gate audit is in [`interference/__init__.py`](interference/__init__.py).

## Quantum numbers

The `E₈ → G₂ × F₄` representation theory derives the Standard Model quantum numbers directly ([`root.py`](interference/root.py)):

```
Q_u  = d₁₀/d₁₁ = 2/3         (up-type electric charge)
Q_d  = 1/d₁₁   = 1/3         (down-type electric charge)
N_c  = d₁₁     = 3           (colour multiplicity)
sin²θ_W = d₁₁/(d₁₀²+d₁₁²) = 3/13   (Weinberg angle, at the embedding scale)
```

The charge trace `Σ Q² N_c = (d₁₀²+1+d₁₁)/d₁₁ = 8/3` gives the Singh ratio `(8/3)·C₂(26) = 16` that converts the `G₂` coupling to the electromagnetic one: `α_EM = α_G₂/16 = π/512`.

**CP violation** traces to the oriented Fano plane: seven directed lines select a preferred chirality. Reversing the orientation flips `η̄`, `J_CKM`, `J_lep`, and `η_B` together while leaving CP-even quantities unchanged ([`tests/probes/orientation_bit.py`](tests/probes/orientation_bit.py)). One global orientation bit, not several independently inserted signs.

**Strong CP**: `π₃(G₂) = Z` and the `G₂ → SU(3)_c` embedding has index 1, so instanton number is preserved. `E₈(1)` has a unique integrable vacuum module, hence `θ_G₂ = 0` and `θ_QCD = 0` exactly. No axion, no new symmetry. A measured neutron EDM attributable to `θ_QCD ≠ 0` kills it ([`couplings.py`](interference/couplings.py)).

## Derivation flow

The scorecard is an executable dependency pipeline, run in order by [`interference/run.py`](interference/run.py):

1. [`wiring.py`](interference/wiring.py) scans the finite candidate menus (lattices, conformal pairs, levels) and finds the unique closure. The depth-3 searchlight reads `d₁₀ = 2`, `d₁₁ = 3`, `n₇ = 7`, `n₂₆ = 26`, and `α_G₂(M_Pl) = 1/(24π)` off the winner. [`root.py`](interference/root.py) consumes that vocabulary, and `run.py` asserts that both stages share the same wiring record.
2. [`root.py`](interference/root.py) builds the recursive web, its couplings and scales, and derives `M_Pl` by inverting the electron anchor.
3. [`words.py`](interference/words.py) derives the generation walk counts. [`masses.py`](interference/masses.py) derives the nine fermion masses and `v_EW`.
4. [`mixing.py`](interference/mixing.py) computes the WZW data, exact D⁽⁶⁾* connection, and mixing scorecards. [`couplings.py`](interference/couplings.py) computes `α_s` and `α(0)`.
5. [`gravity.py`](interference/gravity.py) computes the induced-gravity, baryogenesis, and neutrino outputs; [`higgs_ew.py`](interference/higgs_ew.py) computes the Higgs quartic boundary and electroweak observables (SM RGE transport, M_W, M_Z). [`dark_sector.py`](interference/dark_sector.py) computes the cosmological ratios.

The substrate is the physical realization of this chain. Its probes test emergence. They do not supply scorecard inputs.

## Architecture

The derivation lives in [`interference/`](interference/) and is flat, because it reads and runs as one chain rather than as a package. The scrutiny lives in [`tests/`](tests/): checking whether the theory is flawless is testing, so the freezes pin the numbers, the probes in `tests/probes/` cover the mechanisms and the rabbit holes, and the substrate and the audit instruments sit there too, on the scrutiny side. Each derivation module exposes one `derive(...) → dict`, and [`run.py`](interference/run.py) calls them in dependency order. The full narrative follows that order in [`interference/__init__.py`](interference/__init__.py).

| Module | Role |
|:---|:---|
| `run.py` | End-to-end runner: `python3 interference/run.py` |
| `root.py` | Embedding, Casimirs, central charges, the recursive echo web, `M_Pl` from the `m_e` anchor. The four integers are the depth-3 searchlight's vocabulary, read from `wiring.py` |
| `masses.py` | All nine masses: Brannen Z₃ leptons, F₄ quarks, `v_EW` |
| `mixing.py` | WZW/D⁽⁶⁾ boundary data, CKM parameter map, PMNS |
| `couplings.py` | `α_s` (WZW cancellation at μ*), `α(0)` from bridge self-interference |
| `gravity.py` | Conditional induced-gravity closure, `η_B`, neutrinos |
| `higgs_ew.py` | SM RGE transport, Higgs quartic boundary, M_W, M_Z, import audit |
| `dark_sector.py` | `Ω_DM/Ω_b = 2π − 1` from bridge venting |
| `words.py` | Generation word lemma: nimrep walk counts → quark mass bases |
| `wiring.py` | Closure-residual scan over the candidate menus; the four integers are never an input |
| `reference.py` | Comparison data, single home, quarantined by `tests/test_blindness.py` |
| `searchlight/` | The octonions: `G₂/SU(3)` CG cross-check, embedding uniqueness through all six gates, protected forgetting. The one subdirectory, and the engine never imports it |
| `tests/substrate/` | Z₃-NLS substrate: BdG, stability, Madelung sourcing, and the same equation in 3D |
| `tests/audit/` | The committed grammar, and the grammar run as a generator: blind gates first, data contact once |


## Physics FAQ

The derivation probes in [`tests/probes/`](tests/probes/) answer standard physics objections. Each states its own verdict and kill condition. The table below is the discoverability map.

| Question | Answer | Code | Status |
|:---|:---|:---|:---|
| Why 3+1 dimensions? | The `Z₃` relative sector is `C²`, so the order parameter space is `S³`. `π₃(S³) = Z` first permits winding in 3 spatial dimensions. Topology-changing zeros of a four-real-component field have codimension 4, making them point-events exactly in `D = 4` | [`dimensionality.py`](tests/probes/dimensionality.py) | Derived |
| Why are solitons stable? | The bare NLS unwinds (intentional negative result in [`knot_charge.py`](tests/probes/knot_charge.py)). The Faddeev-Skyrme quartic stabilising term is then constructively derived from the framework's own energy functional | [`effective_action.py`](tests/probes/effective_action.py) | Derived |
| Is there confinement? | Vortex proliferation measured. Area-law fit beats perimeter-law decisively. Quantum topological order (gap + degeneracy + box scaling) is still open | [`confined_phase.py`](tests/probes/confined_phase.py) | Classical signature derived |
| Does gravity produce GR-like physics? | Common acoustic metric is the induced spin-2 channel. `c_gw/c = 1` identically, PPN `γ = 1` exactly, equivalence principle universal (`PvP = 0` algebraic identity). Scalar q=0 channel Yukawa-screened at the healing length | [`gravity.py`](interference/gravity.py) | Derived |
| What about chirality/CP? | One Fano orientation bit propagates: CP-odd quantities flip together, CP-even quantities unchanged | [`orientation_bit.py`](tests/probes/orientation_bit.py) | Derived (V-A vertex open) |
| Are EW quantities secretly fitted? | Import audit: top contribution from the framework's own `G_F`, `m_t`. Hadronic VP classified as measurement data, higher-order pieces as loop mathematics. No free theory parameters enter | [`ew_internal.py`](tests/probes/ew_internal.py) | Audited |
| Is the algebraic wiring arbitrary? | `(7,26)` bridge is the unique nontrivial boson of `G₂(1) × F₄(1)`. Face channels confined. Sector graph and cycle inventory derived from the modular tensor category | [`wiring_theorem_probe.py`](tests/probes/wiring_theorem_probe.py) | Derived |
| Can comparison data leak into predictions? | `test_blindness.py` perturbs reference values and requires predictions to remain bit-identical. AST audit of `reference.py` access patterns | [`test_blindness.py`](tests/test_blindness.py) | Enforced |
| Does the knot form in 3D? | Relative-mode RMS grows ×15.6 under `G₂` coupling, flat under `g₁ = 0` control. Common-mode shadow forms at the core | [`skyrmion_3d.py`](tests/probes/skyrmion_3d.py) | Derived |

## Where the questions end

The engine takes one number, and even that number is a unit name rather than physics: the ratio `m_e/M_Pl` is itself derived, so the dimensionless output takes no input at all. How strongly the substrate interferes is not a second input either. It is the thing being derived, one channel at a time. What the axiom cannot reach is why the amplitude is nonzero at all. Self-presence is universal, not conditional: the empty world is present to itself too, and precisely because it is nothing the self-touching generates no amplitude. That is not an exception to the axiom but its trivial fixed point, the zero orbit of `x = b + W(x)`, zero in and zero out. Nothing does not escape the axiom. It just does not close down the road.

## What is derived vs identified

Every quantity in this framework carries an explicit epistemic status. The distinction matters: the claim is not that everything follows from pure mathematics, but that the status of every link is stated, and the identification hypotheses are enumerable.

- **Thm**: the discrete structure is an executable theorem in this repo. The six gates, the wiring scan, the circulant algebra, the conformal weights, the walk counts, the orientation bit.
- **IH**: an identification hypothesis is load-bearing. A mathematical object is matched with a physical one. Example: identifying the `Z₃` centre that crystallises when `G₂ → SU(3)` with three particle generations, or carrying the amplitude ratio `B/A = √2` from the derived octonionic CG into the substrate coupling. Each IH has a stated consumer and is separately falsifiable.
- **conv**: a stated scheme convention enters. Comparing a derived mass reading with PDG requires naming a coordinate (pole, MS-bar, self-scale). The convention is declared, not absorbed.
- **obs**: an arithmetic observation whose derivation is not claimed. `Q₀ = d₁₀/d₁₁` is an identity at the closing values; the derivation takes the form `Q₀ = 1/3 + d₁₀/6` instead.

The look-elsewhere denominator is printed, not asserted: [`tests/audit/closure_search.py`](tests/audit/closure_search.py) runs the committed grammar as a generator and reports 504 of 2,370,816 blind-viable webs hitting the 1% scorecard band. The word counts are resolved uniquely; the sub-percent corrections are not. Canon is one of 504 at that band, not uniquely selected there.

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
