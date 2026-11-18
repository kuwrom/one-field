# one-field

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20323890.svg)](https://doi.org/10.5281/zenodo.20323890)
[![tests](https://github.com/kuwrom/one-field/actions/workflows/tests.yml/badge.svg)](https://github.com/kuwrom/one-field/actions/workflows/tests.yml)

Anything that exists already affects itself by being there. To be is, minimally, to be self-present. Interference is the physical image of that self-presence. A single photon passing through both slits interferes with itself. To be something is to carry non zero interference weight. The framework derives what emerges from that interference.

**A substrate self-interferes. Exactly one family of standing wave patterns closes to give our universe: the Standard Model, gravity, and the dark sector.**

Concretely, the code starts from one dimensional anchor (the electron mass) and the closure demand, whose unique solution becomes visible at depth three as one algebra. From these it produces the other eight fermion masses, the CKM and PMNS scorecards, the electroweak scale, the fine structure and strong couplings, Newton's constant, and the dark matter fraction. It has no free parameters and no structural choices. The framework has no Lagrangian. It derives each result by observing the interference, so we are calling every result an identification, read with the exhaustion of its committed menu and the reason only this lane closes.

The anchor sets the unit. Cited imports carry an engine reading to the coordinate a measurement is quoted in (`M_W` and `M_Z` in MS-bar, `α_s` at `M_Z`) ([`test_blindness.py`](tests/audit/test_blindness.py)). The runner prints the scorecard and the derivation graph. The substrate simulations test whether the structure emerges.

## Run it

```bash
git clone https://github.com/kuwrom/one-field.git
cd one-field
pip install -r requirements.txt
python3 interference/run.py   # full derivation, prints scorecard and derivation graph
pytest -q                     # 302 fast tests (pytest -m slow runs the 18 slow ones, 320 in all)
```

Python 3.10+, NumPy, SciPy. The tests add pytest, SymPy and python-flint (all in `requirements.txt`). The runner takes under a minute on one core, the fast tests a few minutes.

## Masses (one anchor in, eight out, max error 0.4%)

| Particle | Predicted | PDG | Error | Method | Status |
|---:|---:|---:|---:|:---|:---|
| e | 0.5109990 MeV | 0.5109990 MeV | anchor | CODATA, sole dimensional input | input |
| μ | 105.6584 MeV | 105.6584 MeV | +0.0000003% | Brannen Z₃ (predicted) | identified |
| τ | 1776.909 MeV | 1776.930 MeV | −0.0012% | Brannen Z₃ (predicted) | identified |
| u | 2.158 MeV | 2.16 MeV | −0.1% | (38/9) m_e | identified, conv |
| d | 4.713 MeV | 4.70 MeV | +0.3% | (83/9) m_e | identified, conv |
| s | 93.8 MeV | 93.5 MeV | +0.4% | Q₀+h₁₀/K³ + bridge | identified, conv |
| c | 1273.8 MeV | 1273 MeV | +0.06% | (217/18) m_μ | identified, conv |
| b | 4193.8 MeV | 4183 MeV | +0.3% | Q = Q₀ + h₁₁/K³ | identified, conv |
| t | 172.51 GeV | 172.57 GeV | −0.036% | (1165/12) m_τ | identified, conv |

<details>
<summary><b>Fine print: mass coordinates and status labels</b></summary>

**Status labels.** identified = a mathematical object observed as a physical one, with the exhaustion of its committed menu and the reason no other lane closes, separately falsifiable. Its theorems are executable in this repo. conv = the stated comparison coordinate of the external reference. input = the anchor, which sets the unit.

The framework produces one completed reading per fermion, with every echo included. It does not run to a scale, so each comparison names the coordinate of the external reference. Leptons are compared at the propagator pole, and the top at the PDG direct mass reference. Confined heavy quarks (c, b) are compared at the self-scale `m(m)`. Light quarks (u, d, s) are compared at the PDG `MS-bar(2 GeV)` coordinate. The charm pole reference near 1.67 GeV is a different coordinate, not a test of the reading. The Koide relation is an on-shell statement: transported to `M_Z`, the measured lepton masses give 0.17%. The claim is always the completed reading at its stated coordinate, never a relation at every scale. The coordinate free content of the light sector is its RG invariant ratios: `m_u/m_d = 38/83` (−0.9σ vs PDG 0.473(17)), `m_s/m_ud = 27.318` (+0.2σ vs 27.30(8)), `Q_ellipse = 22.383` (+0.4σ dispersive, −0.2σ FLAG lattice). See `masses.py` for the full treatment.

</details>

## What closes

A substrate self-interferes, and many patterns form. Most do not close. Some lanes (candidate closure routes) make electron like objects whose ratios do not carry on to the next sector. Others make quark like objects that never bind. Six consistency gates eliminate them. Only one lane passes all six: protected `G₂` forgetting applied exactly once. The gates and their scan are the searchlight's bookkeeping, not a mechanism. The substrate either closes or does not. The scan is the computed record that no other lane passes. It runs on Lie table invariants alone, and it prints its one real contender, `E₈ ⊃ E₇ × A₁`, beside the winner ([`embedding_uniqueness.py`](interference/searchlight/embedding_uniqueness.py)).

That surviving lane becomes visible at depth three as one algebra:

```
14  +  52  +  182  =  248  =  dim(E₈)          E₈(1) ⊃ G₂(1) × F₄(1)
gauge   matter   bridge                         c: 14/5 + 26/5 = 8, exact
        + Higgs  → gravity                      248 → (14,1) ⊕ (1,52) ⊕ (7,26)
```

Every degree of freedom of `E₈` has a physical job. Nothing is left over. The audit is in [`interference/__init__.py`](interference/__init__.py).

The same rule works at every scale. Particles, atoms and orbits last because they closed ([`spinning_lump.py`](tests/probes/spinning_lump.py), [`dimension_by_closure.py`](tests/probes/dimension_by_closure.py)), and solar systems and galaxies are the structures that closed and held together. What we see is what closed.

## The searchlight

Self-referencing can run arbitrarily deep, but mathematics can only track it through three layers of compositional history. A real number records the amplitude of a wave. A complex number records the amplitude and its phase. A quaternion records the phase of a phase. An octonion adds one more layer: the phase of a phase of a phase. These are not what the wave "is made of". They are the accounting ledger for tracking how interference patterns reference their own history. The octonion algebra is the deepest ledger that stays consistent: a fourth Cayley–Dickson doubling would introduce zero divisors and the bookkeeping would contradict itself (Hurwitz's theorem). [`searchlight/octonions.py`](interference/searchlight/octonions.py) verifies bilinearity, norm multiplicativity, the explicit sedenion zero divisor, and the Fano/MTC channel map.

The searchlight is the octonionic algebra used to locate the structure. The octonions reduced blind interference fluctuations to seven directions. Through the unique `(7,26)` bridge they reduced the search to a finite 7×26 one. Within a sector the first closure is the three cycle `a → b → c → a`. That bound is the territory's. The map's bound is Hurwitz: three is the deepest the octonionic searchlight resolves. The two bounds are independent, and that they meet at three is why the closing loop could be found at all. The searchlight limits what can be tracked, not what exists. The engine never imports it.

Searchlight depth is the octonionic tracking resolution (at most three). The octonionic searchlight resolves the closed `a → b → c → a` loop at depth 3. The scalar recursive web can then refine it at depth 4 and beyond through registered derived edges. Deeper self-interference is carried by the recursive web, not by extending the Cayley–Dickson tower. The canonical Web carries an explicit depth 3 `e↔q` loop ledger. Every depth 4 backreaction names that ledger as its parent and reads it.

The constants are stored nowhere. They are the fixed point the recursion reaches, as a Chladni figure is the pattern a plate settles into. The Chladni plate analogy holds as a method: the pattern first becomes visible at finite resolution, and later recursive readings refine it. The recursive web is open to later resolution but closed against hindsight. Every new depth is registered with its derivation before comparison.

## What mass is

A particle is a closed standing wave in the field, whose phase circulates at one fixed frequency. It needs nothing around it to hold it together: it is a self-trapped lump ([`lump_in_empty_field.py`](tests/probes/lump_in_empty_field.py)). That frequency is the mass, `E = ħω`. A heavier particle is "heavy" only in that its internal clock is faster.

The mass construction has three clock settings because it uses three child labels with `Z₃` cyclic symmetry (`a → b → c → a`). Three is not a choice. The closing lane's `G₂ ⊃ SU(3)` supplies a `Z₃` centre, and a single sector changing insertion first returns to its sector at third order (`tr S = tr S² = 0`, `tr S³ = 3`, [`masses.py`](interference/masses.py)), so within a sector the first closure is the three cycle. The Web's self loop and two cycle join the G₂ and F₄ sectors. When the bias that singles out one child label is erased (protected forgetting), the resulting `Z₃`-circulant has three eigenstates identified with the electron, muon, and tau generations. The circulant fixes the mass ratios. The electron fixes their dimensional unit. `m_μ` and `m_τ` are then determined by `m_e` alone.

**Lepton mass chain.** The calculation uses Brannen's 2006 circulant form. The `G₂ → SU(3)` Fano plane calculation gives the Clebsch-Gordan ratio `√2`. The closing lane's `SU(3)₃` modular data have canonical positive spherical dimension `d(1,0) = 2`. The channel reading map combines these data with two cycle orientations and gives `B/A = √2`. The circulant identity then gives `Q₀ = 2/3`. The lightest triality changing primary has Sugawara weight `h = 2/9`, and the mass phase map sets `θ = h`. This is the only positive, nondegenerate value in the committed angle menu. With `m_e` as the dimensional anchor, the calculation reproduces `m_μ` to 0.0000003% and `m_τ` to 0.0012%.

Any configuration that closes to give the leptons also produces quarks and gravity, from the same structure. The universe does not have to get lucky twice. The lepton standing wave closes without involving the embedding and sits at the confinement scale. The quark standing wave closes through the embedding and sits at the electroweak scale. The three generations of each come from the same `Z₃`, not from triality. The mass sector is one generator plus a grammar: the lepton circulant is the generator, and the quark masses are arithmetic on it. The full mass chain is in [`interference/__init__.py`](interference/__init__.py).

## Quantum numbers

The Standard Model quantum numbers are read off the four integers (`d₁₀ = 2`, `d₁₁ = 3`, `n₇ = 7`, `n₂₆ = 26`) by one dictionary, an identification stated in [`root.py`](interference/root.py):

```
Q_u  = d₁₀/d₁₁ = 2/3         (up-type electric charge)
Q_d  = 1/d₁₁   = 1/3         (down-type electric charge)
N_c  = d₁₁     = 3           (colour multiplicity)
sin²θ_W = d₁₁/(d₁₀²+d₁₁²) = 3/13   (Weinberg angle, zero edge base)
```

The charge trace `Σ Q² N_c = (d₁₀²+1+d₁₁)/d₁₁ = 8/3` gives the Singh ratio `(8/3)·C₂(26) = 16`. This ratio converts the `G₂` coupling at the end of its lever arm, `α_G₂(μ*) = π/32`, to the zero edge electromagnetic base (the value before any echo edge): `α_EM = α_G₂/16 = π/512`.

`α_EM = π/512` and `α(0) = 1/137.035999050` are two coordinates of one coupling. `α(0)` is the zero momentum reading: the bridge self loop gives the factor `(1 − 1/(2π))`, and the depth 4 recursive reading of the electron↔quark carrier (two orientations) adds `−2(α/2π)²`, so `1/α(0)` is the real root of `x³ = (512/π)[(1 − 1/(2π))x² − 1/(2π²)] = 137.035999050`. The depth 1 truncation alone, `(2⁹/π)(1 − 1/(2π)) = 137.036439`, is an intermediate identity, not the prediction.

**CP violation** is one global orientation bit. It becomes visible at depth three as the oriented Fano plane: its seven directed lines select one orientation sign. Under the map from Fano orientation to phase sign, reversing the orientation flips `η̄`, `J_CKM` and `J_lep` together, and leaves CP-even quantities unchanged ([`tests/probes/orientation_bit.py`](tests/probes/orientation_bit.py)).

**Strong CP**: the index one `G₂ → SU(3)_c` embedding preserves winding. In the framework's ontology, which has no Lagrangian, the physical angle is the complete matter coupled character `θ̄`. The full record rule makes that character act on the vacuum record algebra `A = 1 + B`. Multiplication leaves only the trivial phase automorphism, so `θ̄ = 0 mod 2π`. The rule is the records primitive ontology's own. A measured neutron EDM attributable to `θ̄ ≠ 0` kills it ([derivation](interference/couplings.py), [supporting probe](tests/probes/strong_cp.py)).

## Mixing, couplings, gravity, dark sector

| Sector | Result | Reference | Status |
|:---|:---|:---|:---|
| Electroweak scale `v_EW` | 246.219645 GeV | 246.219640 GeV (from G_F) | +0.02 ppm |
| Fine structure `1/α(0)` | 137.035999050 | Berkeley Cs 137.035999046(27) | +0.13σ vs Cs, −14σ vs Rb. The framework commits to the Cs side of the 5.5σ dispute. A resolution toward Rb kills it |
| CKM (15 comparison coordinates) | diagonal pull score/coordinate = 0.21 | PDG 2024 | max marginal pull 0.7σ. Without a joint covariance model this is a descriptive score, not χ². Where mixing lives is itself an executable theorem. Two `Z₃`-equivariant transports commute and give `V = 1` exactly, so mixing is not a sector transport. It enters as four readings of the solved `D⁽⁶⁾*` boundary, one per face. These are fixed engine constants with no continuous freedom ([`mixing_z3_dichotomy.py`](tests/probes/mixing_z3_dichotomy.py)). The per face grading is in [`mixing.py`](interference/mixing.py) |
| PMNS (3 correlated angles) | `sin²θ₁₂ = 0.30783`, `sin²θ₂₃ = 0.56139`, `sin²θ₁₃ = 0.02190`. Diagonal pull score/angle = 0.00 | Frozen NuFit 6.0 NO snapshot, IC19 without SK atmospheric data | Identified conjugation orbit TBM map with `U_e = R₂₃(−φ)R₁₃(θ_C,δ)`. The construction uses no new engine numbers. Symmetric marginal errors and no covariance make this descriptive rather than χ². The rephasing invariant `δ_CP = 275.74°` lies inside all NuFIT 6.0 3σ ranges. Two data blind gates fix its sign: the lepton sector's own pairing and the charge dictionary. The NuFIT ranges confirm it ([`pmns_sign_bits.py`](tests/probes/pmns_sign_bits.py)) |
| Strong coupling `α_s(M_Z)` | 0.1187 | PDG 0.1180(9) | +0.8σ (π/32 exact at μ* = 253.5 GeV) |
| Weinberg angle `sin²θ_W` | 0.231285 | PDG 0.23129(4) | pull −0.11σ |
| W mass `M_W` | 80.356 GeV | 80.3625(77) world avg, without CDF II | pull −0.86σ |
| Z mass `M_Z` | 91.189 GeV | 91.188 GeV | +0.0010% |
| Newton's constant `G` | 6.673936×10⁻¹¹ | CODATA 6.67430(15)×10⁻¹¹ | −2.43σ under CODATA's multiplicative combination, within 1σ under an added spread or the unweighted median ([`g_measurements.py`](tests/probes/g_measurements.py) G6). `M_Pl` comes from the anchor inversion `m_e → M_Pl` (+22 ppm). `G` reads it through the common mode (the field component that carries no label): `G M_b² = 3Q_l − 1`, the lepton circulant's Koide moment at depth four. One of nine readings of this count survives ([`common_mode_count.py`](tests/probes/common_mode_count.py)). Identified |
| Dark / baryon ratio `Ω_DM/Ω_b` | `2π − 1` = 5.2832 | Planck 5.364 | pull −1.26σ (ACT DR6 +1.24σ) |

**Neutrino rank and the vacuum term.** The searchlight forces the shape of the neutrino sector at depth three: Dirac masses, the two F₄ singlets as the right handed halves, the lightest neutrino massless, and no neutrinoless double beta decay ([`neutrino_mass_term.py`](tests/probes/neutrino_mass_term.py)). An observed neutrinoless double beta decay or a third nonzero mass kills it. Either ordering is allowed (`m₁ = 0` for normal, `m₃ = 0` for inverted). Separately, `(3/8π) M_Pl² H₀²` is the Friedmann critical density, not `ρ_Λ`, and the dark energy density today is `Ω_Λ` times it. The vacuum term has one closing form, the traceless one, with Λ an integration constant ([`lambda_lane.py`](tests/probes/lambda_lane.py)). Cohen-Kaplan-Nelson (CKN) saturation is an outside bound (12.2 times the observed density), not a derivation of it.

## What gravity is

Mass, above, is a clock. Gravity is the effect of that ticking on the shared background, seen from outside.

A standing wave touches the background at its edge, and the background starts moving too. The background is not something else. It is the same field, the part all standing waves share, called the common mode. Every standing wave moves it, and they move it together, which means less friction than moving it alone. Other matter follows that shared motion. Having mass and sourcing gravity are the same event.

**Gravity is not a fifth force added on top of matter. It is the motion that matter gives the background, and that other matter follows.** The full treatment is in [`gravity.py`](interference/gravity.py) and [`interference/__init__.py`](interference/__init__.py).

The strength of gravity is the anchor inversion read through the common mode. The gravity lane has six gates: sign, universality, magnitude, one cone, tensor carrier, and no preferred location. Each names the tests that derive it and the case in which it fails ([`gravity.py`](interference/gravity.py), [`test_gravity_gates.py`](tests/common_mode/test_gravity_gates.py)).

## Where the questions end

The engine takes one number, and even that number is a choice of unit, not physics. The ratio `m_e/M_Pl` is itself derived, so the dimensionless output takes no input at all. How strongly the substrate interferes is not a second input either. It is what is being derived, one channel at a time. What the axiom cannot explain is why the amplitude is nonzero at all. Self-presence is universal, not conditional: the empty world is present to itself too, and precisely because it is nothing the self-touching generates no amplitude. That is not an exception to the axiom but its trivial fixed point: the zero orbit of `x = b + W(x)`, zero in and zero out. Nothing does not escape the axiom. It simply never closes into anything further.

## Architecture

The derivation lives in the [`interference/`](interference/) package as one dependency chain. Package relative imports give one canonical module graph and one Web singleton. The checks live in [`tests/`](tests/), one folder per branch (foundations, masses, couplings, mixing, common mode, substrate, audit). Freezes pin the numbers, probes cover mechanisms and falsifiers, and the substrate and audit instruments provide numerical checks. Each derivation module exposes one `derive(...) → dict`, and [`run.py`](interference/run.py) calls them in dependency order. The full narrative follows that order in [`interference/__init__.py`](interference/__init__.py).

| Module | Role |
|:---|:---|
| `root.py` | Embedding, Casimirs, central charges, the recursive echo web, `M_Pl` from the `m_e` anchor. The four integers are the depth 3 searchlight's vocabulary, read from `wiring.py` |
| `masses.py` | All nine masses: Brannen Z₃ leptons, F₄ quarks, `v_EW` |
| `mixing.py` | WZW/D⁽⁶⁾ boundary data, braiding phase, CKM parameter map, PMNS conjugation orbit TBM construction |
| `couplings.py` | `α_s` (WZW cancellation at μ*), `α(0)` from bridge self-interference |
| `gravity.py` | Induced gravity ledger, the six gravity gates, neutrinos |
| `higgs_ew.py` | SM RGE transport, Higgs quartic boundary, M_W, M_Z, import audit |
| `dark_sector.py` | `Ω_DM/Ω_b = 2π − 1` from bridge venting |
| `wiring.py` | Closure residual scan over the candidate menus. The four integers are never an input |
| `searchlight/` | The octonions: `G₂/SU(3)` CG verification, embedding uniqueness through six coherence conditions on parent/factor pairs, protected forgetting. The one subdirectory, and the engine never imports it |

## Contributing

Verification, critique, tooling, and candidate derivations meeting the promotion bar (a committed edge menu and stated kill conditions) are welcome. Corrections with anything adjustable in them are not, because there is nothing here to tune.

## Paper

- The Octavian Coherence Gate: The Four Irreducible Integers of the E₈⁽¹⁾ ⊃ G₂⁽¹⁾ × F₄⁽¹⁾ Conformal Embedding ([10.5281/zenodo.20493955](https://doi.org/10.5281/zenodo.20493955)). The four integers and the uniqueness of the embedding. The chain from the integers to the scorecard is this repository.

Related records: The Innocent Lepton ([10.5281/zenodo.19899091](https://doi.org/10.5281/zenodo.19899091)), One Substrate, Three Generations ([10.5281/zenodo.20069456](https://doi.org/10.5281/zenodo.20069456)), The Echo of Standing Waves ([10.5281/zenodo.20144381](https://doi.org/10.5281/zenodo.20144381)).

Cite as: Kahsay, Kibrom Kidane (2026). One-field. https://github.com/kuwrom/one-field, with the Octavian Coherence Gate paper. See also [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).
