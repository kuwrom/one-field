"""
Masses: from (d₁₀, d₁₁, n₇, n₂₆, π) + the electron anchor.
m_e is the sole dimensional input (M_Pl is derived in root.py).
m_μ, m_τ and the six quark masses are predictions.

All derivations are self-contained, every formula traces to the four
irreducible numbers via the WZW/MTC structures below.

Gate 4 (closure) lives here.  The Brannen circulant Δ_k = A + B cos(θ +
2πk/3) is the UNIQUE form compatible with the cyclic constraint
[Δ, S] = 0, where S is the Z₃ shift on the three child-centre labels.
The direct finite-G₂ path retains the harmonic v = (1, −1/2, −1/2) and
does NOT close.  The cyclic branch erases v and closes on circulant
algebra alone, this is the survivor filter that selects the unique
lane.  See __init__.py for the full six-gate chain.

═══════════════════════════════════════════════════════════════════════
WHAT MASS IS
═══════════════════════════════════════════════════════════════════════

The substrate self-interferes.  Against vanishing odds, interference
loops emerge: a -> b -> c -> a.  A loop carries cyclic Z₃ symmetry,
because nothing inside it distinguishes which position came first.
Most loops just cycle.

Some, though, drop the bias that singles out one of their three
positions.  The asymmetric label v = (1, −1/2, −1/2) inherited from
the parent G₂ structure is erased, and all three positions become
equivalent.  A forgetting loop can no longer drift between its
positions.  It gets stuck in one of three eigenstates.  That stuck
pattern is a standing wave, and the three eigenstates are the three
lepton generations: electron, muon, tau.  Same loop, three places it
can lock in.

Mass is that pattern.

═══════════════════════════════════════════════════════════════════════
THE SUBSTRATE PICTURE (lepton paper, self-contained summary)
═══════════════════════════════════════════════════════════════════════

The effective local model is the Z₃-symmetric coupled nonlinear
Schrödinger system for three child condensates ψ_k:

    iħ ∂ψ_k = −D∇²ψ_k + Σ_j g_kj |ψ_j|² ψ_k,   k = 1,2,3
    g_kk = g₀,  g_kj = g₁ (k≠j),  with the G₂ CONSTRAINT g₁ = g₀/√d₁₀.

A lepton is a Z₃-equivariant Bogoliubov-de Gennes (BdG) standing
wave on this system.  Closure imposes cyclicity [Δ, S] = 0 on the
BdG mass-amplitude matrix (the "forgetting" that erases which of
the three phases the loop came from). The surviving spectrum is the
circulant gap Δₖ below, and the mass is the gap squared: mₖ = Δₖ².

BdG/KOIDE IDENTITY (theorem): for any real Z₃-equivariant positive
triplet Δₖ = A + B cos(θ + 2πk/3) with mₖ = Δₖ²,
    Q ≡ Σmₖ/(Σ√mₖ)² = 1/3 + B²/(6A²)  exactly.
On the positive branch Q = 2/3 ⟺ |B/A| = √2.  So the Koide
relation is an algebraic identity of the branch.

═══════════════════════════════════════════════════════════════════════
LEPTON MASSES, Z₃ circulant → Koide → masses
═══════════════════════════════════════════════════════════════════════

Step 1: Z₃ circulant structure (derives the Koide relation).
    The SU(3)₃ MTC has a Z₃ center symmetry (simple currents).
    Any hermitian matrix commuting with the Z₃ shift S is a circulant:
        Δ = c₀I + c₁S + c₂S†
    Reality forces c₂ = c₁*, giving three real eigenvalues:
        Δₖ = A + |B| cos(θ + 2πk/3),   k = 0,1,2
    On the positive-gap branch (all Δₖ > 0), the BdG standing-wave
    picture gives √mₖ = Δₖ, so:
        √mₖ = A[1 + (|B|/A) cos(θ + 2πk/3)]
    Using Σ cos(θ+2πk/3) = 0 and Σ cos²(θ+2πk/3) = 3/2:
        Q ≡ Σmₖ / (Σ√mₖ)² = 1/3 + |B/A|²/6
    This is a mathematical identity of ANY Z₃-equivariant spectrum.
    Credit: Koide (1982) for the empirical observation. Brannen (2006)
    for the parameterisation.  The derivation is from Z₃ circulant algebra.

Step 2: B/A = √d₁₀ (derives the Koide value Q₀ = 2/3).
    Two independent structures fix B/A:
    (a) Octonionic CG coefficients from G₂ = Aut(O):
        G₂ breaks to SU(3).  The 7 → 3 ⊕ 3̄ ⊕ 1.
        Octonionic cross product gives:
          singlet channel:  |C₁| = 1
          triplet channel:  |C₃̄| = √2  (from Fano-plane structure)
        The ratio |C₃̄|/|C₁| = √2 = √d₁₀ is convention-independent
        (Schur's lemma).
    (b) Quantum Schur suppression from SU(3)₃ MTC:
        The off-diagonal (charged) fusion channel propagates through
        the fundamental representation with quantum dimension d(3) = d₁₀ = 2.
        The quantum Schur lemma gives a 1/d(fund) = 1/d₁₀ suppression.
    Assembly: B/A = 2|C₃̄| / (|C₁| · d₁₀) = 2√2/2 = √2 = √d₁₀
    Then Q₀ = 1/3 + d₁₀/6 = 1/3 + 1/3 = 2/3 = d₁₀/d₁₁.

Step 3: θ = h₁₀ = d₁₀/d₁₁² = 2/9 (derives the phase).
    The SU(3)₃ Sugawara construction gives conformal weight:
        h(R) = C₂(R) / (k + h∨)
    For the fundamental: C₂(fund, SU(d₁₁)) = (d₁₁²−1)/(2d₁₁) = 4/3
    At level k = d₁₁ = 3, altitude K = k + h∨ = d₁₀d₁₁ = 6:
        h₁₀ = (4/3) / 6 = 2/9 = d₁₀/d₁₁²
    The mass insertion is a local two-point OPE coupling, so the
    effective potential V_eff(θ) = −Σ d(λ)² cos(q(θ − h(λ))) has its
    minimum at θ = h(fund) = 2/9 (the lightest charged primary dominates).

Step 4: Confinement scale A = ½ M_Pl exp(−S₀), S₀ = d₁₁²π²/2,
    with M_Pl DERIVED from the electron anchor (root.py inversion).
    The absolute scale chain (lepton paper): α_G₂(M_Pl) = |Z₃|/(2π
    D²_tot) = 1/(24π) from the SU(3)₃ modular tensor category
    (executable in root.py), trace-lift normalisation c_eff = 1/2.
    S₀ = h∨(F₄)·π²/2 is the F₄ instanton action.  The quark layer
    then vents the action to S = 9π²/2 − 6 + 15/512 (Casimir vent
    −C₂(26), 30-mode vertex echo +15/512. root.py ledger), giving
    v_EW = M_Pl·exp(−S), the same exponent in both readings.

Step 5: Pole mass = tree × (1 − α(0)/(2π)).
    QED vertex back-reaction at the mass shell (standard QFT calls this
    the "vertex correction". Here it is the EM channel's interference
    with the propagator pole, forced by the field content).

═══════════════════════════════════════════════════════════════════════
UP-TYPE QUARKS, WZW emergence
═══════════════════════════════════════════════════════════════════════

    m_q(n) = [base(n) + corr(n)] × m_lep(n)

    base ladder:  d₁₀²  →  d₁₀²d₁₁  →  d₁₁⁴+d₁₀⁴
    corr ladder:  h₁₀   →  δ          →  1/(2K)
                  2/9       1/18          1/12
    Ratios: corr₁/corr₂ = d₁₀² = 4,  corr₂/corr₃ = Q₀ = 2/3  (alternating)

═══════════════════════════════════════════════════════════════════════
DOWN-TYPE QUARKS, Koide with WZW representation back-reactions
═══════════════════════════════════════════════════════════════════════

The Koide relation Q = 2/3 is exact for leptons (identity rep, h=0).
For quark triplets, each closes through a specific SU(3)₃ representation λ,
producing a sub-leading back-reaction:
    Q = Q₀ + h(λ)/K³

FIRST-INVARIANT-ORDER THEOREM (why K³): the insertion is the Z₃
sector-changer, and the confined sector admits only Z₃-invariant
operators ([Δ,S] = 0, selector theorem).  One or two insertions do
not return to the identity sector (tr S = tr S² = 0). Three do
(tr S³ = 3).  Orders 1 and 2 are forbidden, not small. The first
admissible order is cubic.  In OPE language: Σ√m receives its first
Z₃-invariant back-reaction at third order in the conformal weight
expansion.  The same theorem gives bridge² its K² altitude: a paired
S S† insertion is invariant already at order 2.  (Witness:
tests/test_coverage.py. Derivation program 1, registry.)
ENUMERATION over all eight structures Q₀ + h(λ)/Kᵖ (λ ∈ {fund, adj},
p ∈ {1..4}): only the adjoint weight at cubic altitude closes on the
bottom mass (+0.26%). The nearest alternative misses by +2.16%, all
others by >5%.  The (s,c,b) assignment (h₁₀, K³) carries the same
representation content at the same depth.

    (c,b,t) triplet: adjoint (1,1), h₁₁ = 1/d₁₀ = 1/2
        → Q = d₁₀/d₁₁ + (1/d₁₀)/(d₁₀d₁₁)³ = 289/432
    (s,c,b) triplet: fundamental (1,0), h₁₀ = d₁₀/d₁₁² = 2/9
        → Q = d₁₀/d₁₁ + (d₁₀/d₁₁²)/(d₁₀d₁₁)³ = 649/972

    bridge² = Q₀² · d₁₀³/d₁₁ = 32/27
    (Albert algebra trace norm × Dynkin Z₂ orientation factor)

═══════════════════════════════════════════════════════════════════════
LIGHT QUARKS, d₁₀² ↔ d₁₁² swap
═══════════════════════════════════════════════════════════════════════

    m_u/m_e = d₁₀² + h₁₀ = 38/9   (up-type: d₁₀² = quantum dim²)
    m_d/m_e = d₁₁² + h₁₀ = 83/9   (down-type: d₁₁² = quantum dim²)
    The swap d₁₀² ↔ d₁₁² distinguishes up from down via the triality
    of the SO(8) subalgebra: up quarks see the fundamental channel d₁₀²,
    down quarks see the adjoint channel d₁₁².

═══════════════════════════════════════════════════════════════════════
THE MASS COORDINATE, one rule over three dynamical classes
═══════════════════════════════════════════════════════════════════════

A running mass m(μ) is a coordinate on an RG orbit. A mass prediction
is complete only with its coordinate.  The table above is stated in
one rule:

  LEPTONS (e, μ, τ): the PROPAGATOR POLE, the only scheme-independent
    mass an asymptotic state has.  The QED mass-shell factor
    (1 − α(0)/2π) above IS the tree → pole conversion.

  TOP QUARK (t): also the PROPAGATOR POLE.  The top decays before
    hadronising (Γ_t ≈ 1.4 GeV >> Λ_QCD ≈ 0.2 GeV), so its pole
    mass is as well-defined as a lepton's: renormalon ambiguity is
    Λ_QCD/m_t ≈ 0.12%, safely below measurement error.  The entry
    matches the PDG kinematic mass.

  CONFINED HEAVY quarks (c, b): the SELF-SCALE m(m), the unique
    fixed point of μ ↦ m(μ).  No pole exists below confinement
    (renormalon ambiguity Λ_QCD/m ~ 4.8% for b, 15.7% for c), so
    the fixed-point coordinate is the canonical one, and it is the
    PDG reference coordinate.

  SCHEME DEPENDENCE NOTE: choosing the wrong coordinate inflates
    residuals catastrophically.  Example: the charm pole mass
    m_c^pole ≈ 1.67 GeV shifts the Koide residual from +0.06%
    (at m(m)) to roughly −24%.  Every mass in this table is stated
    in the coordinate the physics selects.

  LIGHT quarks (u, d, s): no perturbative self-scale exists
    (m < Λ_QCD).  The scheme-free content is the RG-INVARIANT
    RATIOS (QCD running is flavour-blind, so it cancels in m_q/m_q'):
        m_u/m_d = 38/83 = 0.45783        PDG 2024: 0.473(17)  (−0.9σ)
        m_s/m_ud = 27.318                PDG 2024: 27.30(8)   (+0.2σ)
        Q_ellipse = 22.383               η→3π dispersive: 22.1(7)
                                         (+0.4σ). Lattice 23.4(6)
                                         (−1.7σ. The two data
                                         determinations disagree,
                                         PDG review Sec. 60)
    The absolute entries are quoted in the PDG MS-bar(2 GeV)
    coordinate, a declared dictionary entry.

The algebra does not run to a scale. Its output IS the prediction.
A reader who wants a different convention applies standard RGE
transport with the framework's own α_s(M_Z) = 0.1184 (couplings.py)
-- both
endpoints are algebraic, so zero freedom enters -- but that is the
SM's coordinate change, not a step in the prediction.
"""

import math
from fractions import Fraction


def derive(R):
    """
    Derive all 9 fermion masses.

    R : dict from root.derive()
    """
    d10, d11 = R['d10'], R['d11']
    h10, h11, delta = R['h10'], R['h11'], R['delta']
    K, Q0 = R['K'], R['Q0']

    print("\n" + "=" * 78)
    print("  MASSES: 9 predictions from (d₁₀, d₁₁) = (2, 3)")
    print("=" * 78)

    # ── Leptons (Z₃ circulant → Koide → masses) ────────────────────
    #
    # √mₖ = A[1 + √d₁₀ cos(h₁₀ + 2πk/3)]
    #
    # B/A = √d₁₀ = √2:
    #   octonionic CG |C₃̄|/|C₁| = √2 from Fano cross product,
    #   × quantum Schur suppression 2/d(fund) = 2/d₁₀ = 1
    #   → B/A = 2√2/(1·2) = √2  (three independent proofs: surgery,
    #     correlator decomposition, partial quantum trace)
    #
    # θ = h₁₀ = d₁₀/d₁₁² = 2/9:
    #   Sugawara h(R) = C₂(R)/(k+h∨), fund: C₂ = (d₁₁²−1)/(2d₁₁) = 4/3
    #   → h₁₀ = (4/3)/(d₁₀d₁₁) = d₁₀/d₁₁² = 2/9
    #
    # A = ½ M_Pl exp(−S₀), S₀ = d₁₁²π²/2 (instanton action, see root.py)
    #
    # Q₀ = 1/3 + |B/A|²/6 = 1/3 + d₁₀/6 = 2/3 = d₁₀/d₁₁
    #   (Z₃ circulant identity)
    #
    from root import (M_Pl_MeV, BA_ratio, S_lepton, PDG_MASSES, pct,
                      WEB, Ledger)

    # Each lepton is born as a LEDGER: the knot's circulating amplitude
    # plus its echo stack.  Charged knots vent through the EM channel, so
    # they inherit the lepton_EM echo from the web, the old explicit
    # "× QED_factor" is gone. The echo law produces it.
    lepton_tree, lepton_pred = {}, {}
    for k_idx, name in [(1, 'e'), (2, 'mu'), (0, 'tau')]:
        gap = 1.0 + BA_ratio * math.cos(h10 + 2*math.pi*k_idx/3)
        # live ratio base: the SAME standing wave read from the solved
        # ruler node, so the whole mass sector is one coupled recursion
        led = Ledger(f"m_{name}",
                     (lambda s, g=gap:
                      0.5 * s["M_Pl_MeV"] * math.exp(-S_lepton) * g * g),
                     "mul", "MeV")
        for t in WEB["lepton_EM"].terms:
            led.echo(t.path, t.factor, t.depth, t.status, t.note, t.kind)
        # FORCED depth-3 vent (vertex composition rule): the mass node's
        # own channel weight W(fund) = 1/d₁₀ (the B/A Schur projector)
        # composes with the Albert traversal dim J₃(O) = 27 across the
        # e↔q loop → 27/2, venting ∝ the knot's own amplitude Δₖ.
        # The vertex rule gives the multiplicity exactly: 27/2.
        led.echo(["e↔q(J₃(O))"],
                 (lambda s, g=gap:
                  -(d11**3 / d10)
                  * ((1.0/s["inv_alpha"]) / (2*math.pi))**2 * g),
                 3, "FORCED",
                 "vertex rule: W(fund)·dim J₃(O) = 27/2, ∝ Δₖ")
        WEB[f"m_{name}"] = led
        lepton_tree[name] = 0.5 * M_Pl_MeV * math.exp(-S_lepton) * gap**2
        lepton_pred[name] = led.value()

    m_e, m_mu, m_tau = lepton_pred['e'], lepton_pred['mu'], lepton_pred['tau']

    # the electron is the ANCHOR (sole dimensional input): its ledger
    # value reproduces the input exactly by the M_Pl inversion in root
    from root import M_E_ANCHOR_MEV
    assert abs(lepton_pred['e']/M_E_ANCHOR_MEV - 1.0) < 1e-12

    print(f"\n  Leptons: √mₖ = A[1 + √{d10} cos({Fraction(d10, d11**2)} + 2πk/3)]")
    print(f"    m_e   = {lepton_pred['e']:12.7f} MeV  ANCHOR (sole "
          f"dimensional input, CODATA ±0.3 ppb)")
    for n in ['mu', 'tau']:
        print(f"    m_{n:<3s} = {lepton_pred[n]:12.7f} MeV  PREDICTED "
              f"(PDG: {PDG_MASSES[n]:12.6f},  {pct(lepton_pred[n], PDG_MASSES[n]):+.4f}%)")

    # ── Up-type quarks (WZW emergence) ───────────────────────────────
    #
    # Unified: m_q(n) = [base(n) + corr(n)] × m_lep(n)
    #
    #   base ladder:  d₁₀²  →  d₁₀²d₁₁  →  d₁₁⁴+d₁₀⁴
    #   corr ladder:  h₁₀   →  δ          →  1/(2K)
    #                 2/9      1/18          1/12
    #
    #   Ratios: corr₁/corr₂ = d₁₀² = 4,  corr₂/corr₃ = Q₀ = 2/3 (alternating)

    # ECHO-LAW FORM (exact rational identities, same numbers, one law):
    # each quark knot binds to its generation lepton and VENTS through
    # its closure channel. What standard QFT calls additive "corrections"
    # are the vents
    # read multiplicatively:
    #   m_u = d₁₀²·(1 + δ)    · m_e    [vent = OPE gap δ = 1/18, atom]
    #   m_c = d₁₀²d₁₁·(1 + K⁻³)· m_μ   [vent = 1/K³, SAME suppression
    #                                   as the down-type Koide terms]
    #   m_t = base₃·(1 + 1/(2K·base₃))·m_τ  [base₃ = d₁₁⁴+d₁₀⁴ = 97.
    #                                   vent rule pending base(n) task]
    # Proof of exactness: d₁₀²(1+δ) = 4·19/18 = 38/9.
    # d₁₀²d₁₁(1+1/216) = 12·217/216 = 217/18. 97(1+1/1164) = 1165/12.
    base_u, base_c, base_t = d10**2, d10**2 * d11, d11**4 + d10**4
    # live ratio bases: quark = word-count integer × the lepton node,
    # read from the state (the ratio edge of the one graph)
    led_u = Ledger("m_u", lambda s: base_u * s["m_e"], "mul", "MeV").echo(
        ["WZW(δ)"], delta, 1,
        "FORCED", "vent = OPE gap δ (dictionary atom): h₁₀/d₁₀² = δ")
    led_c = Ledger("m_c", lambda s: base_c * s["m_mu"], "mul", "MeV").echo(
        ["WZW(K³)"], 1.0 / K**3, 1,
        "FORCED", "vent = K⁻³ (same suppression as Koide Q-terms)")
    # base(n) = boundary-walk counts of generation words (Lemma
    # 'Generation words', words.py): w₁=ff (4),
    # w₂=ffa (12, step letter forced by the neutral lane C₀ = C_ad),
    # w₃=f⁴⊕a⁴ (97, Z₂-orbit terminus at J³=1).  The terminus is
    # DERIVED (knot definition + conversion lemma. Module-category
    # form verified in words.py). Mixed words (36, 78) are vent-side
    # and independently killed by experiment.
    led_t = Ledger("m_t", lambda s: base_t * s["m_tau"], "mul", "MeV").echo(
        ["WZW(2K)"], 1.0 / (2 * K * base_t), 1,
        "FORCED", "vent = corr₃/base₃. words lemma, terminus DERIVED")
    WEB["m_u"], WEB["m_c"], WEB["m_t"] = led_u, led_c, led_t

    m_u, m_c, m_t = led_u.value(), led_c.value(), led_t.value()

    coeff_u = Fraction(d10**2) * (1 + Fraction(1, d10*d11**2))
    coeff_c = Fraction(d10**2 * d11) * (1 + Fraction(1, (d10*d11)**3))
    coeff_t = Fraction(d11**4 + d10**4) + Fraction(1, 2*d10*d11)
    assert coeff_u == Fraction(38, 9) and coeff_c == Fraction(217, 18)

    print(f"\n  Up-type quarks (WZW emergence):")
    print(f"    m_u = ({coeff_u}) m_e   = {m_u:.4f} MeV  ({pct(m_u, PDG_MASSES['u']):+.1f}%)")
    print(f"    m_c = ({coeff_c}) m_μ  = {m_c:.1f} MeV  ({pct(m_c, PDG_MASSES['c']):+.2f}%)")
    print(f"    m_t = ({coeff_t}) m_τ = {m_t/1e3:.2f} GeV  ({pct(m_t, PDG_MASSES['t']):+.3f}%)")
    print(f"    Back-reaction ladder: h₁₀={Fraction(d10, d11**2)} →÷d₁₀²→ δ={Fraction(1, d10*d11**2)} →÷Q₀→ 1/2K={Fraction(1, 2*d10*d11)}")
    print(f"      ratios: d₁₀² = {d10**2},  Q₀ = {Fraction(d10, d11)}")

    # ── Down-type quarks (Koide with WZW representation back-reactions) ─
    #
    # The Koide relation Q = Q₀ = d₁₀/d₁₁ is exact for leptons (identity
    # rep, h = 0).  For quark triplets, each closes through a specific
    # SU(3)₃ WZW representation λ, producing a sub-leading back-reaction:
    #     Q = Q₀ + h(λ)/K³
    # (Standard QFT would call these "representation corrections". Here
    #  each h(λ)/K³ is a forced interference term.)
    #
    # K³ denominator: the altitude cubed arises because the Koide sum rule
    # involves Σ√m, and √m ∝ Δ (gap) receives an OPE back-reaction at third
    # order in the conformal weight expansion.  The leading and quadratic
    # terms vanish by Z₃ symmetry, leaving the cubic as the first non-trivial
    # contribution.
    #
    # (c,b,t): adjoint (1,1), h₁₁ = 1/d₁₀ = 1/2
    #   → Q = 2/3 + (1/2)/216 = 289/432
    # (s,c,b): fundamental (1,0), h₁₀ = d₁₀/d₁₁² = 2/9
    #   → Q = 2/3 + (2/9)/216 = 649/972
    #
    # Echo-law form: the Koide observable of each triplet is a ledger,
    # base Q₀ (the universal circulant value) plus the triplet's vent
    # through its closure representation, suppressed by K³ (the same K³
    # as the charm vent above: one depth law for up- and down-type).
    led_Qcbt = Ledger("Q(c,b,t)", Q0, "add").echo(
        ["WZW(adj)"], h11 / K**3, 1, "FORCED", "adjoint vent h₁₁/K³")
    led_Qscb = Ledger("Q(s,c,b)", Q0, "add").echo(
        ["WZW(fund)"], h10 / K**3, 1, "FORCED", "fundamental vent h₁₀/K³")
    WEB["Q_cbt"], WEB["Q_scb"] = led_Qcbt, led_Qscb
    Q_cbt = led_Qcbt.value()             # 2/3 + (1/d₁₀)/(d₁₀d₁₁)³ = 289/432
    Q_scb = led_Qscb.value()             # 2/3 + (d₁₀/d₁₁²)/(d₁₀d₁₁)³ = 649/972

    # Bottom from Q(c,b,t)
    sc, st = math.sqrt(m_c), math.sqrt(m_t)
    a = Q_cbt - 1
    b = 2*Q_cbt*(sc + st)
    c_coeff = (Q_cbt - 1)*(m_c + m_t) + 2*Q_cbt*sc*st
    disc = b**2 - 4*a*c_coeff
    x1 = (-b + math.sqrt(disc)) / (2*a)
    x2 = (-b - math.sqrt(disc)) / (2*a)
    m_b = x1**2 if m_c < x1**2 < m_t else x2**2

    # Strange from Q(s,c,b) + bridge
    sc_r, sb_r = math.sqrt(m_c), math.sqrt(m_b)
    a_r = 1 - Q_scb
    b_r = 2*Q_scb*(sc_r + sb_r)
    c_r = m_c + m_b - Q_scb*(sc_r + sb_r)**2
    disc_r = b_r**2 - 4*a_r*c_r
    x1_r = (-b_r + math.sqrt(disc_r)) / (2*a_r)
    x2_r = (-b_r - math.sqrt(disc_r)) / (2*a_r)
    x_phys = x1_r if x1_r > 0 else x2_r
    # ── COMPRESSION: the 8/3 trinity ─────────────────────────────────
    # One object, three appearances.  The dictionary identity
    #     d₁₀³ = d₁₁² − 1 = d₁₀² + d₁₁ + 1   (= 8)
    # (equivalent to the Sugawara lock d₁₁²−1 = 4d₁₀ at d₁₀ = 2) makes
    # d₁₀³/d₁₁ = (d₁₀²+1+d₁₁)/d₁₁ = CHARGE TRACE = 8/3 (root.py), so:
    #   • Singh ratio       = charge_trace × C₂(26) = 16
    #   • strange bridge²   = Q₀² × charge_trace × (1 + h₁₁/K²) = 32/27 × 73/72
    #   • top back-reaction  = h₁₀ / charge_trace    = corr₃ = 1/12
    from root import charge_trace, WEB
    assert d10**3 == d11**2 - 1 == d10**2 + d11 + 1
    assert Fraction(d10**3, d11) == charge_trace
    assert Fraction(2, 9) / charge_trace == Fraction(1, 12)  # corr₃
    bridge_sq = WEB.state["bridge_sq"]   # web ledger: 32/27 × (1 + h₁₁/K²)
    m_s = x_phys**2 * math.sqrt(bridge_sq)

    # Light quarks: d₁₀² ↔ d₁₁² swap gives up vs down
    #   m_u/m_e = d₁₀² + h₁₀ = 38/9   (up-type base)
    #   m_d/m_e = d₁₁² + h₁₀ = 83/9   (down-type base)
    # m_d = d₁₁²(1 + d₁₀/d₁₁⁴) m_e exactly: 9·(1+2/81) = 83/9.
    # Down-type light vent = d₁₀/d₁₁⁴ (the d₁₀²↔d₁₁² triality swap
    # acting on the same h₁₀ vent: h₁₀/d₁₁² = d₁₀/d₁₁⁴).
    led_d = Ledger("m_d", lambda s: d11**2 * s["m_e"], "mul", "MeV").echo(
        ["WZW(h₁₀)"], d10 / d11**4, 1,
        "FORCED", "vent = h₁₀/d₁₁² = d₁₀/d₁₁⁴ (triality-swapped)")
    WEB["m_d"] = led_d
    m_d = led_d.value()                  # (83/9) m_e

    # ── m_b, m_s as CONSTRAINT NODES of the web (live ratio bases) ───
    # The same Koide quadratics as the closed-form computation above,
    # read live from the state: the middle root between m_c and m_t,
    # and the positive root scaled by the bridge.  The kernel holds
    # them while the state is away from the fixed point. At the fixed
    # point they equal the closed forms exactly (asserted after the
    # final solve).
    def _m_bottom_of(s):
        sc_, st_ = math.sqrt(s["m_c"]), math.sqrt(s["m_t"])
        Q = s["Q_cbt"]
        a_ = Q - 1.0
        b_ = 2.0 * Q * (sc_ + st_)
        c_ = (Q - 1.0) * (s["m_c"] + s["m_t"]) + 2.0 * Q * sc_ * st_
        d_ = b_*b_ - 4.0*a_*c_
        if d_ < 0.0:
            raise ValueError("koide(b): no real root yet")
        for x in ((-b_ + math.sqrt(d_)) / (2.0*a_),
                  (-b_ - math.sqrt(d_)) / (2.0*a_)):
            if s["m_c"] < x*x < s["m_t"]:
                return x * x
        raise ValueError("koide(b): root not bracketed yet")

    def _m_strange_of(s):
        Q = s["Q_scb"]
        sc_, sb_ = math.sqrt(s["m_c"]), math.sqrt(s["m_b"])
        a_ = 1.0 - Q
        b_ = 2.0 * Q * (sc_ + sb_)
        c_ = s["m_c"] + s["m_b"] - Q * (sc_ + sb_)**2
        d_ = b_*b_ - 4.0*a_*c_
        if d_ < 0.0:
            raise ValueError("koide(s): no real root yet")
        x1 = (-b_ + math.sqrt(d_)) / (2.0*a_)
        x2 = (-b_ - math.sqrt(d_)) / (2.0*a_)
        x = x1 if x1 > 0 else x2
        return x * x * math.sqrt(s["bridge_sq"])

    WEB["m_b"] = Ledger("m_b", _m_bottom_of, "mul", "MeV")
    WEB["m_s"] = Ledger("m_s", _m_strange_of, "mul", "MeV")

    Q_cbt_frac = Fraction(d10, d11) + Fraction(1, d10) / Fraction(d10*d11)**3
    Q_scb_frac = Fraction(d10, d11) + Fraction(d10, d11**2) / Fraction(d10*d11)**3

    print(f"\n  Down-type quarks (Koide):")
    print(f"    Q(c,b,t) = Q₀ + h₁₁/K³ = {Q_cbt_frac} = {float(Q_cbt_frac):.10f}")
    print(f"    m_b = {m_b:.1f} MeV  ({pct(m_b, PDG_MASSES['b']):+.1f}%)")
    print(f"    Q(s,c,b) = Q₀ + h₁₀/K³ = {Q_scb_frac}")
    print(f"    bridge² = Q₀²·d₁₀³/d₁₁ = {Fraction(d10**2, d11**2) * Fraction(d10**3, d11)}")
    print(f"    m_s = {m_s:.1f} MeV  ({pct(m_s, PDG_MASSES['s']):+.1f}%)")
    print(f"    m_d = ({Fraction(d11**2) + Fraction(d10, d11**2)}) m_e = {m_d:.4f} MeV  ({pct(m_d, PDG_MASSES['d']):+.1f}%)")
    print(f"    Note: m_u/m_e = d₁₀²+h₁₀,  m_d/m_e = d₁₁²+h₁₀  (d₁₀²↔d₁₁² swap)")

    # ── Scheme-invariant light-quark ratios (no mass coordinate) ─────
    # QCD running is flavour-blind: m_q(μ)/m_q'(μ) is μ- and
    # scheme-independent.  These three numbers are predictions with
    # NO coordinate choice at all.
    mu_over_md = m_u / m_d                       # = 38/83 exactly
    assert abs(mu_over_md - 38.0/83.0) < 1e-12
    m_ud = 0.5 * (m_u + m_d)
    ms_over_mud = m_s / m_ud
    Q_ellipse = math.sqrt((m_s**2 - m_ud**2) / (m_d**2 - m_u**2))
    print(f"\n  Scheme-invariant ratios (RG-invariant, no coordinate):")
    print(f"    m_u/m_d  = 38/83 = {mu_over_md:.5f}   PDG 2024: 0.473(17)"
          f"   ({(mu_over_md-0.473)/0.017:+.1f}σ)")
    print(f"    m_s/m_ud = {ms_over_mud:.3f}            PDG 2024: 27.30(8)"
          f"    ({(ms_over_mud-27.30)/0.08:+.1f}σ)")
    print(f"    Q_ellipse = {Q_ellipse:.3f}           dispersive 22.1(7)"
          f" ({(Q_ellipse-22.1)/0.7:+.1f}σ), lattice 23.4(6)"
          f" ({(Q_ellipse-23.4)/0.6:+.1f}σ)")
    print(f"    (the two Q determinations disagree. PDG review Sec. 60)")

    # ── Summary ──────────────────────────────────────────────────────

    n1 = sum(1 for n in PDG_MASSES if abs(pct(
        {'e':m_e,'mu':m_mu,'tau':m_tau,'u':m_u,'d':m_d,'s':m_s,'c':m_c,'b':m_b,'t':m_t}[n],
        PDG_MASSES[n])) <= 1)
    max_err = max(abs(pct(
        {'e':m_e,'mu':m_mu,'tau':m_tau,'u':m_u,'d':m_d,'s':m_s,'c':m_c,'b':m_b,'t':m_t}[n],
        PDG_MASSES[n])) for n in PDG_MASSES)

    print(f"\n  Summary: {n1}/9 within 1%,  max error {max_err:.1f}%")

    # ── the weak scale as a graph node (live ratio base) ─────────────
    WEB["v_EW_GeV"] = Ledger(
        "v_EW_GeV",
        lambda s: s["M_Pl_MeV"] / 1e3 * math.exp(-s["S_quark"]),
        "mul", "GeV")

    WEB.solve()          # masses nodes join the solved web state
    v_EW_MeV = M_Pl_MeV * math.exp(-R['S_quark'])

    # the constraint nodes land exactly on the closed forms
    assert abs(WEB.state["m_b"] / m_b - 1.0) < 1e-12
    assert abs(WEB.state["m_s"] / m_s - 1.0) < 1e-12
    assert abs(WEB.state["v_EW_GeV"] * 1e3 / v_EW_MeV - 1.0) < 1e-12

    return {
        'lepton_pred': lepton_pred,
        'm_e': m_e, 'm_mu': m_mu, 'm_tau': m_tau,
        'm_u': m_u, 'm_d': m_d, 'm_s': m_s,
        'm_c': m_c, 'm_b': m_b, 'm_t': m_t,
        'v_EW_pred_MeV': v_EW_MeV,
        'v_EW_pred_GeV': v_EW_MeV / 1e3,
        'bridge': math.sqrt(bridge_sq),
        'h_fund': h10, 'delta_OPE': delta,
        'Q_val': Q_cbt, 'Q_s': Q_scb,
        'mu_over_md': mu_over_md, 'ms_over_mud': ms_over_mud,
        'Q_ellipse': Q_ellipse,
        'n_sub1': n1, 'max_err': max_err,
    }
