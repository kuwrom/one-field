"""
Layer 2 -- Lepton Masses from the G₂ sector.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

The Brannen Z₃ formula with B/A = √2, θ = 2/9:

    √mₖ = A [1 + √2 cos(2/9 + 2πk/3)]

yields all three charged-lepton pole masses from the
lepton confinement scale Λ_conf = ½ M_Pl exp(−9π²/2),
with the QED pole-mass factor (1 − α(0)/(2π)) where
α(0) = π²/(256(2π−1)) is the physical fine-structure constant,
an algebraic identity from the conformal embedding.

Physical motivation (readable without the papers)
──────────────────────────────────────────────────

WHY B/A = √2 (the amplitude ratio):

  In the Z₃ self-interfering Bose–Einstein condensate (the substrate),
  a standing-wave knot has three components with a gap vector
  v = (1, −½, −½).  The Z₃ cyclic projector averages over the three
  sectors.  The Bogoliubov-de Gennes (BdG) mass eigenvalues take
  the Brannen form √mₖ = A + B cos(θ + 2πk/3), and the Koide ratio
  Q = (m₁+m₂+m₃)/(√m₁+√m₂+√m₃)² = 1/3 + B²/(6A²).

  The octonionic Clebsch–Gordan calculation (see octonions.py) shows
  that the G₂/SU(3) branching has two channels:
    • singlet 3⊗3̄→1:  |C₁| = 1
    • antisymmetric 3⊗3→3̄:  |C₃̄| = √2
  The amplitude ratio B/A = 2|C₃̄|·W(3̄)/|C₁| = 2√2·(1/2)/1 = √2.
  Since |B/A|² = 2, the Koide identity gives Q₀ = 2/3 exactly.
  This is a property of the octonions, not a fit.

WHY θ = 2/9 (the Brannen phase):

  The phase θ is the conformal weight h(1,0) of the fundamental
  primary (1,0) in the SU(3)₃ WZW model:
    h(1,0) = C₂(fund, SU(3)) / (k + h∨)
           = (4/3) / (3 + 3) = 2/9
  where k=3 is the WZW level and h∨=3 the dual Coxeter number.
  This weight sets the scaling dimension of the lightest charged
  field in the SU(3)₃ boundary CFT, and therefore fixes the phase
  of the Z₃ gap triplet.  Both k=3 and h∨=3 are determined by the
  conformal embedding E₈(1) ⊃ G₂(1) × F₄(1).

WHY Λ_conf = ½ M_Pl exp(−9π²/2):

  The lepton confinement scale is set by the instanton action
  S₀ = 2π/(b₀·α_{G₂}) = 9π²/2, derived from one-loop dimensional
  transmutation in the G₂ sector with β-function coefficient
  b₀ = 32/3 and coupling α_{G₂}(M_Pl) = 1/(24π).  Substituting:
  S₀ = 2π × 24π / (32/3) = 2π × 24π × 3/32 = 9π²/2 ≈ 44.41.
  This is the scale where standing-wave knots first become
  energetically stable.  The factor ½ is the substrate healing
  length ξ₀ = ℓ_Pl/2.

WHY the QED correction factor (1 − α(0)/(2π)):

  α(0) = π²/(256(2π−1)) = 1/137.036 is an algebraic identity of
  the conformal embedding (derived in constants.py).  The factor
  (1 − α/(2π)) is the one-loop QED vertex correction that converts
  tree-level (bare) masses to pole masses.  This is not an external
  input: both the coupling and the correction form are determined
  by the algebraic structure.
"""

import math
from .constants import (
    M_Pl_MeV, QED_FACTOR, PDG_MASSES,
    THETA_BRANNEN, BA_RATIO,
)
from .formatting import H, S, pct


def derive(scale: dict):
    """
    Derive the three charged-lepton masses.

    Parameters
    ----------
    scale : dict from scale.derive()

    Returns
    -------
    dict with keys:
        lepton_tree   : {e, mu, tau} tree-level masses in MeV
        lepton_pred   : {e, mu, tau} pole masses in MeV
        m_e, m_mu, m_tau : shorthand MeV values
    """

    S_0 = scale['S_0']

    H("LAYER 2:  LEPTON MASSES -- G₂ SECTOR")

    S("2.1  Brannen formula with B/A = √2, θ = 2/9")

    # Gate 4 (closure) lives here.  The Brannen circulant Δ_k = A +
    # B cos(θ + 2 π k / 3) is the unique form compatible with the
    # cyclic constraint [Δ, S] = 0, where S is the Z_3 shift on the
    # three child-centre labels.  The direct finite-G_2 path retains
    # the harmonic v = (1, -1/2, -1/2) and does NOT close on this
    # algebra; the cyclic branch erases v and closes.  See __init__.py
    # for the full gate chain.
    #
    # The Brannen functional form (Brannen 2006) was originally an empirical
    # numerological ansatz.  In this framework B/A and theta are NOT fitted.
    # They are derived from the SU(3)_3 modular tensor category and the
    # octonionic G_2/SU(3) Clebsch-Gordan data:
    #
    #   theta = h(fund) = C_2(fund,SU(3))/(k+h^v) = (4/3)/(3+3) = 2/9
    #     (Sugawara: conformal weight of the lightest SU(3)_3 primary)
    #
    #   B/A = 2|C_3bar|/(|C_1|*d(fund)) = 2*sqrt(2)/(1*2) = sqrt(2)
    #     (octonionic Clebsch-Gordan: |C_3bar|=sqrt(2), |C_1|=1, d(fund)=2;
    #      see Innocent Lepton paper, Sec. 4)
    #
    # The Koide value Q=2/3 follows as a corollary of |B/A|=sqrt(2)
    # on the positive branch (BdG identity, Kidane 2026).
    ratio = BA_RATIO                 # = sqrt(2), derived (NOT fit)
    theta_L = THETA_BRANNEN          # = 2/9, derived (NOT fit)

    # Brannen formula: √mₖ = A[1 + (B/A)cos(θ + 2πk/3)]
    # where A = ½ Λ_conf = ½ M_Pl exp(−S₀) is the confinement scale.
    # k=0 → tau (heaviest), k=1 → electron (lightest), k=2 → muon.
    lepton_tree = {}
    lepton_labels = [(1, "e"), (2, "mu"), (0, "tau")]
    for k, name in lepton_labels:
        gap = 1.0 + ratio * math.cos(theta_L + 2.0 * math.pi * k / 3.0)
        lepton_tree[name] = 0.5 * M_Pl_MeV * math.exp(-S_0) * gap * gap

    # Read the pole-frame mass: tree value times QED vertex factor.
    # QED_FACTOR = 1 − α(0)/(2π) with α(0) = 1/137.036.
    #
    # α(0) is itself a reading of the same coupling at the IR pole-mass
    # scale, after the (7,26) bridge self-interference (see constants.py
    # and alpha_bridge.py for the derivation).  The bridge is a marginal
    # (h_bridge = 1) primary with unit coupling (D²_local = 1) in a
    # topological coset (c_coset = 0); the one-loop self-interference
    # reads h_bridge / (2π) = 1/(2π) exactly.  α(0) and α_em = π/512
    # are two readings of one coupling at two emergence layers, not
    # a bare value and a dressed value.
    #
    # The framework has no dimensionless tunable inputs, so the
    # vertex factor is not a tuning lever: every algebraic input
    # (B/A, θ, S₀, α_em, α(0)) is fixed by the same embedding.  Two
    # legitimate readings of the lepton spectrum coexist:
    #   * tree-level (no vertex factor): uniform +0.12% residual,
    #     reported in "The Innocent Lepton" Table 1;
    #   * pole frame (vertex factor applied): residual 0.002%,
    #     reported in "One Substrate, Three Generations" Table II.
    # Both are predictions of the same algebra at different reading
    # depths.  The same multiplier hits all three masses, so it
    # cannot move m_μ/m_e or m_τ/m_μ -- Brannen Z₃ closes the ratios
    # on its own.
    lepton_pred = {n: m * QED_FACTOR for n, m in lepton_tree.items()}

    print(f"  √mₖ = A [1 + √2 cos(2/9 + 2πk/3)]")
    print(f"  QED correction: multiply by 1 − α/(2π) = {QED_FACTOR:.9f}")
    print()
    for name in ['e', 'mu', 'tau']:
        tree = lepton_tree[name]
        pred = lepton_pred[name]
        pdg = PDG_MASSES[name]
        err = pct(pred, pdg)
        print(f"  m_{name:<3s} = {pred:12.6f} MeV   (tree: {tree:12.6f}, PDG: {pdg:12.6f},  {err:+.3f}%)")

    return {
        'lepton_tree': lepton_tree,
        'lepton_pred': lepton_pred,
        'm_e': lepton_pred['e'],
        'm_mu': lepton_pred['mu'],
        'm_tau': lepton_pred['tau'],
    }
