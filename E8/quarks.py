"""
Layer 3 -- Quark Masses from the F₄ sector.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

Six quark masses from four derivation mechanisms:
  • Light quarks: F₄ triality + conformal zero-point h(1,0) = 2/9
  • Charm:  WZW emergence  m_c = (d₁₀²d₁₁ + δ) m_μ = (217/18) m_μ
  • Top:    WZW emergence  m_t = (d₁₁⁴+d₁₀⁴+1/(2K)) m_τ = (1165/12) m_τ
  • Bottom: Q(c,b,t) = 2/3 + h₁₁/K³  from Albert algebra + WZW
  • Strange: Rivero inverse with Q(s,c,b) = 2/3 + h₁₀/K³, Albert–Dynkin bridge

Physical motivation (readable without the papers)
──────────────────────────────────────────────────

WHY quark masses are multiples of lepton masses:

  The 26 of F₄ decomposes under Spin(8) ⊂ Spin(9) ⊂ F₄ as
  26 = 8_v ⊕ 8_s ⊕ 8_c ⊕ 1 ⊕ 1.  Via the Pati–Salam embedding
  SU(4)×SU(2) ⊂ Spin(8), the three octets carry SM matter:
  8_v → charged leptons, 8_s → up quarks, 8_c → down quarks.

  The Z₃ centre of SU(3)_c forces all U(1)_em eigenvalues on the 26
  to be multiples of 1/3.  In a Z₃-equivariant standing wave, the
  allowed boundary amplitude for each species is its charge integer,
  and mass is the square of that amplitude.  Normalising by the Z₃
  quantum 1/3 gives the ratio √m_e : √m_u : √m_d = 1 : 2 : 3.

  Squaring: m_u/m_e = 4,  m_d/m_e = 9.  The conformal zero-point
  h(1,0) = 2/9 from the SU(3)₃ WZW fundamental primary adds an
  additive correction, giving m_u = (4 + 2/9)m_e = (38/9)m_e and
  m_d = (9 + 2/9)m_e = (83/9)m_e.

WHY quantum dimensions set the mass hierarchy (the WZW emergence pattern):

  Each up-type quark is shaped by its generation lepton through
  the quantum-dimension data of SU(3)₃.  The structural integers
  in each formula are products of Verlinde quantum dimensions
  d(1,0) = 2 and d(1,1) = 3:

  Gen 1:  base = d(1,0)² = 4         corr = h(1,0) = 2/9
  Gen 2:  base = d(1,0)²·d(1,1) = 12 corr = δ = 1/18   (OPE exponent)
  Gen 3:  base = d(1,1)⁴+d(1,0)⁴ = 97  corr = 1/(2K) = 1/12

  The correction ratios encode exact WZW data:
    corr₁/corr₂ = (2/9)/(1/18) = 4 = d(1,0)²
    corr₂/corr₃ = (1/18)/(1/12) = 2/3 = Q₀ (Koide parameter)

  Physical picture: the WZW fusion algebra controls how many
  independent channels a standing-wave knot can occupy at each
  generation.  The quantum dimension d(λ) counts the effective
  number of fusion channels in representation λ.  Higher generations
  access more channels (higher powers of d), producing larger masses.

WHY the Koide parameter Q governs heavier quarks:

  The Koide ratio Q = (m₁+m₂+m₃)/(√m₁+√m₂+√m₃)² = 1/3 + |B/A|²/6.
  The octonionic CG calculation (octonions.py) gives |B/A|² = 2 →
  Q₀ = 2/3 for all Z₃ Brannen triplets.  The sub-leading WZW
  correction is Q = 2/3 + h(rep)/K³, where the representation λ is
  uniquely selected by each triplet's closure:
    • (e,μ,τ):  identity (0,0), h=0  → Q = 2/3 exactly
    • (c,b,t):  adjoint (1,1), h=1/2 → Q = 289/432
    • (s,c,b):  fundamental (1,0), h=2/9 → Q = 649/972

WHY the bridge factor √(32/27) appears in the strange mass:

  The Albert algebra J₃(𝕆) has a natural trace norm whose spread
  parameter is δ_J = √(3/8) = √(d(1,1)/d(1,0)³).  The Dynkin Z₂
  fixed point (from the triality section) is δ_{Z₂} = 2/3.  The
  ratio bridge² = δ_{Z₂}²/δ_J² = (2/3)²/(3/8) = Q₀²·d(1,0)³/d(1,1)
  = 32/27 connects the Albert algebra geometry to the WZW data.
  The strange quark enters through the antisymmetric ε-channel of
  SU(3), which flips the sign of √m_s in the Koide sum (the Z₂
  outer automorphism of SU(3) sends ε_{ijk} → −ε_{ijk}).
"""

import math
from .constants import PDG_MASSES
from .formatting import H, S, pct


def derive(alg: dict, scale: dict, lep: dict):
    """
    Derive all six quark masses.

    Parameters
    ----------
    alg   : dict from algebra.derive()
    scale : dict from scale.derive()
    lep   : dict from leptons.derive()

    Returns
    -------
    dict with keys:
        m_u, m_d, m_c, m_t, m_b, m_s           : masses in MeV
        m_t_tree                                : Yukawa cross-check v_EW/√2 in MeV
        Q_val                                   : corrected Q-value
        h_fund, delta_OPE                       : conformal data
        bridge                                  : Albert–Dynkin bridge factor
    """

    m_e  = lep['m_e']
    m_mu = lep['m_mu']
    m_tau = lep['m_tau']
    v_EW_MeV = scale['v_EW_pred_MeV']

    H("LAYER 3:  QUARK MASSES -- F₄ SECTOR")

    # ──────────────────────────────────────────────────────────────────
    #  No-tuning note for the entire layer
    # ──────────────────────────────────────────────────────────────────
    # The six quark masses below look like polynomials in small integers
    # over the SU(3)₃ WZW data:
    #
    #     m_u   = (4   + 2/9 ) m_e            light triality + h(1,0)
    #     m_d   = (9   + 2/9 ) m_e            light triality + h(1,0)
    #     m_c   = (12  + 1/18) m_μ            d₁₀²·d₁₁ + δ_OPE
    #     m_t   = (97  + 1/12) m_τ            d₁₁⁴+d₁₀⁴ + 1/(2K)
    #     Q(c,b,t) = 2/3 + h₁₁/K³ = 289/432   adjoint-channel Koide
    #     Q(s,c,b) = 2/3 + h₁₀/K³ = 649/972   fundamental-channel Koide
    #
    # Each of the integers above is a Verlinde quantum dimension or
    # altitude of the unique SU(3)₃ modular tensor category that the
    # conformal embedding E₈(1) ⊃ G₂(1) × F₄(1) at level 1 selects:
    #
    #     d(1,0) = 2          (quantum dim of fundamental)
    #     d(1,1) = 3          (quantum dim of adjoint)
    #     K      = 6          (altitude k + h∨)
    #     h(1,0) = 2/9        (conformal weight of fundamental)
    #     h(1,1) = 1/2        (conformal weight of adjoint)
    #
    # None of {2, 3, 6, 2/9, 1/2} is adjustable — they are operator
    # identities of the WZW model.  And the ladder of corrections
    # satisfies its own self-consistency:
    #
    #     corr₁/corr₂ = (2/9)/(1/18)  = 4   = d(1,0)²
    #     corr₂/corr₃ = (1/18)/(1/12) = 2/3 = Q₀   (the Koide value)
    #
    # so generations 1→2→3 are not three independent ansätze but one
    # ladder, and the ratio of consecutive corrections reproduces the
    # Koide parameter Q₀ that already appears in the lepton sector.
    # There is no free coefficient anywhere in this layer.

    # ══════════════════════════════════════════════════════════════════
    #  3.1  Light quarks: triality + conformal zero-point
    # ══════════════════════════════════════════════════════════════════

    S("3.1  Light quarks: (d² + h(1,0)) × m_e")

    # Conformal weight of the SU(3)₃ WZW fundamental primary (1,0):
    #   h(1,0) = C₂(fund, SU(3)) / (k + h∨)
    #          = (4/3) / (3 + 3) = 2/9
    # where k=3 (WZW level) and h∨=3 (dual Coxeter number of SU(3)).
    # This is the conformal zero-point energy that shifts each mass.
    h_fund = 2.0 / 9.0          # h(1,0) = (4/3)/6 = 2/9

    # OPE sub-leading exponent: the difference between the adjoint
    # conformal weight and twice the fundamental weight.
    delta_OPE = 1.0 / 18.0      # δ = h(1,1) − 2h(1,0) = 1/2 − 4/9

    # The integer coefficients 4 and 9 are NOT fitted.  They come from the
    # Z_3 centre of SU(3)_c, which forces U(1)_em eigenvalues on the 26 of F_4
    # to be integer multiples of 1/3: charges {0, 1/3, 2/3, 1} for {nu, d, u, e}.
    # The Dynkin Z_2 swap (8_v <-> 8_c, i.e., e_R <-> d_R) cross-assigns
    # the charges as sqrt(m) labels: sqrt(m_e):sqrt(m_u):sqrt(m_d) = 1:2:3.
    # Squaring gives the tree ratios m_u/m_e = 4, m_d/m_e = 9.  The h_fund
    # shift adds the SU(3)_3 fundamental conformal weight (2/9).
    # Ref: One Substrate Three Generations, Sec. 5.1 (light-quark triality).
    m_u = (4.0 + h_fund) * m_e     # 4 = (charge int u)^2 = 2^2; total (38/9) m_e
    m_d = (9.0 + h_fund) * m_e     # 9 = (charge int d_R after Z_2 swap)^2 = 3^2; total (83/9) m_e

    print(f"  m_q = (d(λ)² + h(1,0)) × m_e")
    print(f"  h(1,0) = 2/9 = {h_fund:.6f}")
    print(f"  m_u = (4 + 2/9) m_e = (38/9) m_e = {m_u:.4f} MeV   (PDG: {PDG_MASSES['u']:.2f},  {pct(m_u, PDG_MASSES['u']):+.1f}%)")
    print(f"  m_d = (9 + 2/9) m_e = (83/9) m_e = {m_d:.4f} MeV   (PDG: {PDG_MASSES['d']:.2f},  {pct(m_d, PDG_MASSES['d']):+.1f}%)")
    print(f"  m_u/m_d = 38/83 = {38/83:.6f}   (FLAG: 0.459 ± 0.027)")

    # ══════════════════════════════════════════════════════════════════
    #  3.2  Charm mass: WZW emergence
    # ══════════════════════════════════════════════════════════════════

    S("3.2  Charm: (d₁₀²·d₁₁ + δ) × m_μ = (217/18) m_μ")

    # WZW emergence formula for charm (see docstring for full motivation):
    #   base = d(1,0)²·d(1,1) = 4×3 = 12
    #     This product of Verlinde quantum dimensions counts the
    #     independent fusion channels available at generation 2.
    #   correction = δ = h(1,1) − 2h(1,0) = 1/18 (OPE sub-leading exponent)
    #   Each up-type quark mass = (base + correction) × generation lepton
    m_c = (12.0 + delta_OPE) * m_mu           # (217/18) m_μ

    # Cross-check via the Dynkin Z₂ generation ladder:
    #   The Z₂ ⊂ S₃ maps the lepton generation ratio τ/μ to the quark
    #   generation ratio, carried across by the (7,26) bridge of E₈.
    tau_over_mu = m_tau / m_mu
    m_c_ladder = 9.0 * m_e * tau_over_mu**2   # Dynkin cross-check

    print(f"  m_c = (12 + 1/18) m_μ = (217/18) m_μ = {m_c:.2f} MeV   (PDG: {PDG_MASSES['c']:.0f},  {pct(m_c, PDG_MASSES['c']):+.2f}%)")
    print(f"  Correction ratio: corr₁/corr₂ = (2/9)/(1/18) = {(2.0/9)/(1.0/18):.0f} = d₁₀²  ✓")
    print(f"  Dynkin cross-check: m_c(ladder) = 9 m_e × (τ/μ)² = {m_c_ladder:.1f} MeV ({pct(m_c_ladder, PDG_MASSES['c']):+.1f}%)")

    # ══════════════════════════════════════════════════════════════════
    #  3.3  Top mass: WZW emergence (1165/12) m_τ
    # ══════════════════════════════════════════════════════════════════

    S("3.3  Top: (d₁₁⁴ + d₁₀⁴ + 1/(2K)) m_τ = (1165/12) m_τ")

    # Third-generation WZW emergence formula:
    #
    #   m_t = (d(1,1)⁴ + d(1,0)⁴ + 1/(2K)) × m_τ
    #       = (81 + 16 + 1/12) × m_τ  = (1165/12) m_τ
    #
    # Physical picture: the third generation accesses the FULL
    # fusion algebra, both representations at their maximal power.
    # d(1,1)⁴ = 81 counts the adjoint fusion channels (3⁴),
    # d(1,0)⁴ = 16 counts the fundamental channels (2⁴).
    # The correction 1/(2K) = 1/12 is the altitude suppression.
    #
    # The three-generation cascade (see docstring for full derivation):
    #   Gen 1: m_u = (d₁₀²  + h₁₀)  m_e   = (38/9) m_e     corr₁ = 2/9
    #   Gen 2: m_c = (d₁₀²d₁₁ + δ)  m_μ   = (217/18) m_μ   corr₂ = 1/18
    #   Gen 3: m_t = (d₁₁⁴+d₁₀⁴ + 1/(2K)) m_τ = (1165/12) m_τ  corr₃ = 1/12
    #
    # The correction ratios are exact WZW data:
    #   corr₁/corr₂ = (2/9)/(1/18) = 4 = d(1,0)²
    #   corr₂/corr₃ = (1/18)/(1/12) = 2/3 = Q₀ (the Koide parameter)
    d_10 = 2    # quantum dimension of (1,0)
    d_11 = 3    # quantum dimension of (1,1)
    K = 6       # altitude
    corr_3 = 1.0 / (2.0 * K)                              # = 1/12
    top_base = float(d_11**4 + d_10**4)                    # = 97
    m_t = (top_base + corr_3) * m_tau                      # (1165/12) m_τ

    print(f"  m_t = (d₁₁⁴ + d₁₀⁴ + 1/(2K)) m_τ")
    print(f"      = ({d_11}⁴ + {d_10}⁴ + 1/12) × m_τ")
    print(f"      = (1165/12) m_τ")
    print(f"      = {m_t/1e3:.2f} GeV   (PDG: {PDG_MASSES['t']/1e3:.2f},  {pct(m_t, PDG_MASSES['t']):+.3f}%)")
    print(f"  corr₂/corr₃ = (1/18)/(1/12) = 2/3 = Q₀  ✓")

    # Yukawa cross-check: y_t ≈ 1 emerges as a consequence
    m_t_tree = v_EW_MeV / math.sqrt(2.0)
    y_t = m_t / m_t_tree
    print(f"  Yukawa cross-check: v_EW/√2 = {m_t_tree/1e3:.2f} GeV → y_t = {y_t:.6f} (≈ 1 − 1/128)")

    # ══════════════════════════════════════════════════════════════════
    #  3.4  Bottom mass: Q(c,b,t) = 2/3 + h₁₁/K³
    # ══════════════════════════════════════════════════════════════════

    S("3.4  Bottom: Q(c,b,t) = 2/3 + h₁₁/K³")

    # Koide parameter for the (c,b,t) triplet:
    #
    # The triplet (c,b,t) spans the diagonal Z₃ of the Albert algebra
    # J₃(𝕆): each member sits in a different triality sector of
    # F₄ → SU(3)_A × SU(3)_B.  The amplitude ratio B/A = √2 (from
    # octonionic CG, see octonions.py) gives the universal Q₀ = 2/3.
    #
    # The sub-leading correction is the adjoint conformal weight h(1,1)
    # divided by K³ = 216 (the altitude cubed).  This is the unique
    # WZW correction for the adjoint channel through which the (c,b,t)
    # triplet closes (see the unified formula Q = 2/3 + h(rep)/K³).
    h_11 = 0.5                              # conformal weight of (1,1)
    K = 6                                   # altitude = k + h∨ = 3 + 3
    Q_tree = 2.0 / 3.0
    Q_corr = h_11 / K**3                   # h₁₁/K³ = (1/2)/216 = 1/432
    Q_val  = Q_tree + Q_corr

    sc = math.sqrt(m_c)
    st = math.sqrt(m_t)
    a_coeff = Q_val - 1.0
    b_coeff = 2.0 * Q_val * (sc + st)
    c_coeff = (Q_val - 1.0) * (m_c + m_t) + 2.0 * Q_val * sc * st
    disc = b_coeff**2 - 4.0 * a_coeff * c_coeff
    x1 = (-b_coeff + math.sqrt(disc)) / (2.0 * a_coeff)
    x2 = (-b_coeff - math.sqrt(disc)) / (2.0 * a_coeff)
    m_b = x1**2 if m_c < x1**2 < m_t else x2**2

    Q_cbt = (m_c + m_b + m_t) / (math.sqrt(m_c) + math.sqrt(m_b) + math.sqrt(m_t))**2

    print(f"  Q = 2/3 + h₁₁/K³ = 2/3 + (1/2)/216 = 289/432 = {Q_val:.10f}")
    print(f"  m_b = {m_b:.1f} MeV   (PDG: {PDG_MASSES['b']:.0f},  {pct(m_b, PDG_MASSES['b']):+.1f}%)")
    print(f"  Verify: Q = {Q_cbt:.10f}  (target: {Q_val:.10f})  ✓")

    # ══════════════════════════════════════════════════════════════════
    #  3.5  Strange mass: Rivero inverse + Albert–Dynkin bridge
    # ══════════════════════════════════════════════════════════════════

    S("3.5  Strange: Rivero inverse with Q = 2/3 + h₁₀/K³ + √(32/27) bridge")

    # Unified Koide formula: Q = 2/3 + h(rep)/K³
    #   (c,b,t) uses h(1,1) = 1/2  (adjoint)      → Q = 289/432
    #   (s,c,b) uses h(1,0) = 2/9  (fundamental)   → Q = 649/972
    #
    # The (s,c,b) triplet closes through the fundamental representation
    # because the strange quark enters through the antisymmetric
    # ε-channel of SU(3).  The ε_{ijk} tensor is the unique SU(3)
    # invariant in 3⊗3→3̄, and under the Z₂ outer automorphism
    # (complex conjugation: 3↔3̄) it flips sign.  This forces
    # √m_s → −√m_s in the Koide sum, the "Rivero inverse" sign.
    Q_s = 2.0 / 3.0 + h_fund / K**3    # 2/3 + (2/9)/216 = 649/972
    sc_r = math.sqrt(m_c)
    sb_r = math.sqrt(m_b)
    a_r = 1.0 - Q_s
    b_r = 2.0 * Q_s * (sc_r + sb_r)
    c_r = m_c + m_b - Q_s * (sc_r + sb_r)**2
    disc_r = b_r**2 - 4.0 * a_r * c_r
    x1_r = (-b_r + math.sqrt(disc_r)) / (2.0 * a_r)
    x2_r = (-b_r - math.sqrt(disc_r)) / (2.0 * a_r)
    x_phys = x1_r if x1_r > 0 else x2_r
    m_s_bare = x_phys**2

    # Albert–Dynkin bridge factor (see docstring for derivation):
    #
    # The Albert algebra J₃(𝕆) has trace norm spread δ_J = √(3/8),
    # which equals √(d(1,1)/d(1,0)³) in WZW language.  The Dynkin Z₂
    # fixed point from the triality section is δ_{Z₂} = 2/3.
    # The bridge connects these two geometries:
    #   bridge² = δ_{Z₂}²/δ_J² = (2/3)²/(3/8) = (4/9)×(8/3) = 32/27
    # Equivalently: Q₀² × d(1,0)³/d(1,1) = (2/3)² × 8/3 = 32/27.
    d_10 = 2    # quantum dimension of (1,0)
    d_11 = 3    # quantum dimension of (1,1)
    bridge_sq = (2.0/3.0)**2 * d_10**3 / d_11     # = 32/27
    bridge = math.sqrt(bridge_sq)                   # √(32/27) ≈ 1.0887
    assert abs(bridge_sq - 32.0/27.0) < 1e-12
    m_s = m_s_bare * bridge

    print(f"  Q_s = 2/3 + h₁₀/K³ = 2/3 + (2/9)/216 = 649/972 = {Q_s:.10f}")
    print(f"  Unified: Q = 2/3 + h(rep)/K³  [h₁₁ for (c,b,t), h₁₀ for (s,c,b)]")
    print(f"  m_s (bare Rivero) = {m_s_bare:.1f} MeV")
    print(f"  bridge² = Q₀² × d₁₀³/d₁₁ = (4/9)(8/3) = 32/27")
    print(f"  bridge = √(32/27) = {bridge:.6f}")
    print(f"  m_s = {m_s:.1f} MeV   (PDG: {PDG_MASSES['s']:.1f},  {pct(m_s, PDG_MASSES['s']):+.1f}%)")

    # ── 3.6  Nuclear stability ────────────────────────────────────────

    S("3.6  Nuclear stability check")

    delta_ud = m_d - m_u
    print(f"  m_d − m_u = {delta_ud:.3f} MeV")
    print(f"  Hogan window: [0.585, 3.4] MeV")
    print(f"  Status: {'INSIDE ✓' if 0.585 <= delta_ud <= 3.4 else 'OUTSIDE ✗'}")

    return {
        'm_u': m_u, 'm_d': m_d, 'm_c': m_c,
        'm_t': m_t, 'm_b': m_b, 'm_s': m_s,
        'm_t_tree': m_t_tree,
        'Q_val': Q_val, 'Q_cbt': Q_cbt, 'Q_s': Q_s,
        'h_fund': h_fund,
        'delta_OPE': delta_OPE,
        'bridge': bridge,
    }
