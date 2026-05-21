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

    # ══════════════════════════════════════════════════════════════════
    #  3.1  Light quarks: triality + conformal zero-point
    # ══════════════════════════════════════════════════════════════════

    S("3.1  Light quarks: (d² + h(1,0)) × m_e")

    h_fund = 2.0 / 9.0          # h(1,0) = C₂(fund) / (2K) = (8/3)/12
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

    m_c = (12.0 + delta_OPE) * m_mu           # (217/18) m_μ
    tau_over_mu = m_tau / m_mu
    m_c_ladder = 9.0 * m_e * tau_over_mu**2   # Dynkin cross-check

    print(f"  m_c = (12 + 1/18) m_μ = (217/18) m_μ = {m_c:.2f} MeV   (PDG: {PDG_MASSES['c']:.0f},  {pct(m_c, PDG_MASSES['c']):+.2f}%)")
    print(f"  Correction ratio: corr₁/corr₂ = (2/9)/(1/18) = {(2.0/9)/(1.0/18):.0f} = d₁₀²  ✓")
    print(f"  Dynkin cross-check: m_c(ladder) = 9 m_e × (τ/μ)² = {m_c_ladder:.1f} MeV ({pct(m_c_ladder, PDG_MASSES['c']):+.1f}%)")

    # ══════════════════════════════════════════════════════════════════
    #  3.3  Top mass: WZW emergence (1165/12) m_τ
    # ══════════════════════════════════════════════════════════════════

    S("3.3  Top: (d₁₁⁴ + d₁₀⁴ + 1/(2K)) m_τ = (1165/12) m_τ")

    # The third-generation emergence formula:
    #   m_t = (d₁₁⁴ + d₁₀⁴ + corr₃) × m_τ
    #       = (81 + 16 + 1/12) × m_τ
    #       = (1165/12) m_τ
    #
    # This completes the generation cascade:
    #   Gen 1: m_u = (d₁₀²  + h₁₀)  m_e   = (38/9) m_e     corr₁ = 2/9
    #   Gen 2: m_c = (d₁₀²d₁₁ + δ)  m_μ   = (217/18) m_μ   corr₂ = 1/18
    #   Gen 3: m_t = (d₁₁⁴+d₁₀⁴ + 1/(2K)) m_τ = (1165/12) m_τ  corr₃ = 1/12
    #
    # Correction ratios: corr₁/corr₂ = d₁₀² = 4,  corr₂/corr₃ = Q₀ = 2/3
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

    # The Koide correction uses h₁₁ (the adjoint conformal weight 1/2)
    # divided by the altitude cubed K³ = 216.  Equivalently, this is
    # the third-generation WZW correction corr₃ = 1/(2K) = 1/12
    # divided by K²:  Q = 2/3 + corr₃/K² = 2/3 + 1/432 = 289/432.
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

    # Unified formula: Q = 2/3 + h(rep)/K³
    #   (c,b,t) uses h₁₁ (adjoint weight)  → Q = 289/432
    #   (s,c,b) uses h₁₀ (fundamental weight) → Q = 649/972
    # The fundamental weight governs the ε-channel (sign-flip) sector.
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

    # Bridge factor from WZW data:
    #   bridge² = Q₀² × d₁₀³/d₁₁ = (2/3)² × 8/3 = 32/27
    # Equivalently, the Albert algebra spread 3/8 = d₁₁/d₁₀³,
    # so the Dynkin–Albert ratio δ_Z₂²/δ_J² = Q₀²/(d₁₁/d₁₀³) = Q₀² d₁₀³/d₁₁.
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
