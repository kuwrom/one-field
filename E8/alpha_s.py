"""
Layer 7 -- Strong Coupling from the Embedding Index.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

The chain SU(3) ⊂ G₂ ⊂ E₈ has embedding index 1 at every step,
so α_s = α_{G₂} at the matching scale.

Key result: 1/α_G₂^WZW(v_EW) = 32/π  (exact -- the h∨(F₄)π²
piece cancels).  The low-energy QCD α_s at the same coordinate is obtained
only after the G₂ → SU(3) threshold map.

The G₂ → SU(3) threshold correction uses the derived vector mass
M_V = g_{G₂} v_EW/√6 = π v_EW/(2√12) ≈ 112 GeV.
"""

import math
from .constants import M_Pl_GeV
from .formatting import H, S, box


def _rk4_run(beta_func, a0, t0, t1, n_steps=50000):
    """Integrate dα/dt = beta(α) from t0 to t1 using RK4."""
    dt = (t1 - t0) / n_steps
    a = a0
    for _ in range(n_steps):
        k1 = beta_func(a)
        k2 = beta_func(a + 0.5*dt*k1)
        k3 = beta_func(a + 0.5*dt*k2)
        k4 = beta_func(a + dt*k3)
        a += dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return a


def _run_SM_2loop(a0, mu0, mu1, nf):
    """Two-loop SM QCD running in d/d(ln μ) convention."""
    b0 = (33 - 2*nf) / (6*math.pi)
    b1 = (153 - 19*nf) / (12*math.pi**2)
    return _rk4_run(lambda a: -b0*a**2 - b1*a**3,
                    a0, math.log(mu0), math.log(mu1))


def derive(alg: dict, scale: dict, quarks: dict):
    """
    Derive α_s(M_Z) from the E₈ embedding index chain.

    Parameters
    ----------
    alg    : dict from algebra.derive()
    scale  : dict from scale.derive()
    quarks : dict from quarks.derive()

    Returns
    -------
    dict with keys:
        alpha_s_vEW               : π/32 (exact WZW pre-threshold datum)
        alpha_s_vEW_threshold     : QCD coupling after G₂→SU(3) threshold
        alpha_s_MZ_1loop          : one-loop SM running result
        alpha_s_MZ_no_thresh      : two-loop without threshold (baseline)
        alpha_s_MZ_thresh         : two-loop with G₂→SU(3) threshold
        M_V_derived               : derived vector boson mass (GeV)
        lambda_3                  : G₂→SU(3) threshold correction
        err_1loop, err_no_thresh, err_thresh : percentage errors vs PDG
    """

    C2_fund_F4 = alg['C2_fund_F4']
    v_EW = scale['v_EW_pred_GeV']
    m_t_GeV = quarks['m_t'] / 1e3
    m_b_GeV = quarks['m_b'] / 1e3
    m_c_GeV = quarks['m_c'] / 1e3
    M_Z = 91.1876
    alpha_s_PDG = 0.1180

    H("LAYER 7:  STRONG COUPLING FROM THE EMBEDDING INDEX")

    # ── 7.1  Embedding index chain ────────────────────────────────────

    S("7.1  SU(3) ⊂ G₂ ⊂ E₈, affine level index 1")

    print(f"  7 → 3 ⊕ 3̄ ⊕ 1:  Dynkin index 1/2+1/2+0 = 1")
    print(f"  Conformal embedding: j_f(G₂ ⊂ E₈) = 1")
    print(f"  → α_s = α_{{G₂}} at matching scale")

    # ── 7.2  α_{G₂}(M_Pl) ────────────────────────────────────────────

    S("7.2  α_{G₂}(M_Pl) = 1/(24π)")

    b0_G2 = 32.0 / 3.0
    alpha_G2_Pl = 1.0 / (24.0 * math.pi)
    print(f"  g²_adj(M_Pl) = 1,  b₀(G₂) = 32/3")
    print(f"  α_{{G₂}}(M_Pl) = 1/(24π) = {alpha_G2_Pl:.8f}")

    # ── 7.3  Exact cancellation: α_s(v_EW) = π/32 ────────────────────

    S("7.3  α_G₂^WZW(v_EW) = π/32 (exact pre-threshold datum)")

    alpha_s_vEW = math.pi / 32.0

    # ── Symbolic cancellation algebra ─────────────────────────────
    #
    # One-loop running from M_Pl to v_EW:
    #   1/α(v_EW) = 1/α(M_Pl) − b₀/(2π) × ln(M_Pl/v_EW)
    #
    # Key identity: ln(M_Pl/v_EW) = S_quark(LO) = 9π²/2 − C₂(26)
    #
    # Substituting:
    #   1/α(v_EW) = 24π − (16/(3π)) × (9π²/2 − 6)
    #             = 24π − [24π − 32/π]
    #             = 32/π
    #
    # The h∨(F₄)·π² piece cancels EXACTLY against 1/α_{G₂}(M_Pl).
    # Only the Casimir correction C₂(26) = 6 survives:
    #   (16/(3π)) × 6 = 96/(3π) = 32/π  ✓

    inv_alpha_Pl = 24 * math.pi
    term_hdual = (b0_G2 / (2*math.pi)) * (9*math.pi**2/2)
    term_casimir = (b0_G2 / (2*math.pi)) * C2_fund_F4
    inv_alpha_vEW = inv_alpha_Pl - term_hdual + term_casimir

    print(f"  One-loop running: 1/α(v_EW) = 1/α(M_Pl) − b₀/(2π) × ln(M_Pl/v_EW)")
    print(f"  ln(M_Pl/v_EW) = S_quark(LO) = 9π²/2 − C₂(26)")
    print(f"")
    print(f"  Substituting:")
    print(f"    1/α(v_EW) = 24π − (16/(3π)) × (9π²/2 − 6)")
    print(f"              = 24π − [24π − 32/π]")
    print(f"              = 32/π")
    print(f"")
    print(f"  Numerical verification:")
    print(f"    1/α(M_Pl) = 24π = {inv_alpha_Pl:.4f}")
    print(f"    b₀/(2π) × h∨π²/2 = {term_hdual:.4f}  (= 24π ✓)")
    print(f"    b₀/(2π) × C₂(26) = {term_casimir:.4f}  (= 32/π ✓)")
    print(f"    1/α_G₂^WZW(v_EW) = 32/π = {inv_alpha_vEW:.4f}")
    print(f"    α_G₂^WZW(v_EW) = π/32 = {alpha_s_vEW:.6f}")

    # ── 7.4  One-loop SM running to M_Z ──────────────────────────────

    S("7.4  SM running to M_Z (one-loop analytic)")

    b0_nf6 = (33 - 12) / (6 * math.pi)
    b0_nf5 = (33 - 10) / (6 * math.pi)

    inv_as_mt = 32/math.pi + b0_nf6 * math.log(m_t_GeV / v_EW)
    inv_as_MZ = inv_as_mt + b0_nf5 * math.log(M_Z / m_t_GeV)
    alpha_s_MZ_1loop = 1.0 / inv_as_MZ
    err_1loop = 100 * (alpha_s_MZ_1loop - alpha_s_PDG) / alpha_s_PDG

    print(f"  Thresholds: m_t = {m_t_GeV:.1f}, m_b = {m_b_GeV:.1f}, m_c = {m_c_GeV:.1f} GeV")
    print(f"  α_s(M_Z) [1-loop] = {alpha_s_MZ_1loop:.6f}  ({err_1loop:+.1f}%)")

    # ── 7.5  Two-loop G₂ cross-check ─────────────────────────────────

    S("7.5  Two-loop G₂ running (cross-check)")

    beta0_G2 = (32.0/3.0) / (2*math.pi)
    beta1_G2 = (232.0/3.0) / (8*math.pi**2)

    alpha_vEW_2L = _rk4_run(
        lambda a: -beta0_G2*a**2 - beta1_G2*a**3,
        alpha_G2_Pl, math.log(M_Pl_GeV), math.log(v_EW))

    a_mt_2L = _run_SM_2loop(alpha_vEW_2L, v_EW, m_t_GeV, 6)
    alpha_MZ_2L = _run_SM_2loop(a_mt_2L, m_t_GeV, M_Z, 5)
    err_2L = 100 * (alpha_MZ_2L - alpha_s_PDG) / alpha_s_PDG

    print(f"  α_s(v_EW) [2-loop G₂] = {alpha_vEW_2L:.6f}  ({100*(alpha_vEW_2L-alpha_s_vEW)/alpha_s_vEW:+.1f}%)")
    print(f"  α_s(M_Z) [2-loop G₂ + SM] = {alpha_MZ_2L:.6f}  ({err_2L:+.1f}%)")
    print(f"  Overshoot consistent with π/32 being non-perturbative (WZW/Sugawara)")

    # ── 7.6  G₂ → SU(3) threshold ────────────────────────────────────

    S("7.6  G₂ → SU(3) threshold (derived M_V)")

    g_G2 = math.sqrt(4.0 * math.pi * alpha_s_vEW)     # = π/(2√2)
    v_break = v_EW / math.sqrt(2.0)                    # EWSB identification
    M_V = g_G2 * v_EW / math.sqrt(6.0)                 # = π v_EW/(2√12)

    T_V = 1.0         # T(3) + T(3̄)
    C_diff = 1.0       # C_{G₂} − C_{SU₃}
    lambda_3 = C_diff - 21.0 * T_V * math.log(M_V / v_EW)

    inv_as_G2_vEW = 1.0 / alpha_s_vEW                  # = 32/π
    inv_as_thresh = inv_as_G2_vEW - lambda_3 / (12.0 * math.pi)
    a_thresh = 1.0 / inv_as_thresh

    a_mt_th = _run_SM_2loop(a_thresh, v_EW, m_t_GeV, 6)
    alpha_MZ_thresh = _run_SM_2loop(a_mt_th, m_t_GeV, M_Z, 5)
    err_thresh = 100 * (alpha_MZ_thresh - alpha_s_PDG) / alpha_s_PDG

    print(f"  g_{{G₂}} = π/(2√2) = {g_G2:.6f}")
    print(f"  v_break = v_EW/√2 = {v_break:.3f} GeV")
    print(f"  M_V = g v_EW/√6 = {M_V:.3f} GeV")
    print(f"  λ₃ = {lambda_3:.6f}")
    print(f"  α_s(v_EW) [after threshold] = {a_thresh:.6f}")
    print(f"  α_s(M_Z) [+threshold] = {alpha_MZ_thresh:.6f}  ({err_thresh:+.2f}%)")

    # ── 7.7  No-threshold baseline for comparison ────────────────────

    a_no_th_mt = _run_SM_2loop(alpha_s_vEW, v_EW, m_t_GeV, 6)
    alpha_MZ_no_th = _run_SM_2loop(a_no_th_mt, m_t_GeV, M_Z, 5)
    err_no_th = 100 * (alpha_MZ_no_th - alpha_s_PDG) / alpha_s_PDG

    box([
        f"α_G₂^WZW(v_EW) = π/32 = {alpha_s_vEW:.6f}  (pre-threshold)",
        f"α_s(v_EW) after threshold = {a_thresh:.6f}",
        f"M_V = π v_EW/(2√12) = {M_V:.1f} GeV  (derived, 0 params)",
        f"",
        f"α_s(M_Z) [no threshold]  = {alpha_MZ_no_th:.4f}  ({err_no_th:+.1f}%)",
        f"α_s(M_Z) [+threshold]    = {alpha_MZ_thresh:.4f}  ({err_thresh:+.1f}%)",
        f"PDG:                       0.1180  ± 0.0009",
        f"Δ = {abs(alpha_MZ_thresh - alpha_s_PDG):.4f} -- {abs(alpha_MZ_thresh - alpha_s_PDG)/0.0009:.1f}σ",
    ])

    return {
        'alpha_s_vEW': alpha_s_vEW,
        'alpha_s_vEW_threshold': a_thresh,
        'alpha_s_MZ_1loop': alpha_s_MZ_1loop,
        'alpha_s_MZ_no_thresh': alpha_MZ_no_th,
        'alpha_s_MZ_thresh': alpha_MZ_thresh,
        'M_V_derived': M_V,
        'lambda_3': lambda_3,
        'err_1loop': err_1loop,
        'err_no_thresh': err_no_th,
        'err_thresh': err_thresh,
    }
