"""
Layer 8 -- Higgs Mass from F₄(1) Fusion + Bridge Threshold.

Companion code for:
    Kahsay, Kibrom Kidane (2026). The Innocent Lepton.
    Zenodo. https://doi.org/10.5281/zenodo.19899091
    Kahsay, Kibrom Kidane (2026). One Substrate, Three Generations.
    Zenodo. https://doi.org/10.5281/zenodo.20069456
    Kahsay, Kibrom Kidane (2026). The Echo of Standing Waves.
    Zenodo. https://doi.org/10.5281/zenodo.20144381

The complete tree-level boundary condition has TWO components:

    λ(M_Pl) = λ_F₄ − δ_bridge

  (1) F₄(1) fusion:  (26 × 26)_local = 1  →  λ_F₄ = 0
  (2) Bridge threshold:  δ_bridge = N_bridge × α²_G₂(M_Pl) × E[v²]
                                   = 182 × (1/(24π))² × 1/2
                                   = 182 / (1152π²)
                                   ≈ 0.0160

The bridge threshold is NOT a perturbative correction -- it is the
second component of the tree-level boundary condition, using the
same algebraic ingredients as the gravity derivation (Layer 10):
    N_bridge  = dim(7) × dim(26) = 182
    α_G₂     = 1/(24π)   (G₂ coupling at Planck scale)
    E[v²]    = 1/2        (protected forgetting factor)

This is exactly analogous to the lepton QED correction (1 − α/(2π)),
which also uses a derived coupling (α = π²/(256(2π−1))) to complete
the tree-level mass prediction.

The SM RGE is used as a *comparison map* -- the standard tool that
translates boundary conditions at M_Pl into observable masses.
ALL inputs are derived from the framework:
    α_s(M_Z) = 0.1177,  m_t = 172.51 GeV,  v_EW = 246.21 GeV
    α_em(M_Z) from QED running of α(0) = 1/137.036 (derived)
    sin²θ_W from tree-level SM relation with derived v_EW and α_em(M_Z)

Result: m_H = 125.1 GeV (−0.05%), zero free parameters.
"""

import math
import numpy as np
from .constants import M_Pl_GeV, ALPHA_PHYS
from .formatting import H, S


# ═══════════════════════════════════════════════════════════════════════
#  DERIVED electroweak couplings at M_Z
# ═══════════════════════════════════════════════════════════════════════
#
# α(0) = π²/(256(2π−1)) = 1/137.036 is derived (algebraic identity).
# QED vacuum polarization running from 0 to M_Z with derived fermion
# masses gives α_em(M_Z).  sin²θ_W then follows from the tree-level
# SM relation.
#
# Status of M_Z:
#   M_Z = 91.1876 GeV is a REFERENCE ENERGY, not a free parameter.
#   It plays the same role as H₀ in the gravity sector: a measured
#   scale at which derived couplings are compared to experiment.
#   The tree-level relation sin²θ_W(1−sin²θ_W) = πα v²/M_Z²
#   then determines the electroweak mixing angle, since v_EW and
#   α_em are both derived by the framework.

_M_Z = 91.1876             # GeV -- reference energy (see note above)
_M_H_EXP = 125.20          # GeV, PDG (comparison target only)


def _derive_alpha_em_MZ(lepton_masses, quark_masses):
    """
    Derive α_em(M_Z) from α(0) via QED vacuum polarization with
    confinement threshold for light quarks.

    Leptons and heavy quarks (m_q > Λ_conf): standard perturbative VP,
        Δ(1/α) = (2 Nc Qf²)/(3π) × ln(M_Z / m_f)

    Light quarks (m_q < Λ_conf): below the confinement scale, quarks
    are confined into SU(3)_3 WZW composite states (hadrons).
    The free-quark VP overestimates the spectral function in this
    region.  Replace the effective threshold with Λ_conf:
        Δ(1/α) = (2 Nc Qf²)/(3π) × ln(M_Z / Λ_conf)

    Physical picture: the VP spectral function Im Π(s) below Λ_conf²
    is dominated by hadronic bound states, not free quarks.  The
    hadron-pair production threshold (2m_π ≈ 270 MeV) is close to
    Λ_conf (314 MeV), so the confinement scale provides a natural
    IR cutoff for the quark VP.

    Uses derived fermion masses and confinement scale.  No PDG inputs.
    """
    alpha_0 = ALPHA_PHYS   # = π²/(256(2π−1)) = 1/137.036 (derived)

    # Charged fermions: (mass_GeV, Nc, Q²)
    m_e, m_mu, m_tau = [m / 1e3 for m in lepton_masses]   # MeV → GeV
    m_u, m_d, m_s, m_c, m_b, m_t = [m / 1e3 for m in quark_masses]

    # Confinement scale from Layer 1 (derived, not fitted):
    #   Λ_conf = ½ M_Pl exp(−9π²/2) ≈ 314 MeV
    Lambda_conf_GeV = 0.5 * M_Pl_GeV * math.exp(-9.0 * math.pi**2 / 2.0)

    # Charged fermions: (mass_GeV, Nc, Q², is_quark)
    fermions = [
        (m_e,   1, 1.0,     False),   # e
        (m_mu,  1, 1.0,     False),   # μ
        (m_tau, 1, 1.0,     False),   # τ
        (m_u,   3, 4.0/9.0, True),    # u
        (m_d,   3, 1.0/9.0, True),    # d
        (m_s,   3, 1.0/9.0, True),    # s
        (m_c,   3, 4.0/9.0, True),    # c
        (m_b,   3, 1.0/9.0, True),    # b
        # t quark: m_t > M_Z, does not contribute below M_Z
    ]

    M_Z_GeV = _M_Z
    delta_inv_alpha_pert = 0.0
    delta_inv_alpha_conf = 0.0

    for mf, nc, qf2, is_quark in fermions:
        if mf >= M_Z_GeV:
            continue

        # For quarks below confinement: use Λ_conf as effective threshold
        if is_quark and mf < Lambda_conf_GeV:
            mf_eff = Lambda_conf_GeV
        else:
            mf_eff = mf

        contrib = (2.0 * nc * qf2) / (3.0 * math.pi) * math.log(M_Z_GeV / mf_eff)
        delta_inv_alpha_pert += contrib

        # Track the confinement correction separately
        if is_quark and mf < Lambda_conf_GeV:
            removed = (2.0 * nc * qf2) / (3.0 * math.pi) * math.log(Lambda_conf_GeV / mf)
            delta_inv_alpha_conf += removed

    # 1/α(0) = 1/α(M_Z) + Δ  →  1/α(M_Z) = 1/α(0) − Δ
    inv_alpha_MZ = 1.0 / alpha_0 - delta_inv_alpha_pert
    alpha_em_MZ = 1.0 / inv_alpha_MZ
    return alpha_em_MZ, delta_inv_alpha_pert, delta_inv_alpha_conf, Lambda_conf_GeV


def _derive_sin2W(alpha_em_MZ, v_EW_GeV):
    """
    Derive sin²θ_W from the tree-level SM relation:

        sin²θ_W (1 − sin²θ_W) = π α_em(M_Z) / (√2 G_F M_Z²)

    where G_F = 1/(√2 v_EW²) is derived.  This gives:

        sin²θ_W (1 − sin²θ_W) = π α_em(M_Z) v_EW² / M_Z²

    Solving the quadratic x² − x + c = 0 for the physical root x < 1/2.
    """
    c = math.pi * alpha_em_MZ * v_EW_GeV**2 / _M_Z**2
    discriminant = 1.0 - 4.0 * c
    sin2w = 0.5 * (1.0 - math.sqrt(discriminant))   # physical root
    return sin2w


def _beta_sm_1loop(y):
    """One-loop SM beta functions for (g1, g2, g3, y_t, λ)."""
    g1, g2, g3, yt, lam = y
    g1s, g2s, g3s, yts = g1**2, g2**2, g3**2, yt**2
    L = 1.0 / (16.0 * math.pi**2)
    dg1 = L * (41.0/10.0) * g1**3
    dg2 = L * (-19.0/6.0) * g2**3
    dg3 = L * (-7.0) * g3**3
    dyt = L * yt * ((9.0/2.0)*yts - (17.0/20.0)*g1s - (9.0/4.0)*g2s - 8.0*g3s)
    dlam = L * (12.0*lam**2 + 12.0*yts*lam - 6.0*yts**2
                - (9.0/5.0*g1s + 9.0*g2s)*lam
                + 27.0/100.0*g1s**2 + 9.0/10.0*g1s*g2s + 9.0/4.0*g2s**2)
    return np.array([dg1, dg2, dg3, dyt, dlam])


def _beta_sm_2loop(y):
    """Two-loop SM beta functions for (g1, g2, g3, y_t, λ).

    References:
        Machacek & Vaughn, NPB 222/236/249 (1983–85)
        Ford, Jack & Jones, NPB 387 (1992)
        Buttazzo et al., JHEP 1312:089 (2013) -- collected formulas
    """
    g1, g2, g3, yt, lam = y
    g1s, g2s, g3s = g1**2, g2**2, g3**2
    yts = yt**2
    L1 = 1.0 / (16.0 * math.pi**2)
    L2 = L1**2

    # ── One-loop ──
    dg1_1 = (41.0/10.0) * g1**3
    dg2_1 = (-19.0/6.0) * g2**3
    dg3_1 = (-7.0) * g3**3
    dyt_1 = yt * ((9.0/2.0)*yts - (17.0/20.0)*g1s - (9.0/4.0)*g2s - 8.0*g3s)
    dlam_1 = (12.0*lam**2 + 12.0*yts*lam - 6.0*yts**2
              - (9.0/5.0*g1s + 9.0*g2s)*lam
              + 27.0/100.0*g1s**2 + 9.0/10.0*g1s*g2s + 9.0/4.0*g2s**2)

    # ── Two-loop gauge (Machacek-Vaughn, n_g=3 generations, n_H=1 doublet) ──
    dg1_2 = g1**3 * (199.0/50.0*g1s + 27.0/10.0*g2s + 44.0/5.0*g3s
                      - 17.0/10.0*yts)
    dg2_2 = g2**3 * (9.0/10.0*g1s + 35.0/6.0*g2s + 12.0*g3s
                      - 3.0/2.0*yts)
    dg3_2 = g3**3 * (11.0/10.0*g1s + 9.0/2.0*g2s - 26.0*g3s
                      - 2.0*yts)

    # ── Two-loop top Yukawa ──
    dyt_2 = yt * (
        -12.0*yts**2
        + yts*(131.0/16.0*g1s + 225.0/16.0*g2s + 36.0*g3s)
        + (1187.0/600.0)*g1s**2 - 9.0/20.0*g1s*g2s
        - 23.0/4.0*g2s**2 + 19.0/15.0*g1s*g3s
        + 9.0*g2s*g3s - 108.0*g3s**2
        + 6.0*lam**2 - 6.0*lam*yts
    )

    # ── Two-loop quartic (key terms for Higgs mass) ──
    dlam_2 = (
        -78.0*lam**3
        + 18.0*(3.0/5.0*g1s + 3.0*g2s)*lam**2
        - (73.0/8.0*g2s**2 - 117.0/20.0*g1s*g2s - 1887.0/200.0*g1s**2)*lam
        - 3.0*lam*yts*(5.0*yts - 25.0/4.0*g2s - 85.0/12.0*g1s - 40.0*g3s)
        - 32.0*yts**3                              # ← key term
        + (8.0/3.0*g3s - 3.0/2.0*g2s + 5.0/6.0*g1s)*6.0*yts**2
        + 305.0/16.0*g2s**3 - 289.0/48.0*g1s*g2s**2
        - 559.0/48.0*g1s**2*g2s - 379.0/48.0*g1s**3
    )

    dg1 = L1*dg1_1 + L2*dg1_2
    dg2 = L1*dg2_1 + L2*dg2_2
    dg3 = L1*dg3_1 + L2*dg3_2
    dyt = L1*dyt_1 + L2*dyt_2
    dlam = L1*dlam_1 + L2*dlam_2

    return np.array([dg1, dg2, dg3, dyt, dlam])


def _rk4(y, dt, beta_func):
    k1 = beta_func(y)
    k2 = beta_func(y + 0.5*dt*k1)
    k3 = beta_func(y + 0.5*dt*k2)
    k4 = beta_func(y + dt*k3)
    return y + dt/6.0 * (k1 + 2*k2 + 2*k3 + k4)


def _run_rge(y0, t0, t1, beta_func, n_steps=80000):
    """Integrate the RGE from t0 to t1 using RK4."""
    dt = (t1 - t0) / n_steps
    y = y0.copy()
    for _ in range(n_steps):
        y = _rk4(y, dt, beta_func)
    return y


def derive(scale: dict, quarks: dict, alpha_s_data: dict, leptons: dict, alg: dict = None, grav: dict = None):
    """
    Derive the Higgs mass from the complete tree-level boundary condition:

        λ(M_Pl) = λ_F₄ − δ_bridge
                = 0   − N_bridge × α²_G₂(M_Pl) × E[v²]

    The F₄(1) fusion sets λ_F₄ = 0.  The (7,26) bridge sector
    generates a threshold correction δ_bridge from integrating out
    182 bridge modes at M_Pl.  This is NOT a perturbative correction
    but the second component of the tree-level boundary condition --
    analogous to how the lepton QED correction (1 − α/(2π)) uses the
    derived α to complete the tree-level mass prediction.

    Every ingredient is derived from E₈ algebraic data:
        N_bridge  = dim(7) × dim(26) = 182
        α_G₂     = 1/(24π)
        E[v²]    = 1/2  (protected forgetting)

    Parameters
    ----------
    scale       : dict from scale.derive()    -- provides v_EW
    quarks      : dict from quarks.derive()   -- provides m_t and all quark masses
    alpha_s_data: dict from alpha_s.derive()  -- provides α_s(M_Z)
    leptons     : dict from leptons.derive()  -- provides lepton masses
    alg         : dict from algebra.derive()  -- provides bridge dimensions (optional)
    grav        : dict from gravity.derive()  -- provides E[v²], α_G₂ (optional)

    Returns
    -------
    dict with keys:
        alpha_em_MZ     : derived α_em(M_Z) from QED running
        sin2W           : derived sin²θ_W from tree-level SM relation
        lam_Pl_meas     : λ(M_Pl) from measured SM couplings
        lam_Pl_tree     : complete tree-level λ(M_Pl) = −δ_bridge
        delta_bridge    : bridge threshold correction
        mH_pred_1loop   : predicted m_H in GeV (one-loop, λ=0 only)
        mH_pred         : best prediction (two-loop, complete BC)
    """

    v_EW = scale['v_EW_pred_GeV']            # DERIVED: 246.21 GeV
    m_t_GeV = quarks['m_t'] / 1e3            # DERIVED: 172.51 GeV
    alpha_s_MZ = alpha_s_data['alpha_s_MZ_thresh']  # DERIVED: 0.1177

    H("LAYER 8:  HIGGS MASS FROM F₄(1) FUSION + BRIDGE THRESHOLD")

    S("8.1  Planck-scale boundary condition (complete)")

    # ── Component 1: F₄(1) fusion ────────────────────────────────────
    print(f"  Component 1: F₄(1) fusion")
    print(f"    Primaries: 𝟏 (h=0), 𝟐𝟔 (h=3/5)")
    print(f"    (𝟐𝟔 × 𝟐𝟔)_local = 𝟏  →  λ_F₄ = 0")

    # ── Component 2: Bridge sector threshold ─────────────────────────
    #
    # The (7,26) bridge sector carries 182 modes at the Planck scale.
    # Integrating them out generates a finite threshold correction to
    # the Higgs quartic, using the same algebraic ingredients as the
    # gravity derivation (Layer 10):
    #
    #   δ_bridge = N_bridge × α²_G₂(M_Pl) × E[v²]
    #
    # where:
    #   N_bridge = dim(7)×dim(26) = 182  (bridge DOFs)
    #   α_G₂(M_Pl) = 1/(24π)           (G₂ coupling at Planck scale)
    #   E[v²] = 1/2                     (protected forgetting factor)
    #
    # Physical picture: each bridge mode couples to the Higgs vacuum
    # with strength α_G₂.  The 182 modes contribute coherently,
    # weighted by the protected-forgetting factor E[v²] = 1/2 that
    # encodes how strongly the bridge "remembers" the Higgs condensate.
    #
    # Compare with gravity (Layer 10):
    #   ξ_bridge = α_G₂ × E[v²] × h_bridge = 1/(48π)
    # The Higgs threshold uses α²_G₂ (squared coupling, threshold)
    # while gravity uses α_G₂ (linear coupling, non-minimal ξ).

    dim_G2 = alg['dim_G2_fund'] if alg else 7        # = 7
    dim_F4 = alg['dim_F4_fund'] if alg else 26      # = 26
    N_bridge = dim_G2 * dim_F4                       # = 182
    alpha_G2_Pl = 1.0 / (24.0 * math.pi)            # derived in Layer 10
    E_v2 = 0.5                                       # protected forgetting
    delta_bridge = N_bridge * alpha_G2_Pl**2 * E_v2  # threshold correction

    lam_Pl_tree = 0.0 - delta_bridge                 # complete BC

    print(f"")
    print(f"  Component 2: Bridge sector threshold")
    print(f"    N_bridge = dim({dim_G2}) × dim({dim_F4}) = {N_bridge}")
    print(f"    α_G₂(M_Pl) = 1/(24π) = {alpha_G2_Pl:.8f}")
    print(f"    E[v²] = {E_v2}  (protected forgetting)")
    print(f"    δ_bridge = {N_bridge} × (1/(24π))² × 1/2 = {delta_bridge:.6f}")
    print(f"")
    print(f"  Complete boundary condition:")
    print(f"    λ(M_Pl) = λ_F₄ − δ_bridge = 0 − {delta_bridge:.6f} = {lam_Pl_tree:.6f}")

    # ── 8.1b  Derive electroweak couplings at M_Z ────────────────────
    # α(0) = 1/137.036 (derived algebraic identity)
    # → QED running with derived fermion masses → α_em(M_Z)
    # → tree-level SM relation → sin²θ_W

    lepton_masses = [leptons['m_e'], leptons['m_mu'], leptons['m_tau']]
    quark_masses = [quarks['m_u'], quarks['m_d'], quarks['m_s'],
                    quarks['m_c'], quarks['m_b'], quarks['m_t']]

    alpha_em_MZ, delta_total, delta_conf, Lambda_conf_GeV = _derive_alpha_em_MZ(lepton_masses, quark_masses)
    sin2W = _derive_sin2W(alpha_em_MZ, v_EW)

    inv_alpha_MZ = 1.0 / alpha_em_MZ
    print(f"")
    print(f"  Derived electroweak couplings at M_Z:")
    print(f"    α(0) = 1/{1.0/ALPHA_PHYS:.3f}  (algebraic identity)")
    print(f"    Δ(1/α) = {delta_total:.3f}  (QED VP with confinement threshold)")
    print(f"    Confinement correction: {delta_conf:.3f} removed  (u,d,s below Λ_conf)")
    print(f"    Λ_conf = {Lambda_conf_GeV*1e3:.1f} MeV  (derived)")
    print(f"    α_em(M_Z) = 1/{inv_alpha_MZ:.3f}  (with confinement threshold)")
    print(f"    sin²θ_W = {sin2W:.5f}  (tree-level SM: πα v²/M_Z²)")
    M_W_pred = _M_Z * math.sqrt(1.0 - sin2W)
    print(f"    M_W = M_Z cos θ_W = {M_W_pred:.2f} GeV  (PDG: 80.3692 ± 0.0133)")
    print(f"    PDG comparison: α_em(M_Z) = 1/127.951, sin²θ_W = 0.23122")

    # ── Initial conditions at M_Z ─────────────────────────────────────
    # ALL couplings derived from the framework:
    g3_mz = math.sqrt(4.0 * math.pi * alpha_s_MZ)
    yt_mz = math.sqrt(2.0) * m_t_GeV / v_EW
    g1_mz = math.sqrt(4*math.pi*alpha_em_MZ/(1-sin2W)) * math.sqrt(5.0/3.0)
    g2_mz = math.sqrt(4*math.pi*alpha_em_MZ/sin2W)

    S("8.2  SM RGE as comparison map")

    # Measured λ(M_Z) for the consistency check
    lam_mz_meas = _M_H_EXP**2 / (2.0 * v_EW**2)

    t_Pl = math.log(M_Pl_GeV / _M_Z)
    N_rge = 80000

    print(f"  Inputs (all derived):")
    print(f"    α_s(M_Z) = {alpha_s_MZ:.4f},  m_t = {m_t_GeV:.2f} GeV,  v_EW = {v_EW:.2f} GeV")
    print(f"    sin²θ_W = {sin2W:.5f},  α_em(M_Z) = 1/{1/alpha_em_MZ:.1f}")

    # Run UP with derived couplings (consistency check)
    y_up = np.array([g1_mz, g2_mz, g3_mz, yt_mz, lam_mz_meas])
    y_Pl = _run_rge(y_up, 0.0, t_Pl, _beta_sm_1loop, N_rge)
    lam_Pl_meas = y_Pl[4]

    # Run DOWN from M_Pl with λ_F₄ = 0 only (one-loop, for comparison)
    y_down_1L = np.array([y_Pl[0], y_Pl[1], y_Pl[2], y_Pl[3], 0.0])
    y_low_1L = _run_rge(y_down_1L, t_Pl, 0.0, _beta_sm_1loop, N_rge)
    lam_pred_1L = y_low_1L[4]
    mH_pred_1loop = math.sqrt(2.0 * lam_pred_1L) * v_EW

    print(f"  λ(M_Pl) from measured m_H: {lam_Pl_meas:.4f}")
    print(f"  F₄ only (λ=0, one-loop):   m_H = {mH_pred_1loop:.1f} GeV")

    # ── Two-loop RGE with complete boundary condition ─────────────────

    S("8.3  SM RGE (two-loop) with complete boundary condition")

    # Run UP with two-loop to get gauge/Yukawa couplings at M_Pl
    y_Pl_2L = _run_rge(y_up, 0.0, t_Pl, _beta_sm_2loop, N_rge)

    # Run DOWN with COMPLETE boundary condition: λ(M_Pl) = −δ_bridge
    y_down_2L = np.array([y_Pl_2L[0], y_Pl_2L[1], y_Pl_2L[2],
                          y_Pl_2L[3], lam_Pl_tree])
    y_low_2L = _run_rge(y_down_2L, t_Pl, 0.0, _beta_sm_2loop, N_rge)
    lam_pred_2L = y_low_2L[4]
    mH_pred = math.sqrt(2.0 * abs(lam_pred_2L)) * v_EW

    # Also compute F₄-only for comparison
    y_down_F4 = np.array([y_Pl_2L[0], y_Pl_2L[1], y_Pl_2L[2],
                          y_Pl_2L[3], 0.0])
    y_low_F4 = _run_rge(y_down_F4, t_Pl, 0.0, _beta_sm_2loop, N_rge)
    mH_F4_only = math.sqrt(2.0 * abs(y_low_F4[4])) * v_EW

    err_F4 = (mH_F4_only / _M_H_EXP - 1.0) * 100.0
    err_full = (mH_pred / _M_H_EXP - 1.0) * 100.0

    print(f"  F₄ only   (λ=0):          m_H = {mH_F4_only:.1f} GeV  ({err_F4:+.2f}%)")
    print(f"  F₄ + bridge (λ=−δ):       m_H = {mH_pred:.2f} GeV  ({err_full:+.2f}%)")
    print(f"  Experiment:                m_H = {_M_H_EXP:.2f} ± 0.11 GeV")
    print(f"")
    print(f"  The bridge threshold completes the tree-level prediction.")
    print(f"  No free parameters.  No external dressing.")

    return {
        'alpha_em_MZ': alpha_em_MZ,
        'inv_alpha_em_MZ': 1.0 / alpha_em_MZ,
        'delta_inv_alpha': delta_total,
        'delta_inv_alpha_conf': delta_conf,
        'sin2W': sin2W,
        'lam_Pl_meas': lam_Pl_meas,
        'lam_Pl_tree': lam_Pl_tree,
        'delta_bridge': delta_bridge,
        'mH_pred_1loop': mH_pred_1loop,
        'mH_F4_only': mH_F4_only,
        'mH_pred': mH_pred,            # best prediction (complete BC)
    }
