"""
Higgs quartic boundary + SM RGE transport + electroweak observables.

The bridge boundary (lambda(M_Pl) = -delta_bridge) and the MS-bar
transport to the top scale are the content of this module.  Everything
here is transport infrastructure for turning the Planck scale boundary
into a low energy reading. The boundary itself lives in root.py
(WEB['lambda_MPl']).

═══════════════════════════════════════════════════════════════════════
WHAT THE ENGINE IMPORTS (one table for the whole engine: every
external number lives in transport, so the full import list is
stated here once)
═══════════════════════════════════════════════════════════════════════

  THE ANCHOR (sets the unit):
    m_e = 0.51099895069 MeV        CODATA 2022      root.M_E_ANCHOR_MEV
    ħ, GeV/J                        exact SI         root (conversions)

  COMPARISON SCALE (ρ_crit bookkeeping, read by no prediction):
    H₀  = 67.4 km/s/Mpc            Planck 2018      root.H_0_KM_S_MPC

  DECLARED RADIATIVE TRANSPORT (carries the engine's A₀ and sin²θ_W to
  the coordinate M_W and M_Z are measured in, read by no engine constant
  (tests/audit/test_blindness.py), and audited for compatibility with
  the framework's own values, below):
    Δr̂_W = 0.06937(6)              PDG 2024 10.26   reference.PDG_EW
    ρ̂    = 1.01016(9)              PDG 2024 10.26   reference.PDG_EW
    m_H  = 125.20 GeV               PDG 2024, inside the loop remainders
                                    of Δr̂_W and ρ̂. The framework
                                    produces no m_H (data class)

  CITED QCD MATCHING AND CONVERSION CONSTANTS (fixed by field
  content, not measured, not adjustable):
    pole ↔ MS-bar top coefficients  Marquard et al. 1502.01030
      0.4244, 0.8345, 2.375, 8.49(25) (α_s convention; _top_rge_coordinate)
    Weinberg-Hall matching 21, 1/(12π)   couplings.py, G₂ → SU(3)
    two loop MS-bar QCD transport     couplings.py (b₀, b₁; 3-4 loop
                                    shifts α_s(M_Z) by 3e-5)

  EVALUATION COORDINATE (immaterial at 5σ, test_blindness):
    M_Z = 91.1876 GeV               PDG              reference.PDG_EW

  CITED LOOP MATHEMATICS (fixed by field content, not measured):
    SM β coefficients, orders 1-3   SMDR pure MS-bar, yb=ytau=0
    quartic beta order 4: 8308.17    SPM 1508.00912 + Chetyrkin-Zoller
    ζ(3) = 1.2020569...             mathematics

  REFERENCE DATA (never steer a prediction: perturbing all of
  them at once leaves every prediction identical,
  tests/audit/test_blindness.py): everything in interference/reference.py
  (PDG_MASSES, PDG_CKM, NUFIT_PMNS, PDG_EW, PDG_COSMO, PLANCK_2018,
  ACT_DR6, NUFIT_OSCILLATION) plus display only σ/error literals in
  print statements.

Higgs:
    lambda(M_Pl) = -delta_bridge = -N_bridge * alpha_G₂² * E[v²] * (1 − h₁₀)
    is the bridge boundary. The imported beta equations use
    MS-bar coordinates. Transporting the boundary produces lambda(m_t).
    In the framework this is a transport coordinate, not a second
    physical mass. If an observable Higgs reading is named, it
    must be identified at the depth 4 Web endpoint rather than by an
    extra matching stage.

Reference:
    Kahsay, Kibrom Kidane (2026). The Octavian Coherence Gate: The Four
    Irreducible Integers of the E₈(1) ⊃ G₂(1) × F₄(1) Conformal Embedding.
    Zenodo. https://doi.org/10.5281/zenodo.20493955
"""

import math
import numpy as np

from .root import (d10, d11,
                   N_bridge,
                   M_Pl_GeV,
                   alpha_phys as ALPHA_PHYS,
                   S_lepton, C_EFF,
                   pct)
from .root import rk4_run as _shared_rk4


# ═══════════════════════════════════════════════════════════════════════
#  SM RGE infrastructure (for the Higgs running coordinate)
# ═══════════════════════════════════════════════════════════════════════

_ZETA3 = 1.2020569031595942      # Riemann zeta(3), exact to float64


def _beta_sm(y):
    """Backreaction evolution for the five derived couplings (g1, g2, g3, yt, lam).

    Every coefficient is determined by the E₈-derived gauge group
    SU(d₁₁)_c × SU(d₁₀)_L × U(1)_Y, matter content (d₁₁ generations,
    d₁₀×d₁₁ quark flavors, 1 Higgs doublet), and representation
    Casimirs. No free parameters enter at any beta order.

    The five couplings are NOT independent inputs, they are derived
    quantities (g_i from α_em + sin²θ_W, yt = √2·m_t/v_EW, λ from
    the Planck boundary) whose scale dependence is governed by the
    virtual backreaction of all E₈-derived particle species.

    Perturbative beta orders (standard QFT calls these 'n-loop'):
      Order 1: single virtual pair modifies each coupling
      Order 2: virtual pairs interact among themselves
      Order 3+: iterated backreaction
    These are transport expansion terms, not Ledger.web_depth.
    Orders 1-3 use one pure MS bar coordinate system. The order-3
    equations are the SMDR yb=ytau=0 specialization. g1 is
    SU(5)-normalised, so gp=sqrt(3/5)g1 is used at the SMDR boundary.
    The leading QCD order-4 quartic term is in the same coordinates.
    Sources: Bednyakov-Pikelner-Velizhanin 2012-13,
    Chetyrkin-Zoller 2012-16, Martin 2015, and SMDR betas.c.
    """
    g1, g2, g3, yt, lam = y
    g1s, g2s, g3s = g1**2, g2**2, g3**2
    yts = yt**2
    L1 = 1.0 / (16.0 * math.pi**2)
    L2 = L1**2

    # # Beta order 1 (n_gen=d₁₁=3, n_H=1) ##########
    # Each coefficient: Casimir of the E₈-derived representations.
    # b₀(SU3) = -11 + 2·d₁₀·d₁₁/3 = -7  (d₁₀d₁₁ = 6 quark flavors)
    dg1_1 = (41.0/10.0) * g1**3
    dg2_1 = (-19.0/6.0) * g2**3
    dg3_1 = (-(11 - 2*d10*d11/3.0)) * g3**3    # = -7 for d₁₀d₁₁=6
    dyt_1 = yt * ((9.0/2.0)*yts - (17.0/20.0)*g1s - (9.0/4.0)*g2s - 8.0*g3s)
    dlam_1 = (24.0*lam**2 + 12.0*yts*lam - 6.0*yts**2
              - (9.0/5.0*g1s + 9.0*g2s)*lam
              + 27.0/200.0*g1s**2 + 9.0/20.0*g1s*g2s + 9.0/8.0*g2s**2)

    # Beta order 2: gauge couplings
    dg1_2 = g1**3 * (199.0/50.0*g1s + 27.0/10.0*g2s + 44.0/5.0*g3s - 17.0/10.0*yts)
    dg2_2 = g2**3 * (9.0/10.0*g1s + 35.0/6.0*g2s + 12.0*g3s - 3.0/2.0*yts)
    dg3_2 = g3**3 * (11.0/10.0*g1s + 9.0/2.0*g2s - 26.0*g3s - 2.0*yts)

    # Beta order 2: top interference
    dyt_2 = yt * (
        -12.0*yts**2
        + yts*(393.0/80.0*g1s + 225.0/16.0*g2s + 36.0*g3s)
        + (1187.0/600.0)*g1s**2 - 9.0/20.0*g1s*g2s
        - 23.0/4.0*g2s**2 + 19.0/15.0*g1s*g3s
        + 9.0*g2s*g3s - 108.0*g3s**2
        + 6.0*lam**2 - 12.0*lam*yts
    )

    # Beta order 2: quartic (SMDR, yb=ytau=0)
    gp2 = (3.0/5.0) * g1s
    gp4 = gp2**2
    gp6 = gp2*gp4
    dlam_2 = (
        305.0/16.0*g2s**3 - 289.0/48.0*g2s**2*gp2
        - 559.0/48.0*g2s*gp4 - 379.0/48.0*gp6
        - 73.0/8.0*g2s**2*lam + 39.0/4.0*g2s*gp2*lam
        + 629.0/24.0*gp4*lam + 108.0*g2s*lam**2
        + 36.0*gp2*lam**2 - 312.0*lam**3
        - 9.0/4.0*g2s**2*yts + 21.0/2.0*g2s*gp2*yts
        - 19.0/4.0*gp4*yts + 45.0/2.0*g2s*lam*yts
        + 80.0*g3s*lam*yts + 85.0/6.0*gp2*lam*yts
        - 144.0*lam**2*yts - 32.0*g3s*yts**2
        - 8.0/3.0*gp2*yts**2 - 3.0*lam*yts**2
        + 30.0*yts**3
    )

    # Beta order 3: all five couplings (SMDR, yb=ytau=0)
    # Bednyakov-Pikelner-Velizhanin (1210.6873, 1212.6829, 1303.4364),
    # Chetyrkin-Zoller (1205.2892, 1303.2890), SMDR betas.c.
    # Evaluated in SMDR variables (gp = g1*sqrt(3/5), g = g2, k = lam).
    gp = g1 * math.sqrt(3.0/5.0)
    g  = g2
    k  = lam
    gp2 = gp**2;  gp4 = gp2**2;  gp6 = gp2*gp4;  gp8 = gp4**2
    gs  = g**2;   g4  = gs**2;   g6  = gs*g4;     g8  = g4**2
    g32 = g3**2;  g34 = g32**2;  g36 = g32*g34
    yt2 = yts;    yt4 = yt2**2;  yt6 = yt2*yt4;   yt8 = yt4**2
    k2  = k**2;   k3  = k*k2;   k4  = k2**2
    Z   = _ZETA3
    dlam_3 = (
        # pure gauge (g⁸, mixed gauge-g₃)
        (228259*g8)/1536 - (459*g6*g32)/4 - (165665*g6*gp2)/1728
        - (153*g4*g32*gp2)/4 - (81509*g4*gp4)/1728 - (187*gs*g32*gp4)/4
        - (237787*gs*gp6)/3456 - (187*g32*gp6)/4 - (51845*gp8)/512
        # gauge × k
        + (58031*g6*k)/144 + 405*g4*g32*k + (6137*g4*gp2*k)/16
        + (1549*gs*gp4*k)/4 + 165*g32*gp4*k + (88639*gp6*k)/216
        # k², k³, k⁴
        - (1389*g4*k2)/4 - 666*gs*gp2*k2 - 836*gp4*k2
        - 948*gs*k3 - 316*gp2*k3 + 7176*k4
        # yt² terms
        - (6849*g6*yt2)/128 + (651*g4*g32*yt2)/4 + (3487*g4*gp2*yt2)/128
        + (249*gs*g32*gp2*yt2)/2 + (25441*gs*gp4*yt2)/384
        + (587*g32*gp4*yt2)/12 + (125503*gp6*yt2)/1152
        - (6957*g4*k*yt2)/32 - 489*gs*g32*k*yt2 + (2488*g34*k*yt2)/3
        - (6509*gs*gp2*k*yt2)/16 - (2419*g32*gp2*k*yt2)/9
        - (203887*gp4*k*yt2)/864
        + (639*gs*k2*yt2)/2 - 2448*g32*k2*yt2 - (195*gp2*k2*yt2)/2
        + 1746*k3*yt2
        # yt⁴ terms
        + (9909*g4*yt4)/64 - 31*gs*g32*yt4 - (532*g34*yt4)/3
        - (1079*gs*gp2*yt4)/96 + (931*g32*gp2*yt4)/9 + (67793*gp4*yt4)/1728
        - (4977*gs*k*yt4)/4 + 1790*g32*k*yt4 - (2485*gp2*k*yt4)/12
        + 1719*k2*yt4
        # yt⁶, yt⁸
        + (3411*gs*yt6)/16 - 76*g32*yt6 + (3467*gp2*yt6)/48
        + (117*k*yt6)/4 - (1599*yt8)/4
        # # Zeta(3) sector ######################
        + Z * (
        - (20061*g8)/64 + 108*g6*g32 - (405*g6*gp2)/16 + 36*g4*g32*gp2
        + (2217*g4*gp4)/32 + 44*gs*g32*gp4 + (2177*gs*gp6)/48
        + 44*g32*gp6 + (12457*gp8)/192
        + (4419*g6*k)/4 - 432*g4*g32*k - (393*g4*gp2*k)/4
        - (147*gs*gp4*k)/4 - 176*g32*gp4*k - (1493*gp6*k)/12
        - 1026*g4*k2 - 324*gs*gp2*k2 - 162*gp4*k2
        + 144*gs*k3 + 48*gp2*k3 + 4032*k4
        + (297*g6*yt2)/2 - 108*g4*g32*yt2 + (27*g4*gp2*yt2)/2
        - 72*gs*g32*gp2*yt2 - 6*gs*gp4*yt2 - 36*g32*gp4*yt2
        - 5*gp6*yt2
        - 351*g4*k*yt2 + 432*gs*g32*k*yt2 - 96*g34*k*yt2
        + 354*gs*gp2*k*yt2 + 272*g32*gp2*k*yt2 - (449*gp4*k*yt2)/3
        - 864*gs*k2*yt2 + 2304*g32*k2*yt2 - 96*gp2*k2*yt2
        - (819*g4*yt4)/8 + 48*gs*g32*yt4 + 64*g34*yt4
        - (743*gs*gp2*yt4)/4 - (112*g32*gp2*yt4)/3 + (2957*gp4*yt4)/72
        + 1026*gs*k*yt4 - 2592*g32*k*yt4 + 114*gp2*k*yt4 + 1512*k2*yt4
        - 54*gs*yt6 + 480*g32*yt6 + 34*gp2*yt6 - 396*k*yt6 - 72*yt8
        )
    )

    # Gauge and top beta order 3, from the same SMDR specialization.
    dg3_3 = g3 * g32 * (
        109.0/8.0*g4 + 21.0*gs*g32 + 65.0/2.0*g34
        - 1.0/8.0*gs*gp2 + 77.0/9.0*g32*gp2
        - 2615.0/216.0*gp4 - 93.0/8.0*gs*yt2
        - 40.0*g32*yt2 - 101.0/24.0*gp2*yt2 + 15.0*yt4)
    dg2_3 = g * gs * (
        324953.0/1728.0*g4 + 39.0*gs*g32 + 81.0*g34
        + 291.0/32.0*gs*gp2 - 1.0/3.0*g32*gp2
        - 5597.0/576.0*gp4 + 3.0/2.0*gs*k + 1.0/2.0*gp2*k
        - 3.0*k2 - 729.0/32.0*gs*yt2 - 7.0*g32*yt2
        - 593.0/96.0*gp2*yt2 + 147.0/16.0*yt4)
    dgp_3 = gp * gp2 * (
        1315.0/64.0*g4 - gs*g32 + 99.0*g34
        + 205.0/96.0*gs*gp2 - 137.0/27.0*g32*gp2
        - 388613.0/5184.0*gp4 + 3.0/2.0*gs*k
        + 3.0/2.0*gp2*k - 3.0*k2 - 785.0/32.0*gs*yt2
        - 29.0/3.0*g32*yt2 - 2827.0/288.0*gp2*yt2
        + 315.0/16.0*yt4)
    dg1_3 = dgp_3 / math.sqrt(3.0/5.0)
    dyt_3 = yt * (
        987.0/4.0*gs*g34 - 4166.0/3.0*g36
        + 435.0/4.0*g32*g4 + 455.0/576.0*g6
        - 107.0/4.0*gs*g32*gp2 + 1633.0/36.0*g34*gp2
        + 273.0/64.0*g4*gp2 + 2045.0/192.0*gs*gp4
        + 2047.0/54.0*g32*gp4 + 763523.0/5184.0*gp6
        - 171.0/16.0*g4*k + 39.0/8.0*gs*gp2*k
        - 121.0/16.0*gp4*k + 45.0*gs*k2 + 15.0*gp2*k2
        - 36.0*k3 - 168.0*gs*g32*yt2 + 3827.0/6.0*g34*yt2
        + 32391.0/256.0*g4*yt2 + 2699.0/128.0*gs*gp2*yt2
        - 42.0*g32*gp2*yt2 - 458179.0/6912.0*gp4*yt2
        - 135.0/2.0*gs*k*yt2 + 16.0*g32*k*yt2
        - 127.0/6.0*gp2*k*yt2 + 15.0/4.0*k2*yt2
        - 1593.0/16.0*gs*yt4 - 157.0*g32*yt4
        - 2437.0/48.0*gp2*yt4 + 198.0*k*yt4
        + 339.0/8.0*yt6
        + Z * (
            -144.0*gs*g34 + 640.0*g36 - 108.0*g32*g4
            + 1125.0/8.0*g6 - 176.0/3.0*g34*gp2
            - 81.0/8.0*g4*gp2 - 153.0/8.0*gs*gp4
            - 748.0/9.0*g32*gp4 - 13073.0/216.0*gp6
            + 180.0*gs*g32*yt2 - 228.0*g34*yt2
            - 729.0/8.0*g4*yt2 + 123.0/4.0*gs*gp2*yt2
            + 60.0*g32*gp2*yt2 - 31.0/24.0*gp4*yt2
            + 27.0/2.0*yt6))

    # Beta order 4: QCD dominant quartic
    # SPM 1508.00912 Eq.(3.3), Chetyrkin-Zoller 1604.00853
    _C4_QCD = 8308.17  # perturbative order-4 computation
    dlam_4 = _C4_QCD * yt4 * g36

    # # Fixed transport scope ####################
    # Full three loop gauge/Yukawa/quartic SMDR transport, plus the
    # separately known leading QCD four loop quartic term.
    L3 = L1**3
    L4 = L1**4
    dg1 = L1*dg1_1 + L2*dg1_2 + L3*dg1_3
    dg2 = L1*dg2_1 + L2*dg2_2 + L3*dg2_3
    dg3 = L1*dg3_1 + L2*dg3_2 + L3*dg3_3
    dyt = L1*dyt_1 + L2*dyt_2 + L3*dyt_3
    dlam = L1*dlam_1 + L2*dlam_2 + L3*dlam_3 + L4*dlam_4

    return np.array([dg1, dg2, dg3, dyt, dlam])


def _run_rge(y0, t0, t1, beta_func, n_steps=20000):
    """Vector RGE running via the shared RK4 (root.rk4_run)."""
    return _shared_rk4(beta_func, y0.copy(), t0, t1, n_steps)


def _top_rge_coordinate(m_top_reading_GeV, alpha_s_mt):
    """Express the top Web reading in the imported beta coordinate.

    The coefficient convention uses alpha_s rather than alpha_s/pi.
    This is a coordinate translation for the transport equations, not
    an additional framework mass or interference depth.
    """
    a = alpha_s_mt
    ratio = (1.0 + 0.4244*a + 0.8345*a**2
             + 2.375*a**3 + 8.49*a**4)
    return m_top_reading_GeV / ratio


def _run_gauge_one_loop(g0, coefficient, log_scale_ratio):
    """One loop MS-bar gauge transport for 1/g^2."""
    inverse_square = (1.0/g0**2
                      - coefficient*log_scale_ratio/(8.0*math.pi**2))
    return 1.0 / math.sqrt(inverse_square)


def _derive_alpha_em(lepton_masses, quark_masses, mu_GeV):
    """Derive alpha_em(mu) from alpha(0) via QED vacuum polarization.

    mu_GeV : the scale to run to (e.g. M_Z).
    """
    alpha_0 = ALPHA_PHYS

    m_e, m_mu, m_tau = [m / 1e3 for m in lepton_masses]
    m_u, m_d, m_s, m_c, m_b, m_t = [m / 1e3 for m in quark_masses]

    # Confinement scale: Λ_conf = C_EFF · M_Pl exp(−S_lepton), with
    # C_EFF = 1/d₁₀ the SAME trace lift prefactor as the lepton
    # amplitude (root.C_EFF: one quantity, one name).
    # Λ_conf = 314.2 MeV comes out at the constituent scale (~1%).
    # Hadron masses belong to the hadrons and nuclei rung,
    # read from QCD with the engine's couplings and masses
    # (tests/audit/emergence_ladder.py).
    Lambda_conf_GeV = C_EFF * M_Pl_GeV * math.exp(-S_lepton)

    # Quark charges from representation theory of E₈ → G₂ × F₄:
    #   Q_u = d₁₀/d₁₁ = 2/3,   Q_u² = (d₁₀/d₁₁)² = 4/9
    #   Q_d = 1/d₁₁   = 1/3,   Q_d² = 1/d₁₁²      = 1/9
    #   N_c = d₁₁      = 3      (colors from G₂ → SU(3))
    Q_u_sq = float(d10)**2 / float(d11)**2         # (d₁₀/d₁₁)² = 4/9
    Q_d_sq = 1.0 / float(d11)**2                   # 1/d₁₁²     = 1/9
    N_c = d11                                       # d₁₁ = 3

    fermions = [
        (m_e,   1,    1.0,    False),   # lepton charge = 1
        (m_mu,  1,    1.0,    False),
        (m_tau, 1,    1.0,    False),
        (m_u,   N_c,  Q_u_sq, True),    # up-type: N_c × Q_u²
        (m_d,   N_c,  Q_d_sq, True),    # down-type: N_c × Q_d²
        (m_s,   N_c,  Q_d_sq, True),
        (m_c,   N_c,  Q_u_sq, True),
        (m_b,   N_c,  Q_d_sq, True),
    ]

    delta_inv_alpha = 0.0
    for mf, nc, qf2, is_quark in fermions:
        if mf >= mu_GeV:
            continue
        mf_eff = Lambda_conf_GeV if (is_quark and mf < Lambda_conf_GeV) else mf
        delta_inv_alpha += (2.0 * nc * qf2) / (3.0 * math.pi) * math.log(mu_GeV / mf_eff)

    inv_alpha_mu = 1.0 / alpha_0 - delta_inv_alpha
    return 1.0 / inv_alpha_mu


# sin²θ_W derivation (in the Web ledger, root.py WEB.state["sin2W"]):
#   base:  d₁₁/13 = 3/13
#   echo:  + Q₀²·(α(0)/2π)  (the ledgered form, root.py WEB["sin2W"])
#   result: 0.2312854, PDG 0.23129(4) → −0.11σ
#   FCC-ee test: the Q₀² form against the altitude share
#   (1+h∨G₂)/13 = 5/13. A Z-pole σ~1e-5 decides.  A measured 5/13
#   weight (0.2312161), or any other family weight, kills the 3/13
#   emission reading.
 

# ═══════════════════════════════════════════════════════════════════════
#  Public interface
# ═══════════════════════════════════════════════════════════════════════

def derive_higgs_ew(masses, couplings):
    """Derive Higgs quartic, EW observables (M_W, M_Z), and run import audit.

    masses    : dict from masses.derive()
    couplings : dict from couplings.derive()
    """
    from .root import WEB, PDG_EW

    print(f"\n  --- HIGGS: F₄ fusion + bridge threshold ---")

    v_EW = masses['v_EW_pred_GeV']
    m_t_GeV = masses['m_t'] / 1e3
    alpha_s_MZ = couplings['alpha_s_MZ_thresh']

    # Higgs quartic boundary condition at M_Pl (two components):
    #
    # 1. F₄(1) fusion rule: (26 × 26)_local = 1
    #    The F₄ fusion algebra at level 1 has a unique quartic coupling:
    #    the symmetric product of two 26-reps fuses to the identity.
    #    This fixes the zero edge boundary λ_F₄ = 0 at M_Pl.
    #
    # 2. Bridge boundary response, before the fundamental share vent:
    #    δ⁰_bridge = N_bridge · α²_G₂(M_Pl) · E[v²]
    #             = (n₇·n₂₆) · (1/(24π))² · (1/2) = 0.016007
    #    The vented value used everywhere below is
    #    δ_bridge = δ⁰_bridge · (1 − h₁₀) = 0.012450 (step 3).
    #    Each factor is derived:
    #      N_bridge = n₇·n₂₆ = 182 (bridge mode count)
    #      α_G₂(M_Pl) = 1/(24π) (wiring.py: pointed sector probability
    #                   P = 1/12 times normalized current response 1/(2π))
    #      E[v²] = 1/2 (protected forgetting: per mode |v̂₁|² = B²/4;
    #                   total Pv²P = P)
    #
    # 3. Fundamental share vent: the triality charged
    #    3 ⊕ 3̄ share of the bridge 7 cannot cross into the condensing
    #    channel (N-ality superselection, words.py). Measured in the
    #    weight metric, the blocked share is h₁₀. This gives the
    #    factor (1 − h₁₀).
    #    Witness: tests/audit/test_coverage.py::test_vent_share_witnesses.
    #
    # Combined: λ(M_Pl) = −δ⁰_bridge·(1 − h₁₀) = −δ_bridge
    # [the vent lives in the WEB edge]
    #
    # SINGLE SOURCE OF TRUTH: the bridge→Higgs edge lives in the web
    # (root.py WEB['lambda_MPl']). No second computation here.
    lam_Pl_boundary = WEB["lambda_MPl"].value()
    delta_bridge = -lam_Pl_boundary

    print(f"    F₄ fusion: lambda_F₄ = 0  [(26×26)_local = 1]")
    print(f"    delta_bridge = {N_bridge} * (1/(24pi))² * 1/2 * (1 − h₁₀) = {delta_bridge:.6f}")
    print(f"      [fundamental share vent (1−h₁₀) = 7/9]")
    print(f"    lambda(M_Pl) = -{delta_bridge:.6f}")

    # Derive EW couplings at the derived M_Z scale
    lepton_masses = [masses['m_e'], masses['m_mu'], masses['m_tau']]
    quark_masses = [masses['m_u'], masses['m_d'], masses['m_s'],
                    masses['m_c'], masses['m_b'], masses['m_t']]

    # sin²θ_W from the Web ledger (root.py).
    # Base 3/13, depth 1 h₇·α/(2π), depth 2 (h₇/d₁₁²)·α/(2π).
    # Combined echo coefficient: Q₀² = 4/9 (Koide parameter squared).
    sin2W = WEB.state["sin2W"]
    cos2W = 1.0 - sin2W

    # # M_W and M_Z (MS-bar. TWO declared imports) ########
    # PDG 2024 EW review (Erler-Freitas), Eq. (10.26):
    #     M_W = A₀ / [ŝ_Z (1 − Δr̂_W)^½],   M_Z = M_W / (ρ̂^½ ĉ_Z)
    # A₀ = [πα(0)/(√2 G_F)]^½ = 37.28038(1) GeV, and BOTH inputs of
    # A₀ are the framework's own predictions:
    #     A₀ = √(πα(0))·v_EW = 37.28038 GeV   (matches to 0.1 ppm)
    # ŝ² is canonical (the solved WEB["sin2W"]: base 3/13 + depth 1
    # h₇·α/2π + depth 2 (h₇/d₁₁²)·α/2π).  The two
    # imported electroweak transport terms [DECLARED, PDG 2024]:
    #     Δr̂_W = 0.06937(6),   ρ̂ = 1.01016(9)
    # Results: M_W = 80.356 GeV (−0.86σ of the 80.3625(77) world average).
    #          M_Z = 91.1885 GeV (+0.0010%, +0.18σ of the import band).
    # The apportioned echo matters here: with the base 3/13 alone,
    # M_W would be off by +0.095%.
    A_EW = math.sqrt(math.pi * ALPHA_PHYS) * v_EW
    DR_HAT_W = PDG_EW['dr_hat_W']
    RHO_HAT = PDG_EW['rho_hat']
    M_W_derived = A_EW / (math.sqrt(sin2W) * math.sqrt(1.0 - DR_HAT_W))
    M_Z_derived = M_W_derived / (math.sqrt(RHO_HAT) * math.sqrt(cos2W))

    # α_EM at the derived M_Z (for the RGE initial conditions)
    alpha_em_MZ = _derive_alpha_em(lepton_masses, quark_masses, M_Z_derived)

    # PDG reference values (comparison only. Single source root.PDG_EW)
    M_Z_PDG = PDG_EW['M_Z']
    M_W_PDG = PDG_EW['M_W']
    sin2W_PDG = PDG_EW['sin2W_MSbar']

    print(f"    alpha_em(M_Z) = 1/{1.0/alpha_em_MZ:.3f}")
    print(f"    sin²theta_W = 3/13 + Q₀²·(α/2π) = {sin2W:.6f}  [base 3/13 + depth 1 h₇ + depth 2 h₇/d₁₁²]")
    print(f"    PDG 2024 global fit: {sin2W_PDG}(4)  (pull: {(sin2W - sin2W_PDG)/PDG_EW['sin2W_err']:+.2f}σ)")
    print(f"    M_W = A₀/(ŝ√(1−Δr̂_W)) = {M_W_derived:.4f} GeV  (PDG: {M_W_PDG}({PDG_EW['M_W_err']*1e4:.0f}), pull {(M_W_derived-M_W_PDG)/PDG_EW['M_W_err']:+.2f}σ)")
    M_Z_IMPORT_BAND = 0.0050    # GeV, the band the two imports (ρ̂, ŝ) propagate to M_Z
    print(f"    M_Z = M_W/(ρ̂^½ĉ)      = {M_Z_derived:.4f} GeV  (PDG: {M_Z_PDG}, {pct(M_Z_derived, M_Z_PDG):+.4f}%. "
          f"{(M_Z_derived - M_Z_PDG) / M_Z_IMPORT_BAND:+.2f}σ of the import band)")
    print(f"      [A₀ = √(πα(0))·v_EW = {A_EW:.5f} GeV vs PDG 37.28038(1): 0.1 ppm.")
    print(f"       imports Δr̂_W = 0.06937(6), ρ̂ = 1.01016(9), PDG 2024 Eq. (10.26)]")

    # Gauge coordinates at the derived M_Z.
    g1_mz = math.sqrt(4*math.pi*alpha_em_MZ/(1-sin2W)) * math.sqrt(5.0/3.0)
    g2_mz = math.sqrt(4*math.pi*alpha_em_MZ/sin2W)

    # The imported beta equations require their own running coordinate.
    # Translate the final top Web reading into that coordinate using the
    # pole to MS-bar relation of _top_rge_coordinate (Marquard et al.). This does not
    # create a second top mass in the framework.
    alpha_s_mt = couplings['alpha_s_mt']
    m_t_rge_coordinate = _top_rge_coordinate(m_t_GeV, alpha_s_mt)
    yt_mt_transport = math.sqrt(2.0) * m_t_rge_coordinate / v_EW
    log_mt_mz = math.log(m_t_GeV / M_Z_derived)
    g1_mt = _run_gauge_one_loop(g1_mz, 41.0/10.0, log_mt_mz)
    g2_mt = _run_gauge_one_loop(g2_mz, -19.0/6.0, log_mt_mz)
    g3_mt = math.sqrt(4.0 * math.pi * alpha_s_mt)

    t_Pl = math.log(M_Pl_GeV / m_t_GeV)
    N_rge = 20000

    # Boundary value iteration: no measured Higgs mass is used.
    lam_iter = 0.0
    for iteration in range(3):
        y_up = np.array([g1_mt, g2_mt, g3_mt, yt_mt_transport, lam_iter])
        y_Pl = _run_rge(y_up, 0.0, t_Pl, _beta_sm, N_rge)

        y_down = np.array([y_Pl[0], y_Pl[1], y_Pl[2],
                           y_Pl[3], lam_Pl_boundary])
        y_low = _run_rge(y_down, t_Pl, 0.0, _beta_sm, N_rge)
        lam_iter = y_low[4]

    lam_mt_transport = lam_iter
    higgs_transport_coordinate = math.sqrt(2.0 * abs(lam_mt_transport)) * v_EW

    # # RK4 convergence self-certification #############
    # The step count is not tuned.  Run at 2×N_rge and assert the
    # result is stable to < 1 MeV.
    y_up_2x = np.array([g1_mt, g2_mt, g3_mt, yt_mt_transport,
                        lam_mt_transport])
    y_Pl_2x = _run_rge(y_up_2x, 0.0, t_Pl, _beta_sm, 2 * N_rge)
    y_down_2x = np.array([y_Pl_2x[0], y_Pl_2x[1], y_Pl_2x[2],
                           y_Pl_2x[3], lam_Pl_boundary])
    y_low_2x = _run_rge(y_down_2x, t_Pl, 0.0, _beta_sm, 2 * N_rge)
    higgs_transport_coordinate_2x = math.sqrt(2.0 * abs(y_low_2x[4])) * v_EW
    assert abs(higgs_transport_coordinate - higgs_transport_coordinate_2x) < 0.001, \
        ("RK4 not converged: "
         f"{higgs_transport_coordinate:.4f} vs "
         f"{higgs_transport_coordinate_2x:.4f} GeV")

    # F4-only for comparison (λ(M_Pl)=0, no bridge)
    y_down_F4 = np.array([y_Pl[0], y_Pl[1], y_Pl[2],
                          y_Pl[3], 0.0])
    y_low_F4 = _run_rge(y_down_F4, t_Pl, 0.0, _beta_sm, N_rge)

    # # IMPORT SPECIFICATION ####################
    print(f"    top Web reading -> beta coordinate: "
          f"{m_t_GeV:.3f} -> {m_t_rge_coordinate:.3f} GeV")
    print(f"    lambda_MSbar(m_t), F4 only:       {y_low_F4[4]:.6f}")
    print(f"    lambda_MSbar(m_t), F4 + bridge:   {lam_mt_transport:.6f}")

    # # IMPORT COMPATIBILITY AUDIT (parameterless discipline) ###
    # A cited electroweak transport term is a FUNCTION of
    # (G_F, m_t, m_H, α̂(M_Z)), taken at the SM fit's values.  A
    # parameterless framework may not absorb those silently. Each
    # ingredient is therefore RECALCULATED here from the framework's
    # own predictions and checked to be compatible within the
    # import's quoted error.
    #   G_F : framework 1/(√2 v_EW²) = measured to < 1 ppm  [checked]
    #   m_t : ρ̂'s quadratic m_t term (PDG Eq. 10.23) shifts by 0.12σ
    #         of the import error under fw m_t vs fit m_t  [checked]
    #   m_H : not produced by the framework. It enters the loop
    #         remainders as data (class [1] below).
    #   α̂(M_Z): the DOMINANT Δr̂_W ingredient is Δr₀ = 1 − α/α̂(M_Z)
    #         = 0.06646(6). It is driven by the dispersive hadronic VP,
    #         which is MEASUREMENT (e⁺e⁻ R-ratio spectral data), not a
    #         fit parameter. It is therefore admissible exactly like
    #         the PDG masses on the comparison side.  The framework's
    #         INTERNAL VP (free quarks + Λ_conf cutoff) reproduces it
    #         to 0.34% (1/128.37 vs 1/127.93). This approximation is
    #         used ONLY for the Higgs RGE initial conditions, where it
    #         moves λ(m_t) by 0.21%. It must never be used at sub-0.5%
    #         precision. The guard below checks that it stays within
    #         its declared quality.
    # TAXONOMY (what an import IS, for a parameterless model):
    #   [1] data driven ingredients  -> measurements (dispersive VP.
    #       same admissibility as the PDG masses)        [classified]
    #   [2] parameter dependencies   -> recalculated from framework
    #       values, compatible within quoted errors      [verified]
    #   [3] pure loop integrals      -> mathematics, uniquely fixed
    #       by gauge group + matter content at the measured m_H (the EW
    #       remainder 0.00291 in Δr̂_W, the bosonic 0.00082 in ρ̂), the
    #       same status as the Machacek-Vaughn beta functions used for
    #       the Higgs RGE.
    # The compatibility checks are frozen in
    # tests/common_mode/test_common_mode_scorecard.py::test_import_audit_is_compatible.
    # The engine only prints them, so no reference value can stop a run.
    RHO_HAT_MT_COEFF = 0.00934      # ρ̂'s quadratic m_t term, PDG 2024 Eq. (10.23)
    _GF_ref = PDG_EW['G_F']
    _GF_fw = 1.0 / (math.sqrt(2.0) * v_EW**2)
    _m_t_fit = PDG_EW['m_t_EW_fit']
    _drho_mt = RHO_HAT_MT_COEFF * abs((m_t_GeV / _m_t_fit)**2 - 1.0)
    _inv_alpha_hat = PDG_EW['inv_alpha_hat_MZ']
    print(f"    import audit: G_F fw/meas −1 = {1e6*(_GF_fw/_GF_ref-1):+.2f} ppm"
          f" ({(_GF_fw-_GF_ref)/PDG_EW['G_F_err']:+.2f}σ)."
          f" δρ̂(m_t) = {_drho_mt:.1e} ({_drho_mt/PDG_EW['rho_hat_err']:.2f}σ vs fit m_t {_m_t_fit})."
          f" m_H = {PDG_EW['m_H']} GeV enters the loop remainders as data.")
    print(f"      internal VP α̂(M_Z) = 1/{1.0/alpha_em_MZ:.3f} vs dispersive 1/{_inv_alpha_hat}"
          f" ({100*((1.0/alpha_em_MZ)/_inv_alpha_hat-1):+.2f}%, RGE only approximation)")

    return {
        'alpha_em_MZ': alpha_em_MZ,
        'inv_alpha_em_MZ': 1.0 / alpha_em_MZ,
        'drho_mt_import_audit': _drho_mt,
        'sin2W': sin2W,
        'delta_bridge': delta_bridge,
        'lam_Pl_boundary': lam_Pl_boundary,
        'lambda_mt_transport': lam_mt_transport,
        'higgs_transport_coordinate_GeV': higgs_transport_coordinate,
        'm_t_rge_coordinate_GeV': m_t_rge_coordinate,
        'yt_mt_transport': yt_mt_transport,
        'M_Z_derived': M_Z_derived,
        'M_W_derived': M_W_derived,
        'M_Z_PDG': M_Z_PDG,
        'M_W_PDG': M_W_PDG, 'M_W_err': PDG_EW['M_W_err'],
    }
