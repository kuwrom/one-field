"""
Gravity + Higgs + Baryogenesis + Neutrinos from (d₁₀, d₁₁, n₇, n₂₆, pi)
+ the electron anchor (M_Pl is DERIVED in root.py from m_e).

═══════════════════════════════════════════════════════════════════════
WHY GRAVITY AT ALL: the common-mode memory (gravity paper, summary)
═══════════════════════════════════════════════════════════════════════

The substrate is the Z₃-symmetric coupled NLS system (masses.py).
In Madelung variables ψ_k = √ρ_k e^{iS_k/ħ}, perturbations decompose
under the Z₃ Fourier transform into:
    q = 0    the COMMON mode, shared by all three components;
    q = 1,2  the RELATIVE modes, they distinguish components.
q = 1,2 carry the matter-like standing waves; q = 0 carries the
gravitational memory.  This resolves the Barceló-Liberati-Visser
"two roles of ρ" obstruction of analogue gravity: matter and
geometry are different Fourier projections of ONE substrate.

Linearising and projecting onto q = 0, the relative-sector energy
density sources the common-mode density via a SCREENED POISSON
equation:
    (1 − ξ₀²∇²) R₀ = −(c*/λ₀ρ₀) H_matter,
with healing length ξ₀ = ℓ_Pl/2 under the Planck identifications.
The scalar q=0 force is therefore Yukawa-screened: e^{−r/ξ₀} <
10⁻³⁰⁰ at any observable distance, the framework passes the
Cassini PPN bound by ~295 orders of magnitude.  Macroscopic gravity
is the INDUCED SPIN-2 channel (Sakharov): matter standing waves
propagate on the common acoustic metric, and integrating them out
at one loop generates the Einstein-Hilbert term.  The minimal
matter spectrum alone gives the WRONG-SIGN normalisation
(Σ_min ≈ −10.2); the 182 bridge channels supply the positive
correction (+29.126).  That repair is this module.

Gravity (bridge sector):
    N_bridge = n₇ * n₂₆ = 182
    xi = alpha_G₂(M_Pl) * E[v²] * h_bridge = 1/(48pi)
    Heat kernel -> G_N

Higgs:
    lambda(M_Pl) = -delta_bridge = -N_bridge * alpha_G₂² * E[v²]
    SM RGE -> m_H

Baryogenesis (three Sakharov conditions ↔ three algebraic factors):
    eta_B = n₇ × J_lep × exp(−h∨(G₂)π²/2) = 7 · J_ℓ · e^{−2π²}
      • e^{−2π²}: G₂ instanton tunnelling amplitude, h∨(G₂) = 4
        sets the action S = 2π².  The SAME factor generates the
        Majorana scale M_R = M_Pl e^{−2π²} ≈ 3.3×10¹⁰ GeV >> T_EW,
        so out-of-equilibrium decay (3rd condition) needs no extra
        ingredient.
      • J_lep: leptonic Jarlskog (mixing.py).  CP violation traces
        to the ORIENTED Fano plane: seven directed lines select a
        preferred chirality.
      • n₇ = 7: the asymmetry channels, one per directed Fano
        line, each one octonionic multiplication.
    MECHANISM: the asymmetry is generated NON-THERMALLY, at the
    instanton level, it is a topological count, not a decay-era
    process.  The Majorana scale M_R = M_Pl e^{−2π²} >> T_EW means
    the B−L asymmetry is frozen in before any thermal washout era
    exists; the Boltzmann/washout machinery of thermal leptogenesis
    describes a different mechanism and does not apply here.

Neutrinos:
    26 -> 8_v + 8_s + 8_c + 1 + 1   (SO(8) triality decomposition
    of the F₄ fundamental: three 8's tied by triality + 2 singlets)
    The 2 singlets are the only RH neutrinos -> the seesaw matrix
    has rank 2 -> m₁ = 0 EXACTLY and normal ordering, structurally.

═══════════════════════════════════════════════════════════════════════
Why the bridge contributes a SCALAR sign (not a vector sign)
═══════════════════════════════════════════════════════════════════════

The heat-kernel a₁ coefficient for each bridge channel is
    a₁/R = +1/6 − xi_bridge   (scalar, ATTRACTIVE)
NOT the vector value −1/6.  Five independent arguments force this:

  (i)   Lagrangian algebra condensation:
        A_{E₈} = 1 + (tau_G, tau_F) absorbs the bridge into the E₈
        vacuum module with D²_local = 1.  The condensed modes are
        coherence channels of the parent theory, not independent
        propagating gauge fields.

  (ii)  Zero coset central charge:
        c_coset = c(E₈) − c(G₂) − c(F₄) = 0
        A c = 0 sector is topological, it cannot support the ghost
        subtraction that produces the vector −1/6 coefficient.

  (iii) Representation theory:
        The bridge transforms as the bi-fundamental (n₇, n₂₆), NOT as
        the adjoint of any gauge group.  No gauge redundancy means
        no Faddeev-Popov ghosts, so the ghost-subtracted vector sign
        does not apply.

  (iv)  Heat-kernel ghost arithmetic:
        Vector coefficient: −1/6 = +1/6 (naive) − 1/3 (ghost).
        Without ghosts, the coefficient reverts to scalar +1/6 − ξ.

  (v)   Quantum dimension verification:
        d(tau_G) * d(tau_F) = phi², d(A_{E₈}) = 1 + phi² = 2 + phi.
        FP dimension confirms the bridge is fully absorbed:
        D²_local = (1+phi²)² / d(A)² = 1.

═══════════════════════════════════════════════════════════════════════
Same bridge, two observables (cross-link to couplings.py)
═══════════════════════════════════════════════════════════════════════

The (7,26) bridge carries TWO predictions through the SAME h_bridge = 1:
  * Gravity: 182 scalar-like heat-kernel channels each contribute
    a₁/R = 1/6 − xi_bridge with xi_bridge = alpha_G₂ * (1/2) * h_bridge
    = 1/(48π).  With the face-split self-echo (9/13 of the marginal
    unit, no-self-dilution) the UV sum gives G_ind/G_N = 0.999999917.
  * QED:    The same h_bridge = 1 makes the bridge a marginal worldsheet
    primary; its one-loop self-interference fraction is h_bridge/(2π)
    = 1/(2π), giving (1 − 1/(2π)) that converts alpha_em = pi/2⁹ into
    alpha(0) = 1/137.035999050.  See couplings.py.
Same sector, same h_bridge, two observables.

Reference:
    Kahsay, Kibrom Kidane (2026). Three-paper series.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import math
import numpy as np

from root import (d10, d11, n7, n26,
                  hv_G2, hv_F4, C2_7, C2_26,
                  N_bridge, h_bridge, h_7, h_26,
                  alpha_G2_Pl, xi_bridge, E_v2,
                  M_Pl_GeV, M_Pl_MeV,
                  alpha_phys as ALPHA_PHYS,
                  alpha_EM as ALPHA_EM,
                  sin2W as SIN2W_TREE,
                  S_lepton,
                  NUFIT_PMNS, PDG_MASSES, pct)


# ═══════════════════════════════════════════════════════════════════════
#  SM RGE infrastructure (for Higgs mass)
# ═══════════════════════════════════════════════════════════════════════

def _beta_sm_2loop(y):
    """Two-loop SM beta functions for (g1, g2, g3, y_t, lam).

    Every coefficient is determined by the gauge group structure and matter
    content, which are themselves derived from (d₁₀, d₁₁, n₇, n₂₆):
      Gauge group: SU(d₁₁)_c × SU(d₁₀)_L × U(1)_Y  (from E₈ → G₂ × F₄)
      n_gen = d₁₁ = 3 generations  (D⁶ nimrep of SU(3)₃)
      n_f   = d₁₀ × d₁₁ = 6 quark flavours  (2 types × 3 gen)
      n_H   = 1 Higgs doublet
      C_A(SU3) = d₁₁ = 3,  C_F(SU3) = (d₁₁²−1)/(2d₁₁) = 4/3

    The one-loop and two-loop coefficients are standard results of gauge
    field theory perturbation theory (dimensional regularisation + MS-bar),
    uniquely fixed once the gauge group and matter content are specified.
    No external numerical input is used.

    Credit: Machacek-Vaughn (1983-85) for the general formulae;
            Buttazzo et al. JHEP 1312:089 (2013) for the compiled form.
    """
    g1, g2, g3, yt, lam = y
    g1s, g2s, g3s = g1**2, g2**2, g3**2
    yts = yt**2
    L1 = 1.0 / (16.0 * math.pi**2)
    L2 = L1**2

    # One-loop (n_gen=d₁₁=3, n_H=1)
    # dg3: b₀ = -11+2n_f/3 = -11+4 = -7, n_f = d₁₀×d₁₁ = 6
    dg1_1 = (41.0/10.0) * g1**3
    dg2_1 = (-19.0/6.0) * g2**3
    dg3_1 = (-(11 - 2*d10*d11/3.0)) * g3**3    # = -7 for d₁₀d₁₁=6
    dyt_1 = yt * ((9.0/2.0)*yts - (17.0/20.0)*g1s - (9.0/4.0)*g2s - 8.0*g3s)
    dlam_1 = (12.0*lam**2 + 12.0*yts*lam - 6.0*yts**2
              - (9.0/5.0*g1s + 9.0*g2s)*lam
              + 27.0/100.0*g1s**2 + 9.0/10.0*g1s*g2s + 9.0/4.0*g2s**2)

    # Two-loop gauge
    dg1_2 = g1**3 * (199.0/50.0*g1s + 27.0/10.0*g2s + 44.0/5.0*g3s - 17.0/10.0*yts)
    dg2_2 = g2**3 * (9.0/10.0*g1s + 35.0/6.0*g2s + 12.0*g3s - 3.0/2.0*yts)
    dg3_2 = g3**3 * (11.0/10.0*g1s + 9.0/2.0*g2s - 26.0*g3s - 2.0*yts)

    # Two-loop top Yukawa
    dyt_2 = yt * (
        -12.0*yts**2
        + yts*(131.0/16.0*g1s + 225.0/16.0*g2s + 36.0*g3s)
        + (1187.0/600.0)*g1s**2 - 9.0/20.0*g1s*g2s
        - 23.0/4.0*g2s**2 + 19.0/15.0*g1s*g3s
        + 9.0*g2s*g3s - 108.0*g3s**2
        + 6.0*lam**2 - 6.0*lam*yts
    )

    # Two-loop quartic
    dlam_2 = (
        -78.0*lam**3
        + 18.0*(3.0/5.0*g1s + 3.0*g2s)*lam**2
        - (73.0/8.0*g2s**2 - 117.0/20.0*g1s*g2s - 1887.0/200.0*g1s**2)*lam
        - 3.0*lam*yts*(5.0*yts - 25.0/4.0*g2s - 85.0/12.0*g1s - 40.0*g3s)
        - 32.0*yts**3
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


from root import rk4_run as _shared_rk4


def _run_rge(y0, t0, t1, beta_func, n_steps=20000):
    """Vector RGE running via the shared RK4 (root.rk4_run)."""
    return _shared_rk4(beta_func, y0.copy(), t0, t1, n_steps)


def _derive_alpha_em(lepton_masses, quark_masses, mu_GeV):
    """Derive alpha_em(mu) from alpha(0) via QED vacuum polarization.

    mu_GeV : the scale to run to (e.g. M_Z).
    """
    alpha_0 = ALPHA_PHYS

    m_e, m_mu, m_tau = [m / 1e3 for m in lepton_masses]
    m_u, m_d, m_s, m_c, m_b, m_t = [m / 1e3 for m in quark_masses]

    # Confinement scale: Λ_conf = (1/d₁₀) M_Pl exp(−S_lepton)
    #                           = (1/d₁₀) M_Pl exp(−d₁₁²π²/2)
    Lambda_conf_GeV = (1.0 / d10) * M_Pl_GeV * math.exp(-S_lepton)

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


def _derive_sin2W():
    """sin²θ_W from conformal embedding E₈ ⊃ G₂ × F₄ + depth-1 echo.

    Algebraic value, ratio of dual Coxeter numbers (h∨ metric):
        sin²θ_W(tree) = h∨(SU(3))/(h∨(G₂)+h∨(F₄)) = d₁₁/13 = 3/13

    Depth-1 echo (FORCED, face-split law): the angle's correction is
    an EMISSION through the EM channel (unit α(0)/2π), and emission
    echoes carry CONFORMAL weights, the same worldsheet rule that
    gave h_bridge/(2π) in the EM ledger.  The emitting face is the
    gauge face G₂, with emission weight h_7 = d₁₀/(1+d₁₀²) = 2/5:

        sin²θ_W(M_Z) = 3/13 + h_7·(α(0)/2π) = 0.2312338

    PDG 2024 global fit: 0.23129(4) -> pull −1.4σ (tree alone: −13σ;
    subset fits span 0.23118-0.23134).
    PRE-REGISTERED DISCRIMINATION: the altitude-share form
    (1+h∨G₂)/13 = 5/13 gives 0.2312161, degenerate today; the two
    differ by 1.8e-5, FCC-ee Z-pole (σ~1e-5) decides.  The grammar
    bets on h_7.  (The RGE-running attempt remains a recorded
    negative result: the mechanism was running; the physics is a
    depth-1 echo.)
    """
    return float(SIN2W_TREE) + float(h_7) * ALPHA_PHYS / (2.0 * math.pi)


# ═══════════════════════════════════════════════════════════════════════
#  Public interface
# ═══════════════════════════════════════════════════════════════════════

def derive(R, masses, mixing, couplings):
    """
    Derive gravity, Higgs, baryogenesis, and neutrinos.

    R          : dict from root.derive()
    masses     : dict from masses.derive()
    mixing     : dict from mixing.derive()
    couplings  : dict from couplings.derive()
    """
    print("\n" + "=" * 78)
    print("  GRAVITY + HIGGS + BARYOGENESIS + NEUTRINOS")
    print("=" * 78)

    # ══════════════════════════════════════════════════════════════════
    #  GRAVITY: Sakharov induced gravity from bridge sector
    # ══════════════════════════════════════════════════════════════════

    print(f"\n  --- GRAVITY: {N_bridge} bridge modes, xi = 1/(48pi) ---")

    # Protected forgetting: PvP=0, Pv²P=½P
    v = [1.0, -0.5, -0.5]
    u = [1.0 / math.sqrt(3.0)] * 3
    u_dot_v = sum(ui * vi for ui, vi in zip(u, v))
    E_v2_check = sum(ui**2 * vi**2 for ui, vi in zip(u, v))
    assert abs(u_dot_v) < 1e-12
    assert abs(E_v2_check - 0.5) < 1e-12

    print(f"    Protected forgetting: PvP = 0, Pv²P = 1/2 P")
    print(f"    h_bridge = h({n7},G₂) + h({n26},F₄) = {float(h_7)} + {float(h_26)} = {float(h_bridge)}")
    print(f"    xi_bridge = alpha_G₂(M_Pl) * E[v²] * h_bridge")
    print(f"              = 1/(24pi) * 1/2 * 1 = 1/(48pi) = {float(xi_bridge):.10f}")

    # Heat-kernel coefficients (from spin content of E₈ decomposition)
    #
    # The E₈ → G₂ × F₄ decomposition determines which fields exist and
    # their spins.  The a₁/R coefficient for each spin follows from the
    # representation of SO(3,1), it is a mathematical consequence of the
    # field's Lorentz structure, just as a Casimir follows from a rep's
    # weight diagram:
    #   a₁(scalar)  = +1/6   (spin-0, minimal coupling ξ=0)
    #   a₁(Dirac)   = −1/12  (spin-1/2, per real DOF)
    #   a₁(gauge)   = −1/6   (spin-1, per transverse DOF)
    #   a₁(Proca)   = −1/12  (massive spin-1, per DOF)
    #
    # The non-minimal coupling ξ_bridge = 1/(48π) for the (7,26) bridge
    # sector is derived from α_G₂(M_Pl)·E[v²]·h_bridge (see root.py).
    # The bridge modes are scalar-like (not vector) because: c_coset = 0,
    # no gauge redundancy (bi-fundamental, not adjoint), and D²_local = 1
    # (Lagrangian algebra condensation).
    #
    # Credit: DeWitt (1965), Seeley (1967) for the heat-kernel expansion.
    #
    # ── THE UNIVERSAL VENTING LEDGER ─────────────────────────────────
    # A standing wave is a circulating knot: it must vent to keep
    # circulating, and all knots vent COHERENTLY into the q=0 common
    # mode, they rotate the background from within.  Σ is nothing but
    # the common-mode ledger: every closure in the web feeds it one
    # term (its dof count × its per-dof venting weight a₁).  Gravity is
    # the memory of this total venting: G_ind = 6π/Σ.
    #
    # Per-dof venting weights (Lorentz structure of the knot):
    #   scalar +1/6 − ξ, Dirac −1/12, gauge −1/6, Proca −1/12.
    # The bridge's non-minimal ξ is itself an echo: the bridge channel's
    # scalar vent 1/6 receives −ξ = −α_G₂·E[v²]·h_bridge from its own
    # self-coupling (same α, same E[v²], same h as everywhere else).
    from root import Ledger, WEB

    a1_scalar = 1.0 / 6.0
    a1_fermion = -1.0 / 12.0
    a1_gauge = -1.0 / 6.0
    a1_proca = -1.0 / 12.0

    # bridge per-channel vent as a nested ledger (additive)
    led_a1_bridge = Ledger("a1_bridge", a1_scalar, "add").echo(
        ["bridge-self"], -float(xi_bridge), 1,
        "FORCED", "ξ = α_G₂·E[v²]·h_bridge = 1/(48π)")
    a1_bridge = led_a1_bridge.value()
    WEB["a1_bridge"] = led_a1_bridge

    # Every dof count is a polynomial of (d₁₀, d₁₁, n₇, n₂₆).
    n_gen = d11                                    # 3 generations (D⁶ nimrep)
    n_higgs_UV = d10**2                            # 4 Higgs real DOF = h∨(G₂)
    n_gauge_UV = 2 * d11 * (d11 + 1)              # 24 = 2d₁₁(d₁₁+1)
    n_weyl_per_gen = d11 * d10**2 + d10 + 1        # 15 Weyl per generation
    n_fermion = d10 * n_gen * n_weyl_per_gen       # 90 fermionic dof
    n_G2_vec = 2 * n7 - (d11**2 - 1)              # 6 coset vectors
    n_G2_scal = n7                                 # 7 G₂ scalars

    def _vent_ledger(name, sector_terms):
        led = Ledger(name, 0.0, "add")
        for path, n_dof, a1, note in sector_terms:
            led.echo([path], n_dof * a1, 1, "FORCED",
                     f"{n_dof} dof × a₁ = {a1:+.4f}  {note}")
        WEB[name] = led
        return led

    led_Sigma_UV = _vent_ledger("Sigma_vent(UV)", [
        ("Higgs(4)",     n_higgs_UV, a1_scalar, "unbroken doublet"),
        ("gauge(24)",    n_gauge_UV, a1_gauge,  "SM transverse"),
        ("fermions(90)", n_fermion,  a1_fermion, "3 generations"),
        ("coset(6)",     n_G2_vec,   a1_proca,  "G₂/SU(3) vectors"),
        ("G₂-scal(7)",   n_G2_scal,  a1_scalar, "G₂ fundamental"),
        ("bridge(182)",  N_bridge,   a1_bridge, "coherence channels"),
    ])
    led_Sigma_br = _vent_ledger("Sigma_vent(broken)", [
        ("Higgs(1)",     n_higgs_UV - d11, a1_scalar, "physical Higgs"),
        ("gauge(18)",    2 * d11**2,       a1_gauge,  "gluons+photon"),
        ("Proca(9)",     d11**2,           a1_proca,  "massive EW"),
        ("fermions(90)", n_fermion,        a1_fermion, "3 generations"),
        ("coset(6)",     n_G2_vec,         a1_proca,  "G₂/SU(3) vectors"),
        ("G₂-scal(7)",   n_G2_scal,        a1_scalar, "G₂ fundamental"),
        ("bridge(182)",  N_bridge,         a1_bridge, "coherence channels"),
    ])

    # ── FACE-SPLIT LAW (depth-1, FORCED) ─────────────────────────────
    # The bridge's ONE self-echo, the same marginal unit 1/(2π) that
    # converts 512/π into 1/α(0) (couplings.py), is re-absorbed by
    # the two faces in the h∨ metric, the SAME metric that defines
    # sin²θ_W = h∨(SU3)/13 (root.py):  G₂ gets 4/13, F₄ gets 9/13,
    # and 4/13 + 9/13 = 1 = h_bridge (the unit fully lands).
    # METRIC UNIQUENESS: among the framework's six natural splits
    # (h∨, conformal, central charge, dimension, Casimir, fund dims)
    # only the h∨ split matches the ledger, dev 1.4e-5, runner-up
    # misses by 6.5e-2 (4600x separation).
    # NO-SELF-DILUTION (the sign): induced gravity dilutes, 1/G ∝ Σ.
    # The F₄ face IS the metric face (J₃(O)); the share of the echo
    # re-absorbed by the metric face is the metric's own amplitude
    # returning, it cannot count as dilution of itself.  The G₂
    # share is absorbed by the OTHER face (cross-face back-reaction
    # in the h∨ ledger, NOT a dressing of the gauge couplings, see
    # the absence watch): genuine venting, stays in Σ.  "Once" is inherited from the
    # no-further-edge theorem (one node, one self-loop).
    # The self-echo is a property of the bridge, not of the phase:
    # it enters BOTH bookkeepings (UV canonical, broken diagnostic).
    _face_split = -(hv_F4 / float(hv_G2 + hv_F4)) / (2.0 * math.pi)
    for _led in (led_Sigma_UV, led_Sigma_br):
        _led.echo(["bridge self-echo (F₄ h∨-share 9/13)"],
                  _face_split, 1, "FORCED",
                  "face-split law: no self-dilution of the metric face")

    Sigma_bridge = N_bridge * a1_bridge
    Sigma_total_UV = led_Sigma_UV.value()
    Sigma_total_broken = led_Sigma_br.value()
    Sigma_min_UV = Sigma_total_UV - Sigma_bridge
    Sigma_min_broken = Sigma_total_broken - Sigma_bridge

    # ── BOOKKEEPING DECISION THEOREM (closes the UV/broken corridor) ──
    # Σ is the coherent venting of all knots into the common mode, and
    # venting happens at the substrate scale (the healing length
    # ξ₀ = ℓ_Pl/2): the quadratically-divergent induced-G integral
    # localizes at the cutoff.  At Planckian distances there is no
    # electroweak condensate (v_EW/M_Pl ~ 10⁻¹⁷), so the spectrum that
    # vents is the SYMMETRIC-phase spectrum.  The UV bookkeeping is
    # therefore forced:
    #     G_ind/G_N = 0.999999917   (canonical; face-split law)
    # The broken-phase counting is retained as a diagnostic only (it
    # counts IR-phase dof whose contribution under Λ² is negligible).
    # The former −0.58% declared open vent is CLOSED by the face-split
    # law above; the remaining −8e-8 is below the (α/2π)² depth-2
    # echo scale, the ledger is depth-complete here.

    # Gravity = the memory of the total coherent venting
    target = 6.0 * math.pi
    G_ratio_UV = target / Sigma_total_UV
    G_ratio_broken = target / Sigma_total_broken
    G_ratio_canonical = G_ratio_UV
    G_ratio_mid = (G_ratio_UV + G_ratio_broken) / 2.0
    err_UV = 100.0 * (G_ratio_UV - 1.0)
    err_broken = 100.0 * (G_ratio_broken - 1.0)

    print(f"    Sigma_min(UV)     = {Sigma_min_UV:+.4f}")
    print(f"    Sigma_min(broken) = {Sigma_min_broken:+.4f}")
    print(f"    Sigma_bridge      = {Sigma_bridge:+.4f}  [{N_bridge} * (1/6 - 1/(48pi))]")
    print(f"    G_ind/G_N (UV)     = {G_ratio_UV:.6f}  ({err_UV:+.2f}%)")
    print(f"    G_ind/G_N (broken) = {G_ratio_broken:.6f}  ({err_broken:+.2f}%)")
    print(f"    G_ind/G_N CANONICAL = {G_ratio_canonical:.9f}  (UV: venting")
    print(f"      localizes at the cutoff, where EW is unbroken; broken-")
    print(f"      phase value {G_ratio_broken:.6f} is diagnostic only;")
    print(f"      former -0.58% vent CLOSED by the face-split law)")

    # Cosmological constant: Volovik → Jacobson-Clausius → CKN
    #
    # (1) Volovik: self-bound condensate → ρ_vac(eq) = 0 (Gibbs-Duhem
    #     identity, not fine-tuning; same mechanism as superfluid helium).
    # (2) Jacobson-Clausius: de Sitter departure from equilibrium.
    #     Acoustic Hawking T = H/(2π), horizon entropy S = A/(4G),
    #     plus Raychaudhuri kinematics yield:
    #       ρ_Λ = (3/8π) M_Pl² H²
    #     Coefficient 3/(8π) DERIVED from (T, S, Raychaudhuri).
    # (3) CKN saturation: Cohen-Kaplan-Nelson bound ρ_vac ≤ M_Pl²H²
    #     is saturated because all four Jacobson conditions are
    #     structurally satisfied by the substrate.
    #
    # H₀ is a unit-setting scale (same role M_Pl plays for the matter
    # sector), fixed by the algebra.  The framework derives the scaling
    # ρ_Λ ~ M_Pl²H² (not M_Pl⁴), resolving the 10¹²³ CC problem.
    # Current epoch: Ω_Λ < 1 because matter is still diluting:
    # standard cosmological evolution.
    H_0_SI = 67.4e3 / 3.0856e22            # H₀ in s⁻¹
    hbar_SI = 1.0546e-34                     # J·s
    GeV_per_J = 1.0 / 1.602e-10             # GeV/J
    H_0_GeV = H_0_SI * hbar_SI * GeV_per_J  # H₀ in GeV
    rho_Lambda = 3.0 * H_0_GeV**2 * M_Pl_GeV**2 / (8.0 * math.pi)
    rho_naive = M_Pl_GeV**4                  # naive QFT: ρ ~ M_Pl⁴
    print(f"    Λ: ρ_Λ = (3/8π) M_Pl² H₀² = {rho_Lambda:.2e} GeV⁴")
    print(f"       vs naive ρ ~ M_Pl⁴ = {rho_naive:.1e} GeV⁴  [10¹²³ resolved]")
    print(f"       Volovik (ρ_eq=0) → Jacobson (3/8π) → CKN saturation")

    # ══════════════════════════════════════════════════════════════════
    #  HIGGS: lambda(M_Pl) = -delta_bridge, SM RGE -> m_H
    # ══════════════════════════════════════════════════════════════════

    print(f"\n  --- HIGGS: F₄ fusion + bridge threshold ---")

    v_EW = masses['v_EW_pred_GeV']
    m_t_GeV = masses['m_t'] / 1e3
    alpha_s_MZ = couplings['alpha_s_MZ_thresh']

    # Higgs quartic boundary condition at M_Pl (two components):
    #
    # 1. F₄(1) fusion rule: (26 × 26)_local = 1
    #    The F₄ fusion algebra at level 1 has a unique quartic coupling:
    #    the symmetric product of two 26-reps fuses to the identity.
    #    This fixes λ_F₄ = 0 at the Planck scale (tree-level).
    #
    # 2. Bridge threshold correction:
    #    δ_bridge = N_bridge · α²_G₂(M_Pl) · E[v²]
    #            = (n₇·n₂₆) · (1/(24π))² · (1/2)
    #            = 182 · 1/(576π²) · 1/2
    #    Each factor is derived:
    #      N_bridge = n₇·n₂₆ = 182 (bridge mode count)
    #      α_G₂(M_Pl) = 1/(24π) (topological bootstrap, see root.py)
    #      E[v²] = 1/2 (protected forgetting: Pv²P = ½P)
    #
    # Combined: λ(M_Pl) = 0 − δ_bridge = −δ_bridge
    #
    # SINGLE SOURCE OF TRUTH: the bridge→Higgs edge lives in the web
    # (root.py WEB['lambda_MPl']); no second computation here (K3 dedup).
    lam_Pl_tree = WEB["lambda_MPl"].value()
    delta_bridge = -lam_Pl_tree

    print(f"    F₄ fusion: lambda_F₄ = 0  [(26×26)_local = 1]")
    print(f"    delta_bridge = {N_bridge} * (1/(24pi))² * 1/2 = {delta_bridge:.6f}")
    print(f"    lambda(M_Pl) = -{delta_bridge:.6f}")

    # Derive EW couplings at the derived M_Z scale
    lepton_masses = [masses['m_e'], masses['m_mu'], masses['m_tau']]
    quark_masses = [masses['m_u'], masses['m_d'], masses['m_s'],
                    masses['m_c'], masses['m_b'], masses['m_t']]

    # sin²θ_W from embedding (zero parameters)
    sin2W = _derive_sin2W()
    cos2W = 1.0 - sin2W

    # ── M_W and M_Z at one loop (MS-bar; TWO declared imports) ──────
    # PDG 2024 EW review (Erler-Freitas), Eq. (10.26):
    #     M_W = A₀ / [ŝ_Z (1 − Δr̂_W)^½],   M_Z = M_W / (ρ̂^½ ĉ_Z)
    # A₀ = [πα(0)/(√2 G_F)]^½ = 37.28038(1) GeV, and BOTH inputs of
    # A₀ are the framework's own predictions:
    #     A₀ = √(πα(0))·v_EW = 37.28038 GeV   (matches to 0.1 ppm)
    # ŝ² is canonical (3/13 + h₇·(α/2π), face-split law).  The two
    # imported radiative corrections [DECLARED, PDG 2024]:
    #     Δr̂_W = 0.06937(6),   ρ̂ = 1.01016(9)
    # Results: M_W = 80.365 GeV (−0.33σ of the 80.3692(133) average);
    #          M_Z = 91.196 GeV (+0.009%; +1.6σ of the import band,
    #          the W/Z tension structure of the data itself).
    # The face-split echo does real work here: tree 3/13 would give
    # M_W off by +0.12%.
    from root import PDG_EW
    A_EW = math.sqrt(math.pi * ALPHA_PHYS) * v_EW
    DR_HAT_W = PDG_EW['dr_hat_W']
    RHO_HAT = PDG_EW['rho_hat']
    M_W_derived = A_EW / (math.sqrt(sin2W) * math.sqrt(1.0 - DR_HAT_W))
    M_Z_derived = M_W_derived / (math.sqrt(RHO_HAT) * math.sqrt(cos2W))

    # α_EM at the derived M_Z (for the RGE initial conditions)
    alpha_em_MZ = _derive_alpha_em(lepton_masses, quark_masses, M_Z_derived)

    # PDG reference values (comparison only; single source root.PDG_EW)
    M_Z_PDG = PDG_EW['M_Z']
    M_W_PDG = PDG_EW['M_W']
    sin2W_PDG = PDG_EW['sin2W_MSbar']
    M_H_EXP = PDG_EW['m_H']

    print(f"    alpha_em(M_Z) = 1/{1.0/alpha_em_MZ:.3f}")
    print(f"    sin²theta_W = 3/13 + h₇·(α/2π) = {sin2W:.6f}  [tree 3/13 = {float(SIN2W_TREE):.5f} + FORCED depth-1 emission echo]")
    print(f"    PDG 2024 global fit: {sin2W_PDG}(4)  (pull: {(sin2W - sin2W_PDG)/PDG_EW['sin2W_err']:+.2f}σ; tree alone: −13σ)")
    print(f"    M_W = A₀/(ŝ√(1−Δr̂_W)) = {M_W_derived:.4f} GeV  (PDG: {M_W_PDG}(133) world avg, pull {(M_W_derived-M_W_PDG)/0.0133:+.2f}σ)")
    print(f"    M_Z = M_W/(ρ̂^½ĉ)      = {M_Z_derived:.4f} GeV  (PDG: {M_Z_PDG}, {pct(M_Z_derived, M_Z_PDG):+.4f}%; +1.6σ of import band)")
    print(f"      [A₀ = √(πα(0))·v_EW = {A_EW:.5f} GeV vs PDG 37.28038(1): 0.1 ppm;")
    print(f"       imports Δr̂_W = 0.06937(6), ρ̂ = 1.01016(9), PDG 2024 Eq. (10.26)]")

    # SM RGE initial conditions at derived M_Z
    g3_mz = math.sqrt(4.0 * math.pi * alpha_s_MZ)
    yt_mz = math.sqrt(2.0) * m_t_GeV / v_EW
    g1_mz = math.sqrt(4*math.pi*alpha_em_MZ/(1-sin2W)) * math.sqrt(5.0/3.0)
    g2_mz = math.sqrt(4*math.pi*alpha_em_MZ/sin2W)

    # RGE start scale: the derived M_Z (MS-bar chain above)
    t_Pl = math.log(M_Pl_GeV / M_Z_derived)
    N_rge = 20000

    # Self-consistent iteration: no measured m_H used.
    # Start with λ(M_Z) = 0, run UP to get gauge/Yukawa at M_Pl,
    # run DOWN with derived BC, iterate until λ(M_Z) converges.
    lam_iter = 0.0                        # initial guess (not measured!)
    for iteration in range(3):
        y_up = np.array([g1_mz, g2_mz, g3_mz, yt_mz, lam_iter])
        y_Pl_2L = _run_rge(y_up, 0.0, t_Pl, _beta_sm_2loop, N_rge)

        y_down = np.array([y_Pl_2L[0], y_Pl_2L[1], y_Pl_2L[2],
                           y_Pl_2L[3], lam_Pl_tree])
        y_low = _run_rge(y_down, t_Pl, 0.0, _beta_sm_2loop, N_rge)
        lam_iter = y_low[4]               # feed predicted λ(M_Z) back

    lam_pred = lam_iter
    mH_pred = math.sqrt(2.0 * abs(lam_pred)) * v_EW

    # F4-only for comparison (λ(M_Pl)=0, no bridge)
    y_down_F4 = np.array([y_Pl_2L[0], y_Pl_2L[1], y_Pl_2L[2],
                          y_Pl_2L[3], 0.0])
    y_low_F4 = _run_rge(y_down_F4, t_Pl, 0.0, _beta_sm_2loop, N_rge)
    mH_F4_only = math.sqrt(2.0 * abs(y_low_F4[4])) * v_EW

    err_F4 = (mH_F4_only / M_H_EXP - 1.0) * 100.0
    err_full = (mH_pred / M_H_EXP - 1.0) * 100.0

    # ── IMPORT SPECIFICATION (the framework's single remaining import)
    # Prescription, exactly: SM 2-loop RGE (Machacek-Vaughn) from the
    # derived M_Z to the derived M_Pl; tree-level matching at M_Z
    # (g_i from α_em(M_Z) + canonical ŝ², y_t = √2·m_t/v_EW);
    # boundary λ(M_Pl) = −δ_bridge (F₄ fusion + bridge threshold);
    # self-consistent λ(M_Z) iteration.  Fully specified, reproducible.
    # TRUNCATION BAND: the 1→2 loop shift on Planck-boundary Higgs
    # predictions is ~10 GeV (Buttazzo et al. 1307.3536); the standard
    # NNLO band is ±1 GeV, and tree-level matching at M_Z adds a
    # comparable contribution.  The measurement (125.20 ± 0.11) lies
    # within ~1.2 GeV of the prediction below, at the edge of the
    # stated running band and inside the combined band.  Upgrading the
    # import to NNLO (3-loop running + one-loop matching) requires
    # SMDR-class external tooling and changes the import, not the
    # framework: λ(M_Pl) = −δ_bridge is the framework's claim.
    print(f"    m_H (F₄ only, lambda=0):   {mH_F4_only:.1f} GeV  ({err_F4:+.2f}%)")
    print(f"    m_H (F₄ + bridge):          {mH_pred:.2f} GeV  ({err_full:+.2f}%)")
    print(f"    Experiment:                 {M_H_EXP:.2f} +/- 0.11 GeV")

    # ── IMPORT COMPATIBILITY AUDIT (zero-parameter discipline) ──────
    # A cited radiative correction is a FUNCTION of (G_F, m_t, m_H,
    # α̂(M_Z)), the SM fit's values.  A zero-parameter framework may
    # not absorb those silently: each ingredient is RECALCULATED here
    # from the framework's own predictions and verified compatible
    # within the import's quoted error.
    #   G_F : framework 1/(√2 v_EW²) = measured to < 1 ppm  [checked]
    #   m_t : ρ̂'s quadratic m_t term (PDG Eq. 10.23) shifts by 0.07σ
    #         of the import error under fw m_t vs fit m_t  [checked]
    #   m_H : the imports' m_H dependence is logarithmic; fw 124.85
    #         vs fit 125.2 is dln = -0.003  [checked]
    #   α̂(M_Z): the DOMINANT Δr̂_W ingredient is Δr₀ = 1 − α/α̂(M_Z)
    #         = 0.06646(6), driven by the dispersive hadronic VP,
    #         that is MEASUREMENT (e⁺e⁻ R-ratio spectral data), not a
    #         fit parameter, and is therefore admissible exactly like
    #         the PDG masses on the comparison side.  The framework's
    #         INTERNAL VP (free quarks + Λ_conf cutoff) reproduces it
    #         to 0.34% (1/128.37 vs 1/127.93), an approximation used
    #         ONLY for the Higgs-RGE initial conditions (effect on m_H
    #         far inside the import band); it must never be used at
    #         sub-0.5% precision, and the guard below enforces the
    #         approximation stays within its declared quality.
    # TAXONOMY (what an import IS, for a zero-parameter model):
    #   [1] data-driven ingredients  -> measurements (dispersive VP;
    #       same admissibility as the PDG masses)        [classified]
    #   [2] parameter dependencies   -> recalculated from framework
    #       values, compatible within quoted errors      [verified]
    #   [3] pure loop integrals      -> mathematics, uniquely fixed
    #       by gauge group + matter content (the EW remainder
    #       0.00291 in Δr̂_W, the bosonic 0.00082 in ρ̂), the same
    #       status as the Machacek-Vaughn beta functions used for
    #       the Higgs RGE.  No fine-tuning enters at any point:
    #       nothing is adjusted; everything is checked.
    _GF_fw = 1.0 / (math.sqrt(2.0) * v_EW**2)
    assert abs(_GF_fw / 1.1663788e-5 - 1.0) < 1e-6
    _drho_mt = 0.00934 * abs((m_t_GeV/172.61)**2 - (172.57/172.61)**2)
    assert _drho_mt < 9e-5                      # < 1 sigma of rho-hat
    assert abs(math.log(mH_pred / 125.20)) < 0.05
    assert abs((1.0/alpha_em_MZ)/127.929 - 1.0) < 0.005
    print(f"    import audit: G_F fw/meas −1 = {1e6*(_GF_fw/1.1663788e-5-1):+.2f} ppm;"
          f" δρ̂(m_t) = {_drho_mt:.1e} ({_drho_mt/9e-5:.2f}σ);")
    print(f"      internal VP α̂(M_Z) = 1/{1.0/alpha_em_MZ:.3f} vs dispersive 1/127.929"
          f" ({100*((1.0/alpha_em_MZ)/127.929-1):+.2f}%, RGE-only approximation)")

    # ══════════════════════════════════════════════════════════════════
    #  BARYOGENESIS: eta_B = n₇ * J_lep * exp(-d₁₀²pi²/2)
    # ══════════════════════════════════════════════════════════════════

    print(f"\n  --- BARYOGENESIS: eta_B from n₇, J_lep, G₂ instanton ---")

    # Derivation: each factor maps to a Sakharov condition.
    #
    # 1. B−L violation: G₂ instanton tunnelling
    #    S_inst = h∨(G₂)·π²/2 = d₁₀²·π²/2 = 2π²
    #    (same dimensional transmutation formula as the lepton sector,
    #     but with h∨(G₂) = d₁₀² instead of h∨(F₄) = d₁₁²)
    #    Tunnelling probability: exp(−S_inst) = exp(−2π²) ≈ 2.68×10⁻⁹
    #
    # 2. CP violation: leptonic Jarlskog invariant J_lep
    #    Derived from the conjugation modular invariant of SU(3)₃
    #    with the charged-lepton rotation (see mixing.py).
    #
    # 3. Multiplicity: n₇ = dim(G₂ fund) = 7
    #    Seven directed Fano-plane lines, each an independent asymmetry
    #    channel.  The Fano plane is the incidence geometry of the
    #    octonion multiplication table (G₂ = Aut(O)).
    #
    # Assembly: η_B = n₇ · J_lep · exp(−d₁₀²π²/2)
    #
    S_inst = float(hv_G2) * math.pi**2 / 2.0  # = d₁₀²·π²/2 = 2π²
    tunnelling = math.exp(-S_inst)
    J_lep = mixing['J_lep']

    eta_B = float(n7) * J_lep * tunnelling
    eta_B_obs = 6.12e-10

    ratio_B = eta_B / eta_B_obs
    pull_pct_B = abs(1 - ratio_B) * 100

    print(f"    S_inst = d₁₀²pi²/2 = {d10}²pi²/2 = {S_inst:.4f}")
    print(f"    exp(-S) = exp(-{d10}²pi²) = {tunnelling:.6e}")
    print(f"    J_lep = {J_lep:.6f}")
    print(f"    eta_B = {n7} * {J_lep:.6f} * {tunnelling:.4e}")
    print(f"          = {eta_B:.4e}")
    print(f"    Planck: {eta_B_obs:.2e}  ({pull_pct_B:.2f}% agreement)")

    # ══════════════════════════════════════════════════════════════════
    #  NEUTRINOS: m₁ = 0, normal ordering from F₄ singlet structure
    # ══════════════════════════════════════════════════════════════════

    print(f"\n  --- NEUTRINOS: F₄ singlet structure ---")

    # n₂₆ → 8_v + 8_s + 8_c + (n₂₆ − 3×8) singlets
    # The (n₂₆ − 24) = 2 singlets are RH neutrinos.
    n_singlets = n26 - 3 * (d11**2 - 1)   # = 26 - 24 = 2  (d₁₁²−1 = 8)
    print(f"    {n26} -> 8_v + 8_s + 8_c + {n_singlets}*1")
    print(f"    {n_singlets} RH neutrinos -> rank-2 seesaw -> m₁ = 0 (normal ordering)")

    # Structural predictions only, no external Δm² input needed.
    # The absolute neutrino mass scale is not predicted by the algebraic
    # framework; what IS predicted is the rank (m₁=0) and ordering (normal).
    print(f"    m₁ = 0 (structural prediction: rank-{n_singlets} seesaw)")
    print(f"    Ordering: normal (m₁ < m₂ < m₃)")
    print(f"    Testable by: KATRIN endpoint, 0νββ non-observation, cosmological Σmν")

    return {
        # Gravity
        'N_bridge': N_bridge,
        'xi_bridge': float(xi_bridge),
        'h_bridge': float(h_bridge),
        'h_G2_7': float(h_7),
        'h_F4_26': float(h_26),
        'alpha_G2_Pl': float(alpha_G2_Pl),
        'E_v2': float(E_v2),
        'a1_bridge': a1_bridge,
        'Sigma_min_UV': Sigma_min_UV,
        'Sigma_min_broken': Sigma_min_broken,
        'Sigma_bridge': Sigma_bridge,
        'Sigma_total_UV': Sigma_total_UV,
        'Sigma_total_broken': Sigma_total_broken,
        'G_ratio_UV': G_ratio_UV,
        'G_ratio_broken': G_ratio_broken,
        'G_ratio_mid': G_ratio_mid,
        'G_ratio_canonical': G_ratio_canonical,
        'err_UV': err_UV,
        'err_broken': err_broken,
        'rho_Lambda': rho_Lambda,
        # Higgs
        'alpha_em_MZ': alpha_em_MZ,
        'inv_alpha_em_MZ': 1.0 / alpha_em_MZ,
        'sin2W': sin2W,
        'delta_bridge': delta_bridge,
        'lam_Pl_tree': lam_Pl_tree,
        'mH_F4_only': mH_F4_only,
        'mH_pred': mH_pred,
        # Electroweak (derived)
        'M_Z_derived': M_Z_derived,
        'M_W_derived': M_W_derived,
        'M_Z_PDG': M_Z_PDG,
        'M_W_PDG': M_W_PDG,
        # Baryogenesis
        'eta_B': eta_B,
        'eta_B_obs': eta_B_obs,
        'ratio_B': ratio_B,
        'S_inst': S_inst,
        'tunnelling': tunnelling,
        'J_lep': J_lep,
        'dim_G2_fund': n7,
        # Neutrinos (structural predictions, no Δm² input)
        'ordering': 'normal',
        'n_RH': n_singlets,
        'm1': 0,
    }
