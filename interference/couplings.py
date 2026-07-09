"""
Couplings: alpha_s(M_Z) + alpha(0) from (d₁₀, d₁₁, n₇, n₂₆, pi).

alpha_s (full chain, self-contained, quark paper Sec. alpha-s):

  1. UV value: α_G₂(M_Pl) = |Z₃|/(2π D²_tot) = 1/(24π)  (SU(3)₃ MTC,
     executable in root.py).

  2. EXACT WZW CANCELLATION down to μ*.  Layer-1 G₂ back-reaction with
     C_A = h∨(G₂) = 4 and three Dirac fermions in the 7 (T(7) = 1):
         b₀ = (11/3)·4 − (4/3)·3 = 32/3
     over the GAUGE lever arm 9π²/2 − 6 (the lepton action with the
     Casimir vent).  The cancellation therefore terminates at the
     gauge-matching scale
         μ* ≡ M_Pl·exp(−(9π²/2 − 6)) ≈ 253.534 GeV
            (exactly v_EW·e^{15/512 − 16(α/2π)²}, the EM complement
            of the gauge lever arm in the S_quark ledger),
     NOT at v_EW = 246.2 GeV: the vertex-echo terms are EM physics
     below μ*.  They shift v_EW but are not part of the G₂ gauge
     running, so they do not participate in this identity.
         1/α_G₂(μ*) = 24π − (32/3)/(2π)·(9π²/2 − 6) = 32/π  EXACTLY
     (the π² in the action cancels the 1/2π of the loop integral).
     So α_G₂(μ*) = π/32 = 0.09817, and the SM running below starts
     at μ*.  (The α(0) chain is separate: π/512
     is the algebraic coupling of the embedding layer, defined by the
     Singh division of the WZW identity, not a running coupling
     evaluated at 246.2 GeV.)
     STATUS: a NON-PERTURBATIVE WZW IDENTITY.  The layer-2
     coefficient b₁ = (34/3)C_A² − ((20/3)C_A + 4C₂(7))·T(7)·n_f
     = 232/3 > 0 would shift α_s(μ*) up ~14%, consistent with
     π/32 being an exact operator identity (Knizhnik-Zamolodchikov
     conformal weights are exact).
     WHY DIRAC COUNTING (n_f = 3, not 6 Weyl): the framework's
     fermions are BdG quasiparticles (masses.py), and a BdG spectrum
     is particle-hole doubled BY CONSTRUCTION, each generation's
     excitation in the 7 is one full Dirac fermion, not a chiral
     half.  The 7 of G₂ is a real representation, consistent with
     exactly this non-chiral embedding.

  3. DERIVED THRESHOLD at G₂ → SU(3) (Weinberg 1980 / Hall 1981 matching).
     (Standard QFT calls this "threshold correction". Here the matching
      scale M_V is derived. See below.)
     The coset G₂/SU(3) (dim 6) gives six massive vectors in 3 ⊕ 3̄:
         1/α_s(μ*) = 1/α_G₂(μ*) − λ₃/(12π),
         λ₃ = (C_{G₂} − C_{SU(3)}) − 21·T_V·ln(M_V/μ*),
     with C_{G₂} − C_{SU(3)} = 4 − 3 = 1 and T_V = T(3)+T(3̄) = 1.
     M_V is DERIVED: the G₂-breaking scalar lives in the
     7 → 3 ⊕ 3̄ ⊕ 1. The singlet takes the VEV v_break = v_EW/√2
     (canonical complex-scalar normalisation. G₂→SU(3) breaking
     identified with EWSB).  The Casimir sum rule
         Σ_a M_a² = g² v_break² C₂(7) = 2 g² v_break²
     over six degenerate bosons gives
         M_V = g_{G₂} v_EW/√6 = π v_EW/(2√12) ≈ 112 GeV,
     where g_{G₂} = √(4πα_G₂) = π/(2√2).
     WHY THE COSET MODES DON'T APPEAR IN THE SPECTRUM: the D⁽⁶⁾
     boundary structure projects out the non-constant coset label
     before localisation, the off-diagonal nimrep fraction (2/3)
     transmits only scalar-like spectral weight, and no coset
     quantum number survives as a conserved charge below matching.
     M_V enters ONLY as a matching parameter (exactly like heavy
     states in GUT threshold back-reactions), not as a particle.

  4. SM layer-2 QCD running μ* → M_Z with derived thresholds
     → α_s(M_Z) = 0.1184 (+0.33%, +0.43σ of PDG 0.1180(9)).
     Without the threshold, π/32 + SM running alone gives 0.1119
     (−5.1%): the derived 112 GeV matching is load-bearing.

alpha(0):
    Bridge self-interference, full depth-3 form.  1/alpha(0) is the
    real root of the cubic
        x³ = (2⁹/pi)[(1 − 1/(2pi)) x² − 1/(2pi²)]
           = 137.035999050
    (depth-1 self-echo 1/(2pi) + the FORCED depth-3 e↔q two-orientation
    cycle, solved self-consistently in root.py).  Comparison, stated
    completely: Berkeley Cs 137.035999046(27) → +0.13σ. LKB Rb
    137.035999206(11) → −14σ. The two experiments disagree with each
    other by 5.5σ, and the framework's registered commitment is to
    Cs, with the Rb outcome as a registered kill condition.  The
    depth-1 TRUNCATION alone,
        (2⁹/pi)(1 − 1/(2pi)) = 2⁸(2pi−1)/pi² = 137.036439,
    is an intermediate identity, not the prediction.

Two readings, one coupling.  alpha_em = pi/2⁹ and alpha(0) = pi²/(2⁸(2pi-1))
are not "bare alpha and dressed alpha."  They are the SAME coupling's
back-reaction through the bridge sector, read at
two emergence layers of the conformal embedding:

  * alpha_em = pi/2⁹              at the conformal-embedding (Planck) scale.
                                  Used inside the instanton vertex.
  * alpha(0) = pi²/(2⁸(2pi-1))    at the IR pole-mass scale.
                                  Used in the lepton vertex factor.

The relation between them is NOT a perturbative running approximation.
It is a layer-1 identity of the (7,26) bridge sector, fixed by three
independent rational identities in the four integers:

    (A) h_bridge = h(n₇,G₂) + h(n₂₆,F₄) = 2/5 + 3/5 = 1     [MARGINAL]
    (B) D²_local = 1                                         [LAGRANGIAN]
    (C) c_coset  = c(E₈) − c(G₂) − c(F₄) = 0                 [TOPOLOGICAL]

The worldsheet integral of a marginal (h=1) primary with unit coupling
gives the universal h/(2pi) = 1/(2pi). c_coset = 0 forbids higher-layer
back-reactions.  The reading is layer-1 exact.

Same bridge, two observables.  gravity.py reads the same (7,26) sector
as 182 scalar-like heat-kernel channels with non-minimal coupling
ξ_bridge = α_G₂·E[v²]·h_bridge = 1/(48π) and derives Newton's constant.
The conformal weight h_bridge = 1 sets BOTH the QED layer-1 fraction
1/(2π) and the gravitational coupling 1/(48π).  EM renormalisation and
induced gravity are two readings through one bridge.

TRANSPORT.  The RGE running in this module (alpha_s to M_Z, the EW
chain, thresholds) carries already-derived web values to their
measurement scales.  The three stages are the graph (the web in
root.py, masses.py, gravity.py), the projections (mixing.py), and
transport (this module). Transport never feeds back into the graph.

Reference:
    Kahsay, Kibrom Kidane (2026). Three-paper series.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import math

from root import (d10, d11, n7, n26, h10, h11, K,
                  C2_26, hv_G2, hv_F4,
                  h_7, h_26, h_bridge,
                  alpha_G2_WZW, alpha_G2_Pl, alpha_phys, inv_alpha_phys,
                  M_Pl_GeV, S_lepton, N_bridge, xi_bridge, E_v2, pct)


# ═══════════════════════════════════════════════════════════════════════
#  RK4 integrator for RGE running
# ═══════════════════════════════════════════════════════════════════════

from root import rk4_run as _shared_rk4


def _rk4_run(beta_func, a0, t0, t1, n_steps=10000):
    """Scalar RGE running via the shared RK4 (root.rk4_run)."""
    return _shared_rk4(beta_func, a0, t0, t1, n_steps)


def _b0(nf):
    """Layer-1 QCD beta coefficient, derived from SU(d₁₁) Casimirs."""
    C_A = d11
    return (11 * C_A - 2*nf) / (6*math.pi)


def _b1(nf):
    """Layer-2 QCD beta coefficient, derived from SU(d₁₁) Casimirs.

      b₁ = (34·C_A² − (10·C_A + 6·C_F)·nf) / (24π²)
         = (153 − 19·nf) / (12π²)     for SU(3)

    No external numerical input: all numbers trace to d₁₀, d₁₁.
    """
    C_A = d11
    C_F = (d11**2 - 1) / (2.0 * d11)
    return (34*C_A**2 - (10*C_A + 6*C_F)*nf) / (24*math.pi**2)


def _run_SM_2loop(a0, mu0, mu1, nf):
    """Layer-2 SM QCD running in d/d(ln mu) convention.

    nf = d₁₀ × d₁₁ = 6 active flavours (2 types × 3 gen).
    """
    b0, b1 = _b0(nf), _b1(nf)
    return _rk4_run(lambda a: -b0*a**2 - b1*a**3,
                    a0, math.log(mu0), math.log(mu1))


# ═══════════════════════════════════════════════════════════════════════
#  Bridge self-interference: alpha(0)
# ═══════════════════════════════════════════════════════════════════════

def _derive_bridge():
    """Derive the physical fine-structure constant from bridge sector.

    Self-contained derivation from (d₁₀, d₁₁, n₇, n₂₆):

    1. Bridge conformal weights (from Casimirs and dual Coxeter numbers):
         h(7, G₂) = C₂(7)/(1+h∨(G₂)) = d₁₀/(1+d₁₀²) = 2/5
         h(26, F₄) = C₂(26)/(1+h∨(F₄)) = d₁₀d₁₁/(1+d₁₁²) = 3/5
         h_bridge = 2/5 + 3/5 = 1  (marginal operator)

    2. Lagrangian algebra condensation:
         The (7,26) bi-fundamental forms a Lagrangian algebra in E₈(1).
         The local quantum dimension is d_A = 1 + φ² = 2 + φ (golden ratio),
         giving D²_local = d_A²/d_A² = 1.
         This means the bridge self-coupling g² = D²_local = 1 (unit strength).

    3. Coset central charge:
         c_coset = c(E₈) − c(G₂) − c(F₄) = 8 − 14/5 − 26/5 = 0
         Vanishing c_coset makes the bridge sector TOPOLOGICAL:
         no propagating degrees of freedom → layer-1 exact (no higher
         back-reaction layers in the self-interference integral).

    4. Self-interference integral:
         A marginal primary (h=1) with unit coupling (g²=1) on a genus-0
         worldsheet gives the universal back-reaction:
           δ(1/α) = g²·h/(2π) = 1/(2π)
         This is a standard 2D CFT result for the integrated two-point
         function of a dimension-1 primary on the sphere.

    5. Assembly (depth 1, the truncation):
         1/α_alg = 2⁹/π = 512/π  (emergence value, from Singh ratio)
         (512/π)(1 − 1/(2π)) = 2⁸(2π−1)/π² = 137.036439
       The canonical 1/α(0) adds the FORCED depth-3 cycle and is the
       cubic root 137.035999050 (steps continue below).
    """
    # Lagrangian algebra condensation → D²_local = 1
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    d_A = 1.0 + phi**2                    # = 2 + phi (Lagrangian algebra dim)
    D2_local = (1.0 + phi**2)**2 / d_A**2  # = 1 (unit self-coupling)

    g2 = D2_local                          # bridge coupling² = 1
    h = float(h_bridge)                    # = 1 (marginal)
    delta_inv = g2 * h / (2.0 * math.pi)  # = 1/(2π) (layer-1 exact, c_coset=0)
    inv_alpha_alg = 2**9 / math.pi         # = 512/π (emergence value)

    # Depth-3 echo: the mutual electron↔quark loop through the EM channel.
    # ORIENTATION RULE: a closed echo cycle contributes once per
    # orientation (precedent: the published B = 2ρ factor from the
    # Hermitian pair S, S†, two orientations of one Z₃ cycle).  The
    # e↔q loop is a 2-cycle → multiplicity exactly 2.  The coupling
    # inside the loop is the physical α itself (phase of phase of
    # phase): self-consistent fixed point, depth 3 = the last layer the
    # Cayley-Dickson ledger can track (Hurwitz).
    #     1/α = (512/π)(1 − 1/(2π) − 2(α/2π)²)
    # equivalently the real root of
    #     x³ = (512/π)[(1 − 1/(2π))x² − 1/(2π²)]
    # SINGLE SOURCE OF TRUTH: the kernel's solved web state (root.py).
    # The fixed point of x ← b + W(x) already resums this cycle. No
    # second solver is kept here (dedup, task K3).
    inv_alpha = inv_alpha_phys
    alpha = alpha_phys

    # Verify: depth-1 truncation is the published 2⁸(2π−1)/π² identity,
    # and the full value solves the cubic.
    inv_depth1 = 2**8 * (2.0 * math.pi - 1.0) / math.pi**2
    assert abs(inv_alpha_alg * (1.0 - delta_inv) - inv_depth1) < 1e-10
    assert abs(inv_alpha**3 - inv_alpha_alg * ((1.0 - delta_inv)*inv_alpha**2
               - 1.0/(2.0*math.pi**2))) < 1e-6

    return {
        'h_bridge': h,
        'h_G2_7': float(h_7),
        'h_F4_26': float(h_26),
        'g_bridge_sq': g2,
        'd_A': d_A,
        'D2_local': D2_local,
        'c_coset': 0.0,
        'delta_inv_alpha': delta_inv,
        'inv_alpha_alg': inv_alpha_alg,
        'inv_alpha_phys': inv_alpha,
        'alpha_phys': alpha,
    }


# ═══════════════════════════════════════════════════════════════════════
#  Public interface
# ═══════════════════════════════════════════════════════════════════════

def derive(R, masses):
    """
    Derive alpha_s(M_Z) and alpha(0).

    R      : dict from root.derive()
    masses : dict from masses.derive()
    """
    print("\n" + "=" * 78)
    print("  COUPLINGS: alpha_s + alpha(0) from (d₁₀, d₁₁, n₇, n₂₆)")
    print("=" * 78)

    v_EW = masses['v_EW_pred_GeV']
    m_t_GeV = masses['m_t'] / 1e3
    m_b_GeV = masses['m_b'] / 1e3
    m_c_GeV = masses['m_c'] / 1e3
    # ─── PDG COMPARISON VALUES ───
    # M_Z is the PDG conventional scale at which alpha_s is quoted.
    # We evaluate alpha_s at this scale to compare with PDG's measurement.
    # The framework derives M_Z = 91.196 GeV via the MS-bar chain (gravity.py).
    from root import PDG_EW
    M_Z_PDG = PDG_EW['M_Z']               # PDG scale for alpha_s comparison
    alpha_s_PDG = PDG_EW['alpha_s_MZ']    # PDG comparison value only

    # ── alpha_s: exact WZW result at the gauge-matching scale μ* ──
    #   b₀(G₂) = (11/3)C_A − (4/3)ΣT_R = (11/3)(4) − (4/3)(3) = 32/3
    #   C_A = h∨(G₂) = 4,  matter = d₁₁=3 Dirac fermions in 7,  T(7)=1
    # The cancellation terminates where the GAUGE lever arm ends:
    #   μ* = M_Pl·exp(−(9π²/2 − 6)) = v_EW·e^{15/512} ≈ 253.5 GeV.
    # The 15/512 vertex echo is EM physics below μ* (it shifts v_EW,
    # not the G₂ gauge running), so π/32 holds at μ*, and the SM
    # running starts there.
    alpha_s_mu_star = math.pi / 2**5   # pi/32, exact at mu*
    mu_star = M_Pl_GeV * math.exp(-(S_lepton - float(C2_26)))
    # consistency: ln(mu*/v_EW) = 15/512 − 16(α/2π)², the EM physics
    # below mu* (depth-1 vertex echo minus the depth-3 e↔q vent), the
    # exact complement of the gauge lever arm in the S_quark ledger.
    em_below = 15.0/512.0 - 16.0*(alpha_phys/(2.0*math.pi))**2
    assert abs(math.log(mu_star / v_EW) - em_below) < 1e-9

    print(f"\n  alpha_G₂(mu*) = pi/2⁵ = pi/32 = {alpha_s_mu_star:.6f}")
    print(f"    mu* = M_Pl·e^-(9pi²/2-6) = {mu_star:.3f} GeV")
    print(f"        = v_EW·e^(15/512 - 16(a/2pi)²)  (the EM complement)")
    print(f"    (gauge lever arm endpoint. The 15/512 vertex echo is EM")
    print(f"     physics below mu* and does not feed the G₂ running)")
    print(f"    b₀(G₂) = (11/3)·4 − (4/3)·3 = 32/3  [C_A=h∨(G₂)=4, 3 Dirac in 7]")
    print(f"    Cancellation algebra:")
    print(f"    1/alpha(mu*)  = 24pi - (32/3)/(2pi) * (d₁₁²pi²/2 - C₂(26))")
    print(f"                  = 24pi - [24pi - 32/pi]")
    print(f"                  = 32/pi = {32/math.pi:.4f}")

    # ── No-threshold baseline ──
    n_f = d10 * d11                         # 6 quark flavours = 2 types × 3 gen
    a_no_th_mt = _run_SM_2loop(alpha_s_mu_star, mu_star, m_t_GeV, n_f)
    alpha_MZ_no_th = _run_SM_2loop(a_no_th_mt, m_t_GeV, M_Z_PDG, n_f - 1)
    err_no_th = 100 * (alpha_MZ_no_th - alpha_s_PDG) / alpha_s_PDG

    print(f"\n  SM layer-2 running (no threshold):")
    print(f"    alpha_s(M_Z) = {alpha_MZ_no_th:.4f}  ({err_no_th:+.1f}%)")

    # ── G₂ -> SU(3) threshold (matched at μ*) ──
    g_G2 = math.sqrt(4.0 * math.pi * alpha_s_mu_star)  # = pi/(2*sqrt(2))
    # M_V = g_G₂ × v_EW / √C₂(26),  C₂(26) = d₁₀d₁₁ = 6.
    # M_V is tied to the PHYSICAL EWSB scale v_EW (the coset vectors
    # take their mass from the breaking VEV), while the Weinberg-Hall
    # log is referenced to the matching scale μ*.
    M_V = g_G2 * v_EW / math.sqrt(float(C2_26))        # = pi*v_EW/(2*sqrt(12))

    # Dictionary forms (were bare floats): T_V = T(3)+T(3̄) = 2h₁₁ = 1.
    # C_diff = C_A(G₂) − C_A(SU3) = d₁₀² − d₁₁ = 1.
    T_V = 2.0 * float(h11)                 # Dynkin index of coset vectors
    C_diff = float(d10**2 - d11)           # matching constant
    assert T_V == 1.0 and C_diff == 1.0
    # 21 = d₁₁ × n₇ = 3 × 7: coset coefficient for G₂ → SU(d₁₁) threshold
    lambda_3 = C_diff - float(d11 * n7) * T_V * math.log(M_V / mu_star)

    inv_as_G2_mu_star = 1.0 / alpha_s_mu_star          # 32/pi
    inv_as_thresh = inv_as_G2_mu_star - lambda_3 / (12.0 * math.pi)
    a_thresh = 1.0 / inv_as_thresh

    a_mt_th = _run_SM_2loop(a_thresh, mu_star, m_t_GeV, n_f)
    alpha_MZ_thresh = _run_SM_2loop(a_mt_th, m_t_GeV, M_Z_PDG, n_f - 1)
    err_thresh = 100 * (alpha_MZ_thresh - alpha_s_PDG) / alpha_s_PDG

    print(f"\n  G₂ -> SU(3) threshold (matched at mu*):")
    print(f"    g_G₂ = pi/(2sqrt(2)) = {g_G2:.6f}")
    print(f"    M_V = g*v_EW/sqrt(6) = {M_V:.3f} GeV  (derived, 0 params)")
    print(f"    lambda_3 = 1 - 21·ln(M_V/mu*) = {lambda_3:.6f}")
    print(f"    alpha_s(mu*) after threshold = {a_thresh:.6f}")
    print(f"    alpha_s(M_Z) [+threshold] = {alpha_MZ_thresh:.4f}  ({err_thresh:+.2f}%)")
    print(f"    PDG: 0.1180 +/- 0.0009")
    print(f"    Delta = {abs(alpha_MZ_thresh-alpha_s_PDG):.4f}"
          f", {abs(alpha_MZ_thresh-alpha_s_PDG)/0.0009:.1f}sigma")

    # ── RK4 convergence self-certification ──
    # Re-run the full chain at 2× step count. Assert agreement.
    a_mt_2x  = _rk4_run(lambda a: -_b0(n_f)*a**2 - _b1(n_f)*a**3,
                         a_thresh, math.log(mu_star), math.log(m_t_GeV), 20000)
    a_MZ_2x  = _rk4_run(lambda a: -_b0(n_f-1)*a**2 - _b1(n_f-1)*a**3,
                         a_mt_2x, math.log(m_t_GeV), math.log(M_Z_PDG), 20000)
    assert abs(alpha_MZ_thresh - a_MZ_2x) < 1e-6, \
        f"RK4 not converged: {alpha_MZ_thresh} vs {a_MZ_2x}"

    # ── Bridge self-interference: alpha(0) ──
    bridge = _derive_bridge()

    print(f"\n  Bridge self-interference -> alpha(0):")
    print(f"    h_bridge = h({n7},G₂) + h({n26},F₄) = {float(h_7)} + {float(h_26)} = {float(h_bridge)}")
    print(f"    D²_local = {bridge['D2_local']:.1f}  (Lagrangian condensation)")
    print(f"    c_coset = {bridge['c_coset']:.0f}  (topological -> layer-1 exact)")
    print(f"    depth 1: 1/alpha = (2⁹/pi)(1 - 1/(2pi)) = 137.036439")
    print(f"    depth 3: e↔q loop, 2 orientations (S,S† rule), self-consistent:")
    print(f"    1/alpha(0) = (2⁹/pi)(1 - 1/(2pi) - 2(alpha/2pi)²)")
    print(f"               = {bridge['inv_alpha_phys']:.9f}")
    print(f"    Berkeley Cs (2018):  137.035999046(27)   "
          f"{(bridge['inv_alpha_phys']-137.035999046)/27e-9:+.2f}σ")
    print(f"    LKB Rb (2020):       137.035999206(11)   "
          f"{(bridge['inv_alpha_phys']-137.035999206)/11e-9:+.2f}σ")
    print(f"    CODATA blend:        137.035999177(21)   "
          f"{(bridge['inv_alpha_phys']-137.035999177)/21e-9:+.2f}σ")
    print(f"    The two recoil experiments disagree by >5σ (unresolved).")
    print(f"    the framework lands on the cesium side, falsifiable call.")

    return {
        'mu_star_GeV': mu_star,
        'alpha_s_mu_star': alpha_s_mu_star,
        'alpha_s_mu_star_threshold': a_thresh,
        'alpha_s_MZ_no_thresh': alpha_MZ_no_th,
        'alpha_s_MZ_thresh': alpha_MZ_thresh,
        'M_V_derived': M_V,
        'lambda_3': lambda_3,
        'err_no_thresh': err_no_th,
        'err_thresh': err_thresh,
        # Bridge data
        'bridge': bridge,
        'inv_alpha_phys': bridge['inv_alpha_phys'],
        'alpha_phys': bridge['alpha_phys'],
    }
