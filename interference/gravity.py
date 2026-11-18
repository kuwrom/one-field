"""
Gravity, baryogenesis and neutrinos (the Higgs transport is in higgs_ew.py)
from (d₁₀, d₁₁, n₇, n₂₆, pi)
+ the electron anchor (M_Pl is DERIVED in root.py from m_e).

Import inventory: higgs_ew.py (the consolidated table lives there).

═══════════════════════════════════════════════════════════════════════
WHAT GRAVITY IS
═══════════════════════════════════════════════════════════════════════

A matter standing wave is a circulating soliton.  Its circulation
biases the background shared by the substrate's three components.

The exact Z₃ Fourier split separates the q = 1, 2 relative modes from
the q = 0 common mode.  The geometry is the common mode at second
order: the gradient bilinear summed over the components,
    h_μν = Σ_a (∂_μ φ̄_a ∂_ν φ_a + c.c.).
This is the metric of a field whose action carries no background metric
(tests/probes/composite_metric.py).  The DeWitt–Seeley heat kernel
theorem supplies its induced Einstein–Hilbert term.  The bridge ledger
fixes that term's sign and form. Its size moves with the regulator, so
the strength of gravity is the anchor inversion (root.py).

In the empty field a particle is a self-trapped lump with the clock
inside it (tests/probes/lump_in_empty_field.py).  Gravity is the common
geometry induced by matter's circulation.

═══════════════════════════════════════════════════════════════════════
WHY GRAVITY AT ALL: the common mode memory
═══════════════════════════════════════════════════════════════════════

The substrate is the Z₃-symmetric coupled NLS system (masses.py).
In Madelung variables ψ_k = √ρ_k e^{iS_k/ħ}, perturbations decompose
under the Z₃ Fourier transform into:
    q = 0    the COMMON mode, shared by all three components.
    q = 1,2  the RELATIVE modes, they distinguish components.
The q = 1, 2 modes carry matter as standing waves. The q = 0 mode
carries the gravitational memory.  This resolves the Barceló-Liberati-Visser
"two roles of ρ" obstruction of analogue gravity: matter and
geometry are different Fourier projections of ONE substrate.

The metric is the gradient bilinear, a symmetric tensor, so the
conformally flat acoustic metric of (ρ, v) does not arise.  At first
order in the field no tensor exists (tests/probes/orbital_lock.py). The
tensor is the second order bilinear, the first order that PvP = 0 does
not remove.  A lump's own field falls off exponentially outside it, so
the long range channel is the tensor mode of h. A lump sources that
mode with its whole energy, binding energy included
(tests/probes/two_lumps.py).
The ledger applies the Sakharov heat kernel formula to everything that
propagates on h and fixes the sign of the Einstein-Hilbert
coefficient. The minimal matter spectrum alone gives the WRONG SIGN
(Σ_min = −10/3). The 182 bridge channels supply the positive
backreaction (+182/6).  This module computes the sum with the bridge
included.

Gravity (bridge sector):
    N_bridge = n₇ * n₂₆ = 182, minimally coupled (boundary.py)
    Heat kernel -> the sign of 1/G. G is the anchor inversion (root.py)

COMMON MODE CONSEQUENCES.
    The vacuum is the empty field and carries no profile, so every row
    propagates on the carrier's one cone (gate C below).  PvP = 0 is an
    algebraic identity of the Z₃ harmonic (protected_forgetting.py) and
    supplies the universal common mode coupling used by the gravity and
    dark sector branches.

Higgs + EW: higgs_ew.py (SM RGE transport, M_W, M_Z, import audit).

Baryogenesis assembly (the depth 4 G₂ face action):
    eta_B = −n₇ × J_lep × exp(−h∨(G₂)π²/2) = −7 · J_ℓ · e^{−2π²}
      • e^{−2π²}: root.py derives the G₂ face action from the F₄
        dimensional transmutation action: b_face=(32/3)(9/4)=24=
        h∨(E8)-C2(26), hence S_G2=2pi/(24 alpha_G2)=2pi². The same
        derived factor is used for the
        Majorana scale M_R = M_Pl e^{−2π²} ≈ 3.3×10¹⁰ GeV >> T_EW,
        so out of equilibrium decay (3rd condition) needs no extra
        ingredient.
      • J_lep: leptonic Jarlskog (mixing.py).  CP violation traces
        to the ORIENTED Fano plane: seven directed lines select a
        common CP orientation sign.
      • n₇ = 7: the asymmetry channels, one per directed Fano
        line, each one octonionic multiplication.
    The asymmetry is generated NON THERMALLY, at the
    instanton level. It is a topological count, not a process that
    needs a decay era.  The Majorana scale M_R = M_Pl e^{−2π²} >> T_EW means
    the B−L asymmetry is frozen in before any thermal washout era
    exists. The Boltzmann/washout machinery of thermal leptogenesis
    describes a different mechanism and does not apply here.

Neutrinos:
    26 -> 8_v + 8_s + 8_c + 1 + 1   (SO(8) triality decomposition
    of the F₄ fundamental: three 8's tied by triality + 2 singlets)
    In a minimal type-I seesaw with no other light neutrino mass terms,
    identifying the 2 singlets with RH neutrinos bounds the light neutrino
    mass matrix to rank at most 2, hence at least one massless state.
    This conditional rank bound does not select the ordering.

═══════════════════════════════════════════════════════════════════════
Why the bridge contributes a SCALAR sign (not a vector sign)
═══════════════════════════════════════════════════════════════════════

The bridge sector has 182 channels. Their scalar heat kernel reading
follows the framework's route: Lagrangian
algebra condensation absorbs the bridge into the vacuum module. The
channel has no adjoint gauge redundancy or ghost sector. It enters
the metric's stiffness as a scalar common mode row.  The amount carries
no names, and a sector with no names builds a metric from scalar
gradients alone (tests/probes/unnamed_metric.py). Given that routing,
the scalar heat kernel coefficient follows.
Reading the same channels as vectors would contradict this routing and
would give −1/3 per state in the ledger.

The heat kernel a₁ coefficient for each bridge channel is
    a₁/R = +1/6   (scalar, minimally coupled, ATTRACTIVE)
NOT the vector value −1/3 per state. Five linked checks establish this reading
within the framework:

  (i)   Lagrangian algebra condensation:
        A_{E₈} = 1 + (tau_G, tau_F) absorbs the bridge into the E₈
        vacuum module with D²_local = 1.  The condensed modes are
        coherence channels of the parent theory, not independent
        propagating gauge fields.

  (ii)  Zero coset central charge:
        c_coset = c(E₈) − c(G₂) − c(F₄) = 0
        In the unitary conformal embedding this removes an independent
        coset Virasoro sector. The absence of vector ghost subtraction
        comes from the separate scalar/nonadjoint routing in (i), (iii).

  (iii) Representation theory:
        The bridge transforms as the bifundamental (n₇, n₂₆), NOT as
        the adjoint of any gauge group.  No gauge redundancy means
        no Faddeev-Popov ghosts, so the ghost subtracted vector sign
        does not apply.

  (iv)  Heat kernel ghost arithmetic (per massless vector, summed over
        its two states): −2/3 = 4·(1/6) (four components) − 1 (Ricci
        endomorphism of the Hodge Laplacian) − 2·(1/6) (two ghosts).
        A channel with no gauge redundancy has neither the endomorphism
        nor the ghosts, so its coefficient is the scalar +1/6.

  (v)   Quantum dimension verification:
        d(tau_G) * d(tau_F) = phi², d(A_{E₈}) = 1 + phi² = 2 + phi.
        FP dimension confirms the bridge is fully absorbed:
        D²_local = D²(C)/d(A)² = (1+phi²)²/(2+phi)² = 1, the
        Kirillov-Ostrik dimension of the local modules, computed from
        the four simple objects (couplings.py._derive_bridge).

═══════════════════════════════════════════════════════════════════════
Same bridge, two observables (cross link to couplings.py)
═══════════════════════════════════════════════════════════════════════

The (7,26) bridge carries two readings through the same h_bridge = 1:
  * Gravity: 182 scalar heat kernel channels each contribute
    a₁/R = 1/6, minimally coupled (boundary.py). Newton's constant is
    the anchor inversion read through the common mode (root.py).
  * QED:    The same h_bridge = 1 enters the normalized current response
    h_bridge/(2π) = 1/(2π), giving (1 − 1/(2π)) that converts
    alpha_em = pi/2⁹ into alpha(0) = 1/137.035999050.  See couplings.py.
Same sector, same h_bridge, two observables.

Reference:
    Kahsay, Kibrom Kidane (2026). The Octavian Coherence Gate: The Four
    Irreducible Integers of the E₈(1) ⊃ G₂(1) × F₄(1) Conformal Embedding.
    Zenodo. https://doi.org/10.5281/zenodo.20493955
"""

import cmath
import math

from . import boundary as _boundary
from .root import (d10, d11, n7, n26,
                   N_bridge, h_bridge, h_7, h_26,
                   alpha_G2_Pl, xi_bridge, E_v2,
                   M_Pl_GeV, b0_G2,
                   BA_ratio, h10)
from .higgs_ew import derive_higgs_ew


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

    print(f"\n  --- GRAVITY: {N_bridge} bridge modes, the amount minimally coupled ---")

    # Protected forgetting on the canonical Z₃ harmonic
    #   v_k = B cos(θ + 2πk/3),  B = √d₁₀ = √2,  θ = h₁₀ = 2/9
    # (see protected_forgetting.py for the moment bookkeeping).
    #
    #   PvP    = Σv_k/3   = 0        forgetting, holds for every (B, θ)
    #   Pv²P   = Σv_k²/3  = B²/2 = 1 TOTAL second moment (both channels)
    #   E[v²]  = |v̂₁|²    = B²/4 = ½ PER MODE power, one Z₃ channel
    #
    # The unique (7,26) bridge carries this G2/Z3 response tensored
    # with identity on the F4 26 (wiring.py W6, Schur trace), so the
    # ledger uses the per mode power E[v²] = ½, not the total 1.
    # The two are easy to confuse: the notation "Pv²P = ½P" names
    # the total, not the per mode power.
    # The canonical harmonic comes from ONE source (root.BA_ratio =
    # √d₁₀, root.h10 = 2/9), so a drift in either parameter fails the
    # asserts below.
    v = [BA_ratio * math.cos(float(h10) + 2*math.pi*k/3)
         for k in range(3)]
    u = [1.0 / math.sqrt(3.0)] * 3
    u_dot_v = sum(ui * vi for ui, vi in zip(u, v))
    Pv2P_total = sum(ui**2 * vi**2 for ui, vi in zip(u, v))
    v_hat1 = sum(vk * cmath.exp(-2j*math.pi*k/3)
                 for k, vk in enumerate(v)) / 3.0
    E_v2_check = abs(v_hat1)**2
    assert abs(u_dot_v) < 1e-12
    assert abs(Pv2P_total - 1.0) < 1e-12
    assert abs(E_v2_check - 0.5) < 1e-12
    assert abs(Pv2P_total - 2.0 * E_v2_check) < 1e-12

    print(f"    Protected forgetting: PvP = 0, "
          f"Pv²P = 1 (total) = 2·E[v²], E[v²] = 1/2 (per Z₃ channel)")
    print(f"    h_bridge = h({n7},G₂) + h({n26},F₄) = {float(h_7)} + {float(h_26)} = {float(h_bridge)}")
    print(f"    the label's xi_bridge = alpha_G₂(M_Pl) * E[v²] * h_bridge (used only by the labelled bridge)")
    print(f"              = 1/(24pi) * 1/2 * 1 = 1/(48pi) = {float(xi_bridge):.10f}")

    # Heat kernel coefficients (from spin content of E₈ decomposition)
    #
    # The E₈ → G₂ × F₄ decomposition determines which fields exist and
    # their spins.  The a₁/R coefficient for each spin follows from the
    # representation of SO(3,1). It is a mathematical consequence of the
    # field's Lorentz structure, just as a Casimir follows from a rep's
    # weight diagram.  Per physical state, in the supertrace that
    # induces 1/G (fermion loops carry the relative minus sign):
    #   a₁(scalar)  = +1/6   (spin-0, minimal coupling ξ=0)
    #   a₁(Weyl)    = +1/12  (spin-1/2: −1/6 per Weyl pair of states,
    #                          sign flipped by the fermion loop)
    #   a₁(gauge)   = −1/3   (massless spin-1: −2/3 per two states)
    #   a₁(Proca)   = −1/6   (massive spin-1: −1/2 per three states,
    #                          exactly massless vector + eaten scalar)
    # Summed over a spectrum this is the standard (N₀ + N_½ − 4N₁)/6
    # (Visser gr-qc/0204062 Table 1; Vassilevich hep-th/0306138 §4;
    # overall sign: Larsen and Wilczek hep-th/9506066 Eqs. 3.8, 3.17, 3.20).
    #
    # The (7,26) label's nonminimal coupling ξ_bridge = 1/(48π) is
    # α_G₂(M_Pl)·E[v²]·h_bridge (see root.py). The amount, with no share
    # of b₀, does not carry it (boundary.py, THE AMOUNT'S COUPLING).
    # The bridge modes are scalar (not vector) because: c_coset = 0,
    # no gauge redundancy (bifundamental, not adjoint), and D²_local = 1
    # (Lagrangian algebra condensation).
    #
    # Credit: DeWitt (1965), Seeley (1967) for the heat kernel expansion.
    #
    # # THE UNIVERSAL VENTING LEDGER ################
    # A standing wave is a circulating knot. The part that does not
    # cancel locally vents COHERENTLY into the q = 0 common mode and
    # rotates the background from within. Σ is the common mode ledger:
    # every closure in the web adds one term to it, its dof count × its
    # per dof venting weight a₁. Gravity is the memory of this total
    # venting. Σ fixes the sign of gravity, and its strength is the
    # anchor inversion (root.py).
    #
    # Venting weights per state (Lorentz structure of the knot):
    #   scalar +1/6, Weyl +1/12, gauge −1/3, Proca −1/6.
    # The amount vents at the minimal scalar weight. A row with no share
    # of b₀ has no G₂ vertex, so no self-coupling echo ξ is entered on it
    # (boundary.py, THE AMOUNT'S COUPLING). The label's
    # ξ = α_G₂·E[v²]·h_bridge is used only by the labelled bridge.
    from .root import Ledger, WEB

    a1_scalar = 1.0 / 6.0
    a1_fermion = 1.0 / 12.0
    a1_gauge = -1.0 / 3.0

    # the bridge's venting weight per channel, as a nested additive ledger
    led_a1_bridge = Ledger("a1_bridge", a1_scalar, "add", web_depth=1)
    a1_bridge = led_a1_bridge.value()
    WEB["a1_bridge"] = led_a1_bridge

    # Every dof count is a polynomial of (d₁₀, d₁₁, n₇, n₂₆).
    n_gen = d11                                    # 3 generations (D⁶ nimrep)
    n_gauge_UV = 2 * d11 * (d11 + 1)              # 24 = 2d₁₁(d₁₁+1)
    n_weyl_per_gen = d11 * d10**2 + d10 + 1        # 15 Weyl per generation
    n_fermion = d10 * n_gen * n_weyl_per_gen       # 90 fermionic states
    n_G2_vec = d10 * (2 * n7 - (d11**2 - 1))      # 6 coset vectors × 2 states
    n_G2_scal = n7                                 # 7 G₂ scalars

    def _vent_ledger(name, sector_terms):
        led = Ledger(name, 0.0, "add", web_depth=1)
        for path, n_dof, a1, note in sector_terms:
            led.echo([path], n_dof * a1, 1, "FORCED",
                     f"{n_dof} dof × a₁ = {a1:+.4f}  {note}",
                     kind="vent")
        WEB[name] = led
        return led

    # The full ledger's rows are boundary.py's one list with the 7 G₂
    # scalars added as light rows. The running of α_G₂ (b₀ = 32/3,
    # root.py) does not count these scalars. The two lists differ by
    # those scalars alone, 7/6 at every contact. FACT 7's reference
    # uses these rows.
    led_Sigma_UV = _vent_ledger("Sigma_vent(UV)", [
        ("gauge(24)",    n_gauge_UV, a1_gauge,  "SM transverse"),
        ("fermions(90)", n_fermion,  a1_fermion, "3 generations"),
        ("coset(12)",    n_G2_vec,   a1_gauge,  "G₂/SU(3) vectors, massless"),
        ("G₂-scal(7)",   n_G2_scal,  a1_scalar, "G₂ fundamental"),
        ("bridge(182)",  N_bridge,   a1_bridge, "coherence channels"),
    ])
    # Broken phase reading. The ledger carries no elementary Higgs at
    # the regulator scale, because the electroweak scale is transmuted
    # (masses.py). A Proca field counts as a massless vector plus a
    # composite scalar (Visser gr-qc/0204062 eq. 46), so the two ledgers
    # are identical row by row (tests/probes/regulator_free_ratios.py R6).
    led_Sigma_br = _vent_ledger("Sigma_vent(broken)", [
        ("gauge(18)",    2 * d11**2,       a1_gauge,  "gluons+photon"),
        ("W, Z(6)",      2 * d11,          a1_gauge,  "massive EW, transverse at the regulator"),
        ("fermions(90)", n_fermion,        a1_fermion, "3 generations"),
        ("coset(12)",    n_G2_vec,         a1_gauge,  "G₂/SU(3) vectors, massless"),
        ("G₂-scal(7)",   n_G2_scal,        a1_scalar, "G₂ fundamental"),
        ("bridge(182)",  N_bridge,         a1_bridge, "coherence channels"),
    ])

    # Σ is only the sum of its rows. The F₄ h∨ share of the
    # bridge self-echo returns to the circulation of the metric face
    # and is entered on no ledger (root.py, the apportionment law).

    WEB.solve()          # gravity nodes join the solved web state

    Sigma_bridge = N_bridge * a1_bridge
    Sigma_total_UV = led_Sigma_UV.value()
    Sigma_total_broken = led_Sigma_br.value()
    Sigma_min_UV = Sigma_total_UV - Sigma_bridge

    # # The phase of the ledger #######################
    # The induced G integral localizes at the regulator, where there is
    # no electroweak condensate (v_EW/M_Pl ~ 10⁻¹⁷).  With the standard
    # per state weights the symmetric and broken ledgers are equal term
    # by term, asserted below.
    #
    # 6π/Σ is a diagnostic ratio under the cutoff t₀ = d₁₁ℓ_Pl²
    # (tests/probes/sakharov_normalisation.py FACT 7), not a reading of
    # G.  G is root.py's anchor inversion read through the common mode.
    target = 6.0 * math.pi
    G_ratio_UV = target / Sigma_total_UV
    G_ratio_broken = target / Sigma_total_broken
    assert abs(G_ratio_UV - G_ratio_broken) < 1e-12, \
        "Proca = massless vector + eaten scalar: the two ledgers must agree"
    G_ratio_canonical = G_ratio_UV

    # # the diagnostic ratio 6π/Σ as a graph node (a live ratio base) ####
    # Σ is the total of the venting ledger, taken from the state. The
    # ratio is FACT 7's diagnostic on the full ledger, not a reading of G.
    WEB["G_ratio"] = Ledger(
        "G_ratio", lambda s: 6.0 * math.pi / s["Sigma_vent(UV)"], "mul",
        web_depth=1)
    WEB.solve()          # gravity nodes join the solved web state
    assert abs(WEB.state["G_ratio"] / G_ratio_canonical - 1.0) < 1e-12

    Sigma_matter = Sigma_min_UV
    print(f"    Sigma_matter      = {Sigma_matter:+.4f}  [(N0 + N½ − 4 N1)/6, symmetric phase]")
    print(f"    Sigma_bridge      = {Sigma_bridge:+.4f}  [{N_bridge} * 1/6, the amount minimally coupled]")
    print(f"    Sigma_total       = {Sigma_total_UV:.4f}  [full ledger]")

    # # THE GRAVITY LANE: closure gates, not a coefficient ##########
    # Gravity is read the same way as the matter sector, as a closure.
    # The bias that persists on the common mode is the one that passes
    # the gates below. They are the analogue of the L and Q gates of the
    # matter lane. The strength of gravity is the anchor inversion,
    # m_e → M_Pl, read through the common mode at depth four:
    # G M_b² = 3 Q_l − 1 (root.py, THE COMMON MODE).
    #
    # THE LANE RESTS ON three statements.  The action carries no
    # background metric, only one fixed volume element whose conjugate
    # is Lambda.  The couplings have the self-bound sign g₀ < 0 < g₁ at
    # |g₁/g₀| = 1/√d₁₀ (tests/probes/self_bound_vacuum.py).  The
    # vacuum is the empty field, the engine's zero orbit.
    # WHAT FALSIFIES IT. A tensor speed off c (GW170817 bounds 1e-15), a
    # scalar or vector polarization of gravitational waves, or two
    # clocks of different make whose ratio records a bias
    # (tests/probes/gates_as_records.py).
    N0_st = n_G2_scal                         # real scalars
    Nh_st = n_fermion // 2                    # Weyl fermions
    N1_st = (n_gauge_UV + n_G2_vec) // 2      # massless vectors
    ledger_ratio = {}
    for contact in (1.0, 0.0, -1.0):          # vector contact coefficient
        sig_c = ((N0_st + Nh_st) / 6.0 + 2 * N1_st * (1.0 / 6 - contact / 2.0)
                 + Sigma_bridge)
        ledger_ratio[contact] = target / sig_c
    assert abs(ledger_ratio[1.0] - G_ratio_canonical) < 1e-12
    str_one = (N0_st + N_bridge + 2 * N1_st) - n_fermion  # bosonic − fermionic states
    lambda_G = -(math.pi / 2.0) * str_one / Sigma_total_UV**2
    # # ONE LIST: the rows the running uses (boundary.py) ##########
    # The gates use the same list as b₀ (= 32/3, root.py). The bridge
    # vents as amount and does not run. The regulator is solved from
    # the anchored G, t₀ = Σ/(2π) M_b⁻². The bridge with its
    # searchlight label is kept for comparison (boundary.py, rows()
    # with bridge="light" or bridge="boundary").
    one = _boundary.rows(d10, d11, n7, n26)
    assert _boundary.b0(one) == b0_G2, "the running and the stiffness read one list"
    _booked = _boundary.rows(d10, d11, n7, n26, scalars="light")
    assert abs(_boundary.sigma(_booked, 0.0, 1.0, float(xi_bridge)) - Sigma_total_UV) < 1e-12
    one_list = {}
    for contact in _boundary.CONTACTS:
        roots = _boundary.anchored(one, contact, float(xi_bridge))
        s_ = roots[0] if roots else float("nan")          # no attractive solution: S and M fail
        one_list[contact] = dict(
            sigma=s_, t0=s_ / (2.0 * math.pi),
            sigma_at_boundary=_boundary.sigma(one, 1.0, contact, float(xi_bridge)),
            inside=_boundary.inside(one, contact, float(xi_bridge)))
    _absent = _boundary.sigma(_boundary.rows(d10, d11, n7, n26, bridge="absent"), 0.0, 1.0, float(xi_bridge))
    _labelled = _boundary.rows(d10, d11, n7, n26, bridge="boundary")   # the label as charge, at M_b
    _lab_roots = _boundary.anchored(_labelled, 1.0, float(xi_bridge))
    labelled_bridge = dict(
        sigma_at_boundary=_boundary.sigma(_labelled, 1.0, 1.0, float(xi_bridge)),
        t0=_lab_roots[0] / (2.0 * math.pi) if _lab_roots else float("nan"),
        w_need=(2.0 * math.pi - _absent)
        / (N_bridge * (a1_scalar - float(xi_bridge))))    # the weight at M_b the labelled bridge would need to pass S
    lambda_G_one = -(math.pi / 2.0) * (str_one - N0_st) / one_list[1.0]["sigma"]**2  # no G₂ scalars
    # Gate verdicts.  S and M are computed on the one list. U, C, T and
    # P name the tests that derive them, each with the case in which it
    # fails.  tests/common_mode/test_gravity_gates.py fails if a witness is missing,
    # and tests/audit/gate_verdicts.py recomputes all six.
    _gate_S = "PASS" if all(v["inside"] for v in one_list.values()) else "FAIL"
    _gate_M = "PASS" if all(0.6 < v["t0"] < 15.0 for v in one_list.values()) else "FAIL"
    booked_gates = {
        "S sign": "PASS" if min(ledger_ratio.values()) > 0 and Sigma_total_UV > 0 else "FAIL",
        "M magnitude": "PASS" if all(0.2 < r_ < 5.0 for r_ in ledger_ratio.values()) else "FAIL"}
    _cov = one_list[1.0]
    gravity_gates = {
        "S sign": dict(status=_gate_S, tier="derived", witness=("tests/common_mode/test_gravity_gates.py::test_sign_and_magnitude_are_computed_on_the_one_list", "tests/common_mode/test_horizon_amount.py::test_horizon_amount_the_gates_read_only_the_amount_past_the_searchlight", "tests/common_mode/test_common_mode_probes.py::test_composite_metric_tensor_term_is_the_vacuum_term", "tests/common_mode/test_common_mode_probes.py::test_ledger_kernel_rows_are_einstein_with_the_ledger_weights",), note=(f"the nameable rows of the one list give {_absent:+.3f}, repulsive. The "
            f"bridge vents as amount and makes the stiffness attractive, Sigma = "
            f"{_cov['sigma']:.3f}/{one_list[0.0]['sigma']:.3f}/{one_list[-1.0]['sigma']:.3f} "
            f"across the contacts, with its regulator inside the boundary "
            f"(tests/probes/horizon_amount.py)")),
        "U universality": dict(status="PASS", tier="derived", witness=("tests/common_mode/test_saddle_stabilizer.py::test_saddle_stabilizer_every_row_reads_one_cone", "tests/common_mode/test_saddle_stabilizer.py::test_saddle_stabilizer_a_frame_gives_each_row_its_own_speed", "tests/common_mode/test_two_lumps.py::test_two_lumps_attract_with_their_whole_energy_as_weight", "tests/substrate/test_substrate_probes.py::test_lump_in_empty_field_clock_inside_weight_is_energy",), note=("on the empty field's saddle every row couples to the same metric, "
            "and two lumps attract with their whole energy as the weight, "
            "binding energy included (tests/probes/saddle_stabilizer.py, two_lumps.py)")),
        "M magnitude": dict(status=_gate_M, tier="derived", witness=("tests/common_mode/test_gravity_gates.py::test_sign_and_magnitude_are_computed_on_the_one_list", "tests/common_mode/test_horizon_amount.py::test_horizon_amount_the_gates_read_only_the_amount_past_the_searchlight", "tests/common_mode/test_common_mode_probes.py::test_ledger_kernel_rows_are_einstein_with_the_ledger_weights", "tests/common_mode/test_common_mode_probes.py::test_regulator_free_ratios_higgs_is_the_dimension_two_ratio",), note=(f"the anchored regulator t0 = Sigma/(2 pi) is "
            f"{_cov['t0']:.3f}/{one_list[0.0]['t0']:.3f}/{one_list[-1.0]['t0']:.3f} M_b^-2 "
            f"across the contacts: the stiffness is Planckian (boundary.py)")),
        "C one cone": dict(status="PASS", tier="derived", witness=("tests/common_mode/test_saddle_stabilizer.py::test_saddle_stabilizer_every_row_reads_one_cone", "tests/common_mode/test_saddle_stabilizer.py::test_saddle_stabilizer_a_frame_gives_each_row_its_own_speed", "tests/common_mode/test_common_mode_probes.py::test_signature_of_the_saddle_is_the_carriers", "tests/common_mode/test_common_mode_probes.py::test_one_cone_relative_doublet_is_klein_gordon_on_its_own_cone", "tests/common_mode/test_common_mode_probes.py::test_relativistic_carrier_nr_limit_is_exact_and_one_cone", "tests/common_mode/test_common_mode_probes.py::test_carrier_sakharov_tensor_speed_needs_scale_separation", "tests/substrate/test_substrate_probes.py::test_lump_in_empty_field_clock_inside_weight_is_energy",), note=("the vacuum has no preferred frame, so every row shares one light cone "
            "and the induced tensor term propagates at c_T = c (tests/probes/saddle_stabilizer.py, "
            "carrier_sakharov.py)")),
        "T tensor carrier": dict(status="PASS", tier="derived", witness=("tests/common_mode/test_common_mode_probes.py::test_composite_metric_tensor_term_is_the_vacuum_term", "tests/common_mode/test_common_mode_probes.py::test_composite_pole_kernel_is_transverse_and_einstein", "tests/common_mode/test_common_mode_probes.py::test_ledger_kernel_rows_are_einstein_with_the_ledger_weights", "tests/common_mode/test_common_mode_probes.py::test_orbital_lock_has_no_invariant_term_and_no_isotropic_state", "tests/common_mode/test_common_mode_probes.py::test_determinant_vectors_closed_form_is_maxwell_on_the_composite_metric", "tests/common_mode/test_common_mode_probes.py::test_signature_of_the_saddle_is_the_carriers", "tests/common_mode/test_common_mode_probes.py::test_self_consistent_graviton_row", "tests/common_mode/test_induced_cubic.py::test_induced_cubic_term_is_einsteins_with_the_kernels_stiffness", "tests/common_mode/test_kernel_sign.py::test_kernel_sign_closed_form_and_limit", "tests/common_mode/test_induced_cubic_rows.py::test_weyl_and_vector_rows_induce_the_einstein_density_beyond_quadratic_order", "tests/common_mode/test_horizon_amount.py::test_horizon_amount_the_gates_read_only_the_amount_past_the_searchlight",), note=("the composite metric's one loop kernel is the linearized "
            "Einstein-Hilbert operator, 1/6 per real component, with the massless "
            "spin-2 pole (tests/probes/composite_metric.py, composite_pole.py)")),
        "P no preferred location": dict(status="PASS", tier="derived", witness=("tests/common_mode/test_horizon_amount.py::test_horizon_amount_the_gates_read_only_the_amount_past_the_searchlight", "tests/common_mode/test_two_lumps.py::test_two_lumps_attract_with_their_whole_energy_as_weight", "tests/common_mode/test_common_mode_probes.py::test_composite_pole_kernel_is_transverse_and_einstein", "tests/common_mode/test_common_mode_probes.py::test_log_class_stiffness_holds_and_extra_modes_sit_at_the_cutoff", "tests/substrate/test_substrate_probes.py::test_lump_in_empty_field_clock_inside_weight_is_energy",), note=("between conserved sources the kernel is the massless spin-2 pole "
            "alone, gamma = 1, and the stiffness depends on the rows' counts, not "
            "their masses. Read as massive, the amount passes S and M but fails here "
            "(tests/probes/composite_pole.py, log_class.py, horizon_amount.py A7)")),
    }
    print(f"    ledger ratio 6π/Sigma (full ledger, diagnostic, regulator dependent): "
          f"{ledger_ratio[1.0]:.3f} covariant, {ledger_ratio[0.0]:.3f} hard cutoff, "
          f"{ledger_ratio[-1.0]:.3f} opposite contact. No value of G is read from it")
    print(f"    one list (boundary.py): Sigma at the regulator M_b = "
          + "/".join(f"{v['sigma_at_boundary']:+.3f}" for v in one_list.values())
          + f". Anchored Sigma = " + "/".join(f"{v['sigma']:.3f}" for v in one_list.values())
          + f", t0 = " + "/".join(f"{v['t0']:.3f}" for v in one_list.values())
          + " M_b^-2 (vector contact +1 covariant / 0 hard cutoff / -1 opposite)")
    print(f"    GRAVITY LANE GATES (each PASS depends on the lane's three statements):")
    for name, g in gravity_gates.items():
        print(f"      {g['status']:4s} {g['tier']:8s} {name:16s} {g['note']}")
    from .root import G_PRED, COMMON_MODE
    print(f"    G itself: the anchor inversion m_e -> M_Pl (+22 ppm) read through the common mode, "
          f"G M_b^2 = 3 Q_l - 1 = {COMMON_MODE:.8f}: G = {G_PRED:.6e} (root.py)")

    # Cosmological density
    #
    # (1) The vacuum term.  In the traceless form of the composite
    #     metric, a constant in the potential never couples. The vacuum
    #     term does not gravitate, and Lambda is set by the state as an
    #     integration constant (tests/probes/lambda_lane.py).  In the
    #     full form, ΛG = −(π/2)·str(1)/Σ², which is of order one.  Both
    #     forms give w = −1. A confirmed evolving dark energy falsifies
    #     this.
    # (2) Jacobson's thermodynamic argument recovers Einstein dynamics
    #     but leaves Lambda as an integration constant. It does not
    #     determine the observed vacuum density.
    # (3) CKN saturation is an outside bound, not a derivation.
    #
    # H₀ is a comparison scale.  This code evaluates the
    # resulting relation.
    from .root import RHO_CRIT_MAPPED as rho_crit_mapped
    rho_naive = M_Pl_GeV**4                  # naive QFT: ρ ~ M_Pl⁴
    print(f"    Friedmann: ρ_crit = 3H₀²/(8πG) = (3/8π) M_Pl² H₀²/common_mode "
          f"= {rho_crit_mapped:.2e} GeV⁴")
    print(f"       vs naive ρ ~ M_Pl⁴ = {rho_naive:.1e} GeV⁴  [scale comparison]")
    print(f"       full form ratio ΛG = −(π/2)·str(1)/Σ² = {lambda_G:+.3f} on the full "
          f"ledger, {lambda_G_one:+.3f} on the one list. The lane uses the traceless form, "
          f"where the vacuum term does not gravitate and Lambda is set by the state")
    print(f"       Jacobson leaves Lambda free. CKN saturation is an outside bound")

    # ══════════════════════════════════════════════════════════════════
    #  HIGGS + EW: SM RGE transport, M_W, M_Z, import audit (higgs_ew.py)
    # ══════════════════════════════════════════════════════════════════

    higgs_ew_data = derive_higgs_ew(masses, couplings)

    # ══════════════════════════════════════════════════════════════════
    #  BARYOGENESIS: eta_B = -n₇ * J_lep * exp(-d₁₀²pi²/2)
    # ══════════════════════════════════════════════════════════════════

    print(f"\n  --- BARYOGENESIS: n₇, J_lep, and the G₂ face action ---")

    # Assembly: each factor maps to a Sakharov condition. The action is
    # the depth 4 G₂ face reading derived in root.py.
    #
    # 1. B−L violation: G₂ instanton tunnelling
    #    S_inst = h∨(G₂)·π²/2 = d₁₀²·π²/2 = 2π²
    #    The h-dual bridge metric transports the F₄-face response to
    #    the G₂ face: b_face=(32/3)(9/4)=24=h∨(E8)-C2(26), so
    #    S_inst=2pi/(b_face alpha_G2)=(4/9)S_lepton=2pi^2.
    #    Tunnelling probability: exp(−S_inst) = exp(−2π²) ≈ 2.68×10⁻⁹
    #
    # 2. CP violation: leptonic Jarlskog invariant J_lep
    #    Derived from the conjugation modular invariant of SU(3)₃
    #    with the charged lepton rotation (see mixing.py).
    #
    # 3. Multiplicity: n₇ = dim(G₂ fund) = 7
    #    Seven directed Fano plane lines, each an independent asymmetry
    #    channel.  The Fano plane is the incidence geometry of the
    #    octonion multiplication table (G₂ = Aut(O)).
    #
    # Assembly: η_B = −n₇ · J_lep · exp(−d₁₀²π²/2)
    # The minus sign is the sphaleron conversion B = −(28/79)L: the
    # lepton sector CP asymmetry (J_lep < 0) converts to a positive
    # baryon excess.
    #
    from .root import S_baryo
    S_inst = S_baryo  # = 2π², the G₂ face action (root.py)
    tunnelling = math.exp(-S_inst)
    J_lep = mixing['J_lep']

    # JOINT COSMOLOGY: eta_B here and Omega_DM/Omega_b
    # (dark_sector.py) share one engine input, R(h) = h/2π.
    # Omega_b h² = eta_10/273.9 = 0.02256 (+1.2σ), and carrying that
    # through 2π−1 gives Omega_DM h² = 0.1192 (−0.7σ vs Planck
    # 0.1200(12)). The BBN conversion eta_10 = 273.9 Omega_b h² is a
    # cited change of variables applied to a printed number. No step
    # of the chain uses it. The pair is falsified if Planck class data
    # fix both densities at values whose ratio and product
    # (eta_B, 2π−1) cannot match.
    # In the traceless form Lambda is set by the state. root.py
    # computes the Friedmann critical density from H0, and the current
    # dark energy density uses the measured Omega_Lambda (dark_sector.py).
    from .root import PDG_COSMO
    eta_B = float(n7) * (-J_lep) * tunnelling
    eta_B_obs = PDG_COSMO['eta_B']

    ratio_B = eta_B / eta_B_obs
    pull_pct_B = abs(1 - ratio_B) * 100
    pull_sigma_B = (eta_B - eta_B_obs) / PDG_COSMO['eta_B_err']

    print(f"    S_inst = d₁₀²pi²/2 = {d10}²pi²/2 = {S_inst:.4f}")
    print(f"    exp(-S) = exp(-{d10}²pi²/2) = {tunnelling:.6e}")
    print(f"    J_lep = {J_lep:+.6f}")
    print(f"    eta_B = {n7} * {-J_lep:.6f} * {tunnelling:.4e}")
    print(f"          = {eta_B:.4e}")
    print(f"    Planck: {eta_B_obs:.2e}  ({(ratio_B - 1) * 100:+.2f}%, "
          f"{pull_sigma_B:+.2f}σ)  [n₇ Fano lines, J_lep from the PMNS map]")

    # ══════════════════════════════════════════════════════════════════
    #  NEUTRINOS: rank bound from the F₄ singlet identification
    # ══════════════════════════════════════════════════════════════════

    print(f"\n  --- NEUTRINOS: F₄ singlet structure ---")

    # n₂₆ → 8_v + 8_s + 8_c + (n₂₆ − 3×8) singlets
    # The framework identifies the (n₂₆ − 24) = 2 singlets with RH
    # neutrinos. In a minimal type-I seesaw with no additional Weinberg
    # operator, type II/III contribution, or other mass term, two RH
    # neutrinos imply rank(M_nu) <= 2. Rank exactly 2 needs an additional
    # full rank Yukawa assumption.
    n_singlets = n26 - 3 * (d11**2 - 1)   # = 26 - 24 = 2  (d₁₁²−1 = 8)
    print(f"    {n26} -> 8_v + 8_s + 8_c + {n_singlets}*1")
    print(f"    Assumption: minimal type-I seesaw. No other neutrino mass terms")
    print(f"    {n_singlets} identified RH neutrinos -> rank(M_nu) <= {n_singlets}")
    print(f"    Consequence: at least one massless light neutrino")

    # Structural rank bound only. No external Δm² input is needed for it.
    # The absolute neutrino mass scale and ordering are not selected by the
    # algebraic framework. PROVENANCE of the normal ordering benchmark
    # (Σm_ν = 58.78 meV, m_β = 8.81 meV, m_ββ ∈ [1.51, 3.70] meV):
    # it assumes normal ordering and m₁ = 0, and combines that
    # assumption with the MEASURED splittings
    # (reference.NUFIT_OSCILLATION, reference data the engine does not
    # read). It sets m₂ = √Δm²₂₁ and m₃ = √Δm²₃₁, then forms the
    # standard Σ/β/ββ combinations with NuFit angles.  The frozen values in
    # tests/common_mode/test_common_mode_scorecard.py are conditional benchmark arithmetic.
    #
    # WHAT FALSIFIES THIS. Three established nonzero light neutrino
    # masses falsify the rank bound from two RH neutrinos. The ordering
    # alone does not.
    print(f"    Ordering: either (NO: m₁ = 0, IO: m₃ = 0)")

    result = {
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
        'Sigma_bridge': Sigma_bridge,
        'Sigma_total_UV': Sigma_total_UV,
        'Sigma_total_broken': Sigma_total_broken,
        'G_ratio_UV': G_ratio_UV,
        'G_ratio_broken': G_ratio_broken,
        'G_ratio_canonical': G_ratio_canonical,   # ledger diagnostic, covariant regulator
        'ledger_ratio': ledger_ratio,
        'gravity_gates': gravity_gates,
        'one_list': one_list,
        'labelled_bridge': labelled_bridge,
        'booked_gates': booked_gates,
        'rho_crit_mapped': rho_crit_mapped,
        # Baryogenesis
        'eta_B': eta_B,
        'eta_B_obs': eta_B_obs,
        'ratio_B': ratio_B,
        'pull_sigma_B': pull_sigma_B,
        'S_inst': S_inst,
        'tunnelling': tunnelling,
        'J_lep': J_lep,
        'dim_G2_fund': n7,
        # Neutrinos (structural rank bound, no Δm² input)
        'ordering': 'either',
        'n_RH': n_singlets,
        'light_mass_rank_max': n_singlets,
        'massless_light_states_min': d11 - n_singlets,
        'massless_state_by_ordering': {'normal': 'm1', 'inverted': 'm3'},
        'neutrino_mass_mechanism_assumption': 'minimal type-I seesaw only',
    }
    result.update(higgs_ew_data)
    return result
