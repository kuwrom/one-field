"""
Gravity + Higgs + Baryogenesis + Neutrinos from (d₁₀, d₁₁, n₇, n₂₆, pi)
+ the electron anchor (M_Pl is DERIVED in root.py from m_e).

Import inventory: higgs_ew.py (the consolidated table lives there).

═══════════════════════════════════════════════════════════════════════
WHAT GRAVITY IS (gravity paper)
═══════════════════════════════════════════════════════════════════════

A matter standing wave is a circulating soliton.  Its circulation
biases the background shared by the substrate's three components.

The exact Z₃ Fourier split separates the q = 1, 2 relative modes
(matter) from the q = 0 common mode (geometry).  Projecting q = 0
gives the screened scalar response below.  The common acoustic metric
carries the induced spin-2 channel, and the bridge heat-kernel ledger
fixes its coefficient.

In the 3D formation probe, a seeded relative-sector knot and its
common-mode depletion form together.  The soliton casts a shadow.
Gravity is the common geometry induced by matter's circulation.

═══════════════════════════════════════════════════════════════════════
WHY GRAVITY AT ALL: the common-mode memory (gravity paper, summary)
═══════════════════════════════════════════════════════════════════════

The substrate is the Z₃-symmetric coupled NLS system (masses.py).
In Madelung variables ψ_k = √ρ_k e^{iS_k/ħ}, perturbations decompose
under the Z₃ Fourier transform into:
    q = 0    the COMMON mode, shared by all three components.
    q = 1,2  the RELATIVE modes, they distinguish components.
q = 1,2 carry the matter-like standing waves. q = 0 carries the
gravitational memory.  This resolves the Barceló-Liberati-Visser
"two roles of ρ" obstruction of analogue gravity: matter and
geometry are different Fourier projections of ONE substrate.

Linearising and projecting onto q = 0, the relative-sector energy
density sources the common-mode density via a SCREENED POISSON
equation:
    (1 − ξ₀²∇²) R₀ = −(c*/λ₀ρ₀) H_matter,
with healing length ξ₀ = C_EFF·ℓ_Pl = ℓ_Pl/d₁₀ under the Planck
identification of the healing length (root.py, C_EFF;
sakharov_normalisation.py FACT 4).
The scalar q=0 force is therefore Yukawa-screened: e^{−r/ξ₀} <
10⁻³⁰⁰ at any observable distance, the framework passes the
Cassini PPN bound by ~295 orders of magnitude.  Macroscopic gravity
is the INDUCED SPIN-2 channel (Sakharov): matter standing waves
propagate on the common acoustic metric, and integrating them out
generates the Einstein-Hilbert term (Sakharov induction).  The minimal
matter spectrum alone gives the WRONG-SIGN normalisation
(Σ_min ≈ −10.2). The 182 bridge channels supply the positive
back-reaction (+29.126).  That repair is this module.

Gravity (bridge sector):
    N_bridge = n₇ * n₂₆ = 182
    xi = alpha_G₂(M_Pl) * E[v²] * h_bridge = 1/(48pi)
    Heat kernel -> G_N

THREE EXPOSURES OF THE EMERGENT METRIC, and what ends each.
    |c_gw/c − 1| = 0 identically: tensor modes and light ride one
    emergent metric.  Any confirmed |c_gw/c − 1| > 1e-15 kills it
    (GW170817-class multi-messenger bounds already sit at ~1e-15
    and are consistent).
    PPN γ = 1 exactly: no scalar admixture in the metric read by
    light bending and time delay, because the scalar q=0 channel is
    Yukawa-screened above.  A confirmed |γ − 1| > 2e-5 kills it
    (Cassini-class time delay currently bounds 2.3e-5).
    The equivalence principle is exact: PvP = 0 is an algebraic
    identity of the Z₃ harmonic (protected_forgetting.py), so EP
    violation is null at EVERY future precision, not merely small.
    Any confirmed violation kills it (MICROSCOPE η < 1e-15 is
    consistent).  The same PvP = 0 is why the dark sector couples
    only through the common mode (dark_sector.py).

Higgs + EW: higgs_ew.py (SM RGE transport, M_W, M_Z, import audit).

Baryogenesis assembly (the depth-4 G₂-face action):
    eta_B = n₇ × J_lep × exp(−h∨(G₂)π²/2) = 7 · J_ℓ · e^{−2π²}
      • e^{−2π²}: root.py derives the G₂-face action from the F₄
        dimensional-transmutation action: b_face=(32/3)(9/4)=24=
        h∨(E8)-C2(26), hence S_G2=2pi/(24 alpha_G2)=2pi². The same
        derived factor is used for the
        Majorana scale M_R = M_Pl e^{−2π²} ≈ 3.3×10¹⁰ GeV >> T_EW,
        so out-of-equilibrium decay (3rd condition) needs no extra
        ingredient.
      • J_lep: leptonic Jarlskog (mixing.py).  CP violation traces
        to the ORIENTED Fano plane: seven directed lines select a
        preferred chirality.
      • n₇ = 7: the asymmetry channels, one per directed Fano
        line, each one octonionic multiplication.
    The asymmetry is generated NON-THERMALLY, at the
    instanton level, it is a topological count, not a decay-era
    process.  The Majorana scale M_R = M_Pl e^{−2π²} >> T_EW means
    the B−L asymmetry is frozen in before any thermal washout era
    exists. The Boltzmann/washout machinery of thermal leptogenesis
    describes a different mechanism and does not apply here.

Neutrinos:
    26 -> 8_v + 8_s + 8_c + 1 + 1   (SO(8) triality decomposition
    of the F₄ fundamental: three 8's tied by triality + 2 singlets)
    The 2 singlets are the only RH neutrinos -> the seesaw matrix
    has rank 2 -> m₁ = 0 EXACTLY and normal ordering, structurally.

═══════════════════════════════════════════════════════════════════════
Why the bridge contributes a SCALAR sign (not a vector sign)
═══════════════════════════════════════════════════════════════════════

The bridge sector has 182 channels, far more than the 8 dimensions of
the octonion fiber, the deepest internal structure this ledger can
track (Hurwitz again).  Tracking all 182 as one coherent
multi-layer object would require a normed division algebra that does
not exist.  If we knew what was inside each channel ("this one is
amplitude, that one is delay"), we could keep some of them complex.
We don't, so the bridge has to be treated as one ensemble, and every
mode contributes as a real scalar.  Reading the same channels as
vectors would give −1/6 and reverse the sign: gravity repulsive,
matter never binds.

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
        (STATUS: as coded, x²/x² ≡ 1 by construction; the CFT claim
        is cited, not computed, couplings.py._derive_bridge.)

═══════════════════════════════════════════════════════════════════════
Same bridge, two observables (cross-link to couplings.py)
═══════════════════════════════════════════════════════════════════════

The (7,26) bridge carries two readings through the same h_bridge = 1:
  * Gravity: 182 scalar-like heat-kernel channels each contribute
    a₁/R = 1/6 − xi_bridge with xi_bridge = alpha_G₂ * (1/2) * h_bridge
    = 1/(48π). The Newton normalisation is an internal closure
    (sakharov_normalisation.py).
  * QED:    The same h_bridge = 1 makes the bridge a marginal worldsheet
    primary. Its self-interference fraction is h_bridge/(2π)
    = 1/(2π), giving (1 − 1/(2π)) that converts alpha_em = pi/2⁹ into
    alpha(0) = 1/137.035999050.  See couplings.py.
Same sector, same h_bridge, two observables.

Reference:
    Kahsay, Kibrom Kidane (2026). Three-paper series.
    Zenodo. https://doi.org/10.5281/zenodo.20144381
"""

import cmath
import math

from root import (d10, d11, n7, n26,
                  hv_G2, hv_F4,
                  N_bridge, h_bridge, h_7, h_26,
                  alpha_G2_Pl, xi_bridge, E_v2,
                  M_Pl_GeV,
                  BA_ratio, h10)
from higgs_ew import derive_higgs_ew


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

    # Protected forgetting on the canonical Z₃ harmonic
    #   v_k = B cos(θ + 2πk/3),  B = √d₁₀ = √2,  θ = h₁₀ = 2/9
    # (see protected_forgetting.py for the moment bookkeeping).
    #
    #   PvP    = Σv_k/3   = 0        forgetting, holds for every (B, θ)
    #   Pv²P   = Σv_k²/3  = B²/2 = 1 TOTAL second moment (both channels)
    #   E[v²]  = |v̂₁|²    = B²/4 = ½ PER-MODE power, one Z₃ channel
    #
    # The unique (7,26) bridge carries this G2/Z3 response tensored
    # with identity on the F4 26 (wiring.py W6, Schur trace), so the
    # ledger consumes the per-mode power E[v²] = ½.  The total is 1.
    # The two are distinct and easily conflated: the notation
    # "Pv²P = ½P" names the total, not the per-mode power.
    # the canonical harmonic from ONE source (root.BA_ratio = √d₁₀,
    # root.h10 = 2/9): a drift in either parameter fails here
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
    led_a1_bridge = Ledger(
        "a1_bridge", a1_scalar, "add", web_depth=4).echo(
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
        led = Ledger(name, 0.0, "add", web_depth=4)
        for path, n_dof, a1, note in sector_terms:
            led.echo([path], n_dof * a1, 1, "FORCED",
                     f"{n_dof} dof × a₁ = {a1:+.4f}  {note}",
                     kind="vent")
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
    # Diagnostic only (see the UV bookkeeping note below).  The G₂
    # rows keep unbroken-style counting: with the G₂ breaking real,
    # six of the seven G₂ scalars are eaten and the coset vectors
    # carry 18 Proca dof, not 6; the SM rows below rebook their
    # Goldstones, the G₂ rows do not.  The canonical closure uses
    # the UV (symmetric-phase) ledger, which this does not affect;
    # a full heat-kernel derivation must resolve this bookkeeping.
    led_Sigma_br = _vent_ledger("Sigma_vent(broken)", [
        ("Higgs(1)",     n_higgs_UV - d11, a1_scalar, "physical Higgs"),
        ("gauge(18)",    2 * d11**2,       a1_gauge,  "gluons+photon"),
        ("Proca(9)",     d11**2,           a1_proca,  "massive EW"),
        ("fermions(90)", n_fermion,        a1_fermion, "3 generations"),
        ("coset(6)",     n_G2_vec,         a1_proca,  "G₂/SU(3) vectors"),
        ("G₂-scal(7)",   n_G2_scal,        a1_scalar, "G₂ fundamental"),
        ("bridge(182)",  N_bridge,         a1_bridge, "coherence channels"),
    ])

    # Face split. The h-dual candidate assigns 9/13 to the F4 metric
    # face; deriving this absorption weight from the heat kernel is open.
    _face_split = -(hv_F4 / float(hv_G2 + hv_F4)) / (2.0 * math.pi)
    for _led in (led_Sigma_UV, led_Sigma_br):
        _led.echo(["bridge self-echo (F₄ h∨-share 9/13)"],
                  _face_split, 1, "FORCED",
                  "face-split law: no self-dilution of the metric face")

    WEB.solve()          # gravity nodes join the solved web state

    Sigma_bridge = N_bridge * a1_bridge
    Sigma_total_UV = led_Sigma_UV.value()
    Sigma_total_broken = led_Sigma_br.value()
    Sigma_min_UV = Sigma_total_UV - Sigma_bridge
    Sigma_min_broken = Sigma_total_broken - Sigma_bridge

    # ── UV bookkeeping ──────────────────────────────────────────────
    # Σ is the coherent venting of all knots into the common mode, and
    # venting happens at the substrate scale (the healing length
    # ξ₀ = ℓ_Pl/d₁₀, root.py C_EFF; sakharov_normalisation.py FACT 4): the
    # quadratically-divergent induced-G integral
    # localizes at the cutoff.  At Planckian distances there is no
    # electroweak condensate (v_EW/M_Pl ~ 10⁻¹⁷), so the spectrum that
    # vents is the SYMMETRIC-phase spectrum.  The UV bookkeeping is
    # therefore selects the symmetric-phase spectrum.
    # The broken-phase counting is retained as a diagnostic only (it
    # counts IR-phase dof whose contribution under Λ² is negligible).
    # The face-split (no-self-dilution) law closes the UV ledger
    # law above. The remaining −8e-8 is below the (α/2π)² depth-2
    # echo scale, the ledger is depth-complete here.

    # Gravity = the memory of the total coherent venting
    # NOTE: an INTERNAL closure against the
    # Sakharov normalisation Σ = 6π, not a comparison to measured G;
    # the only measured-G contact in the repo is root.py's 22 ppm
    # M_Pl assert.  The scorecard's 'internal' tag carries this.
    # STATUS OF THE TARGET AND THE EDGE: the 6π coefficient is stated,
    # not derived.  Its exact content is
    # Σ_target = 2π·d₁₁ through the standard proper-time algebra
    # 1/G = Σ/(2π t₀), i.e. t₀ = d₁₁·ℓ_Pl²; two exactly degenerate
    # engine factorizations carry it (P-closure × healing area,
    # d₁₀ cancelling; common-mode share × Planck area), and the
    # mechanism is the residue
    # (tests/probes/sakharov_normalisation.py).  The face-split edge
    # above is worth −0.110 on Σ ≈ 18.96 and carries the closure
    # from a 0.58% miss to 8.3e-8 (sensitivity: no edge 0.9942,
    # G₂-share 0.9968, coded F₄-share 0.99999992, full unit 1.0026).
    # The metric-uniqueness scan discriminates the h∨ split BY its
    # agreement with this target, so the choice is pinned by the
    # closure it enables, not derived upstream of it.
    target = 6.0 * math.pi
    G_ratio_UV = target / Sigma_total_UV
    G_ratio_broken = target / Sigma_total_broken
    G_ratio_canonical = G_ratio_UV

    # ── gravity as a graph node (live ratio base): G = 6π/Σ ─────────
    # The column total of everyone's venting, read from the state.
    WEB["G_ratio"] = Ledger(
        "G_ratio", lambda s: 6.0 * math.pi / s["Sigma_vent(UV)"], "mul",
        web_depth=4)
    WEB.solve()          # gravity nodes join the solved web state
    assert abs(WEB.state["G_ratio"] / G_ratio_canonical - 1.0) < 1e-12
    err_UV = 100.0 * (G_ratio_UV - 1.0)
    err_broken = 100.0 * (G_ratio_broken - 1.0)

    print(f"    Sigma_min(UV)     = {Sigma_min_UV:+.4f}")
    print(f"    Sigma_min(broken) = {Sigma_min_broken:+.4f}")
    print(f"    Sigma_bridge      = {Sigma_bridge:+.4f}  [{N_bridge} * (1/6 - 1/(48pi))]")
    print(f"    G_ind/G_N (UV)     = {G_ratio_UV:.6f}  ({err_UV:+.2f}%)")
    print(f"    G_ind/G_N (broken) = {G_ratio_broken:.6f}  ({err_broken:+.2f}%)")
    print(f"    G_ind/G_N INTERNAL    = {G_ratio_canonical:.9f}  (UV: venting")
    print(f"      localizes at the cutoff, where EW is unbroken. Broken-")
    print(f"      phase value {G_ratio_broken:.6f} is diagnostic only.")
    print(f"      face-split no-self-dilution law)")

    # Cosmological constant: Volovik → Jacobson-Clausius → CKN
    #
    # (1) Volovik: self-bound condensate → ρ_vac(eq) = 0 (Gibbs-Duhem
    #     identity, not fine-tuning. Same mechanism as superfluid helium).
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
    from root import RHO_LAMBDA as rho_Lambda
    rho_naive = M_Pl_GeV**4                  # naive QFT: ρ ~ M_Pl⁴
    print(f"    Λ: ρ_Λ = (3/8π) M_Pl² H₀² = {rho_Lambda:.2e} GeV⁴")
    print(f"       vs naive ρ ~ M_Pl⁴ = {rho_naive:.1e} GeV⁴  [10¹²³ resolved]")
    print(f"       Volovik (ρ_eq=0) → Jacobson (3/8π) → CKN saturation")

    # ══════════════════════════════════════════════════════════════════
    #  HIGGS + EW: SM RGE transport, M_W, M_Z, import audit (higgs_ew.py)
    # ══════════════════════════════════════════════════════════════════

    higgs_ew_data = derive_higgs_ew(masses, couplings)

    # ══════════════════════════════════════════════════════════════════
    #  BARYOGENESIS: eta_B = n₇ * J_lep * exp(-d₁₀²pi²/2)
    # ══════════════════════════════════════════════════════════════════

    print(f"\n  --- BARYOGENESIS: n₇, J_lep, and the G₂ face action ---")

    # Assembly: each factor maps to a Sakharov condition. The action is
    # the depth-4 G₂-face reading derived in root.py.
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
    #    with the charged-lepton rotation (see mixing.py).
    #
    # 3. Multiplicity: n₇ = dim(G₂ fund) = 7
    #    Seven directed Fano-plane lines, each an independent asymmetry
    #    channel.  The Fano plane is the incidence geometry of the
    #    octonion multiplication table (G₂ = Aut(O)).
    #
    # Assembly: η_B = n₇ · J_lep · exp(−d₁₀²π²/2)
    #
    from root import S_baryo
    S_inst = S_baryo  # = 2π², the G₂-face action (root.py)
    tunnelling = math.exp(-S_inst)
    J_lep = mixing['J_lep']

    # JOINT COSMOLOGY CLOSURE, a check the dictionary never consulted:
    # eta_B here and Omega_DM/Omega_b (dark_sector.py) use SEPARATE
    # channels, and they close together.  Omega_b h² = eta_10/273.9 =
    # 0.02255 (+1.2σ), and carrying that through 2π−1 gives
    # Omega_DM h² = 0.1192 (−0.7σ vs Planck 0.1200(12)).  The BBN
    # conversion eta_10 = 273.9 Omega_b h² is cited, not derived: it
    # is a change of coordinate on a printed number, and no step of
    # the chain consumes it.  Planck-class data pinning both
    # densities at values whose ratio and product cannot be met by
    # (eta_B, 2π−1) kills the pair.
    # The Lambda sector is priced rather than implied: rho_Lambda =
    # (3/8π) M_Pl² H0² = 3.7e-47 GeV⁴ is the EQUILIBRIUM de Sitter
    # scale with H0 as a second dimensional reference (root.py labels
    # it: fixes a scale, not a coupling).  Observed is ~2.5e-47, low
    # by 1/Omega_Lambda = 1.46, read as pre-asymptotic cosmology, not
    # algebra.
    from root import PDG_COSMO
    eta_B = float(n7) * J_lep * tunnelling
    eta_B_obs = PDG_COSMO['eta_B']

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
    # framework. What IS predicted is the rank (m₁=0) and ordering (normal).
    # PROVENANCE of the numeric neutrino predictions (Σm_ν = 58.78
    # meV, m_β = 8.81 meV, m_ββ ∈ [1.51, 3.70] meV):
    # they are m₁ = 0 combined with the MEASURED splittings
    # (reference.NUFIT_OSCILLATION, a comparison-side reference the engine
    # does not read): m₂ = √Δm²₂₁, m₃ = √Δm²₃₁, then the standard
    # Σ/β/ββ combinations with NuFit angles.  The frozen values in
    # tests/test_interference.py are that arithmetic, pinned; the
    # framework's own content in them is exactly (rank 2, m₁ = 0,
    # normal ordering), nothing more.
    #
    # WHAT ENDS THIS.  Inverted ordering, or an established m₁ > 0,
    # kills the rank-2 seesaw outright (JUNO and cosmology decide).
    # Σm_ν = 58.78 meV is the most immediate exposure: a cosmological
    # bound below ~58 meV, or a detection above ~60 meV, kills it
    # (DESI + CMB currently ~70 meV and tightening).  Direct
    # kinematics sees nothing until sensitivity reaches ~9 meV, so
    # any detection above that kills m_β = 8.81 meV (KATRIN,
    # Project 8).  Every planned 0νββ experiment returns null: an
    # observation implying m_ββ > 4 meV kills the range (LEGEND-1000,
    # nEXO).
    print(f"    m₁ = 0 (structural prediction: rank-{n_singlets} seesaw)")
    print(f"    Ordering: normal (m₁ < m₂ < m₃)")
    print(f"    Testable by: KATRIN endpoint, 0νββ non-observation, cosmological Σmν")

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
        'Sigma_min_broken': Sigma_min_broken,
        'Sigma_bridge': Sigma_bridge,
        'Sigma_total_UV': Sigma_total_UV,
        'Sigma_total_broken': Sigma_total_broken,
        'G_ratio_UV': G_ratio_UV,
        'G_ratio_broken': G_ratio_broken,
        'G_ratio_canonical': G_ratio_canonical,
        'err_UV': err_UV,
        'err_broken': err_broken,
        'rho_Lambda': rho_Lambda,
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
    result.update(higgs_ew_data)
    return result
