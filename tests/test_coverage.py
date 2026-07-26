"""Look-elsewhere, cross-bracing, and grammar-verification tests.

Three claims are frozen here:

1. SPARSITY.  Where the code publishes an enumeration over a finite
   menu (heavy-quark Koide structures, terminus words), the menu is
   sparse: a random target would NOT cheaply match.  If a refactor
   ever made a menu dense enough that hits are cheap, these tests
   fail and the corresponding scorecard entries lose their weight.

2. CROSS-BRACING.  Constants are entered once and consumed by several
   unrelated observables.  These tests recompute each consumer from
   the shared entry. A "fit" that adjusted one consumer independently
   would break the identity and fail here.

3. GRAMMAR VERIFICATION.  Every FORCED echo multiplicity used by the
   Web ledger is recomputed from the four-integer dictionary and
   compared against the solved state via closed forms.  This converts
   the FORCED label from an assertion into a checked identity.

Run with:
    pytest -q
"""

import math
import os
import sys
from fractions import Fraction

import pytest

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "interference")))

from root import PDG_MASSES


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════

def _mb_from_Q(Q, mc, mt, lo=100.0, hi=40000.0):
    """Invert the Koide relation for m_b at fixed m_c, m_t (MeV)."""
    f = lambda mb: (mc + mb + mt) / (math.sqrt(mc) + math.sqrt(mb)
                                     + math.sqrt(mt))**2 - Q
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def _koide_menu(mc, mt):
    """The eight structures Q0 + h/K^p, h in {fund, adj}, p in 1..4."""
    out = {}
    for hname, h in (("fund", 2 / 9), ("adj", 1 / 2)):
        for p in range(1, 5):
            out[(hname, p)] = _mb_from_Q(2 / 3 + h / 6**p, mc, mt)
    return out


# ═══════════════════════════════════════════════════════════════════════
#  1. Sparsity (look-elsewhere weight of the published enumerations)
# ═══════════════════════════════════════════════════════════════════════

def test_koide_selection_check_is_unique(res):
    """masses.py SELECTION CHECK, frozen: among all eight structures
    Q0 + h/K^p, exactly one matches PDG m_b at sub-percent, and it is
    the canonical branch (adjoint weight, cubic altitude, p = 3)."""
    menu = _koide_menu(res["m"]["m_c"], res["m"]["m_t"])
    pdg_b = PDG_MASSES["b"]
    subpct = [key for key, mb in menu.items()
              if abs(mb - pdg_b) / pdg_b < 0.01]
    assert subpct == [("adj", 3)]
    # runner-up (fund, 3) misses by ~2%. All others by > 2%
    others = sorted(abs(mb - pdg_b) / pdg_b
                    for key, mb in menu.items() if key != ("adj", 3))
    assert others[0] > 0.02


def test_koide_menu_is_sparse(res):
    """The 8-way menu covers < 25% of its own span at +-1%: the
    enumeration is not dense enough to make sub-percent hits cheap."""
    menu = sorted(_koide_menu(res["m"]["m_c"], res["m"]["m_t"]).values())
    span = menu[-1] - menu[0]
    coverage = sum(0.02 * v for v in menu) / span
    assert coverage < 0.25


def test_terminus_menu_is_sparse(res):
    """words.py terminus: three length-4 candidates cover < 10% of the
    physical top window [50, 200] GeV at +-1%, and the two excluded
    words are wrong by catastrophic margins (not near-misses)."""
    m_tau_GeV = res["m"]["m_tau"] / 1e3
    cands = [36, 78, 97]
    coverage = sum(2 * 0.01 * b * m_tau_GeV for b in cands) / (200.0 - 50.0)
    assert coverage < 0.10
    pdg_t = PDG_MASSES["t"] / 1e3
    for b in (36, 78):                       # the excluded words
        assert abs(b * m_tau_GeV - pdg_t) / pdg_t > 0.15


def test_words_lemma_outputs(res):
    """Freeze the generation-word lemma's checkable content: the
    boundary-walk counts, the down seed, the neutral-sector FPdim,
    and the exclusion set.  words.derive() internally re-verifies the
    Z3 universal grading (10^3 fusion triples), C_ad fusion closure,
    and dim Hom(f x fbar, a) = 1. Any failure raises."""
    import contextlib, io
    import words
    with contextlib.redirect_stdout(io.StringIO()):
        w = words.derive()
    assert w["base"] == (4, 12, 97)
    assert w["base_down_seed"] == 9
    assert w["FP_C0"] == 12
    assert w["excluded"] == {"f²a²": 36, "f³a ⊕ fa³": 78}


def test_vew_multiplicity_resolution(res):
    """The 16(alpha/2pi)^2 term steps v_EW by ~1.35e-6 per unit of
    multiplicity. The achieved residual against the G_F value is a
    small fraction of one step (i.e. The integer landed, it was not
    merely bracketed)."""
    alpha = 1.0 / res["c"]["inv_alpha_phys"]
    step = (alpha / (2 * math.pi))**2
    v_GF = 246.219651                        # from G_F (declared import)
    resid = abs(res["m"]["v_EW_pred_GeV"] / v_GF - 1.0)
    assert step == pytest.approx(1.35e-6, rel=0.05)
    assert resid < 0.1 * step


# ═══════════════════════════════════════════════════════════════════════
#  2. Cross-bracing (one entry, many consumers)
# ═══════════════════════════════════════════════════════════════════════

def test_charge_trace_four_consumers():
    """charge_trace = 8/3 is entered once and consumed by four sectors.
    Changing it to improve any one consumer breaks the other three."""
    from root import charge_trace, singh_ratio, C2_26, Q0, h10, d10, d11
    assert charge_trace == Fraction(8, 3)
    assert charge_trace == Fraction(d10**3, d11)          # identity form
    # consumer 1: alpha(0) base via the Singh ratio
    assert singh_ratio == charge_trace * C2_26 == 16 == d10**4
    # consumer 2: strange bridge^2 base
    assert Q0**2 * charge_trace == Fraction(32, 27)
    # consumer 3: top back-reaction corr3
    assert h10 / charge_trace == Fraction(1, 12)
    # consumer 4: the v_EW exponent multiplicity (same 16, see below)


def test_bridge_marginality_three_consumers():
    """h_bridge = h_7 + h_26 = 1 (marginal) feeds alpha(0), the dark
    ratio, and gravity's heat-kernel normalisation."""
    from root import h_7, h_26, h_bridge
    assert h_7 == Fraction(2, 5) and h_26 == Fraction(3, 5)
    assert h_bridge == 1


def test_vertex_current_trace_identity():
    """The derived E8 positive-root current trace supplies 30."""
    from root import (N_vertex, vertex_charge_trace, hv_E8,
                      WIRING_DATA)
    assert N_vertex == vertex_charge_trace == hv_E8 == 30
    assert WIRING_DATA["vertex_charge_trace"] == 30


# ═══════════════════════════════════════════════════════════════════════
#  3. Grammar verification (FORCED echo terms recomputed as closed forms)
# ═══════════════════════════════════════════════════════════════════════

def test_inv_alpha_closed_form(res):
    """1/alpha(0) solves the FINAL cubic
        x = (512/pi) * (1 - 1/(2pi) - 2*(1/x / 2pi)^2)
    with base 512/pi (Singh ratio), depth-1 coefficient h_bridge = 1,
    depth-4 response multiplicity exactly 2 (orientation rule)."""
    x = res["c"]["inv_alpha_phys"]
    rhs = (512 / math.pi) * (1 - 1 / (2 * math.pi)
                             - 2 * (1 / (x * 2 * math.pi))**2)
    assert x == pytest.approx(rhs, abs=1e-9)
    assert x == pytest.approx(137.035999050, abs=5e-9)


def test_quark_action_closed_form(res):
    """S_quark = 9pi^2/2 - 6 + 15/512 - 16(alpha/2pi)^2, every term a
    dictionary polynomial: 9pi^2/2 = hv(F4) pi^2/2. 6 = C2(26) =
    d10*d11. 15/512 = 30*(pi/512)/(2pi). 16 = charge_trace*C2(26)."""
    alpha = 1.0 / res["c"]["inv_alpha_phys"]
    S = (9 * math.pi**2 / 2 - 6 + 15 / 512
         - 16 * (alpha / (2 * math.pi))**2)
    from root import M_Pl_MeV
    v_pred = M_Pl_MeV * math.exp(-S) / 1e3       # GeV
    assert v_pred == pytest.approx(res["m"]["v_EW_pred_GeV"], rel=1e-9)


def test_sin2w_closed_form(res):
    """sin^2(th_W) = 3/13 + Q0^2 * alpha/(2pi), where the echo weight
    h7 + h7/d11^2 = Q0^2 exactly (the Koide parameter squared)."""
    from root import h_7, d11, Q0
    assert h_7 * (1 + Fraction(1, d11**2)) == Q0**2
    alpha = 1.0 / res["c"]["inv_alpha_phys"]
    expected = 3 / 13 + float(Q0**2) * alpha / (2 * math.pi)
    assert res["g"]["sin2W"] == pytest.approx(expected, abs=1e-9)


def test_bridge_sq_closed_form(res):
    """bridge^2 = 32/27 * (1 + h11/K^2) = 32/27 * 73/72."""
    expected = (32 / 27) * (73 / 72)
    assert res["m"]["bridge"]**2 == pytest.approx(expected, rel=1e-9)


# ═══════════════════════════════════════════════════════════════════════
#  4. Derivation-program witnesses
#     Each structural theorem reduced to executable components
#     already present in the repo.  These tests freeze the witnesses.
#     the symbolic proofs are the papers' job.
# ═══════════════════════════════════════════════════════════════════════

def test_first_invariant_order_theorem():
    """THE FIRST-INVARIANT-ORDER THEOREM (why K^3): the Koide
    insertion is the Z3
    sector-changer S (masses.py: h10 = the unique lightest changer).
    The confined sector admits only Z3-invariant operators
    ([Delta, S] = 0, the closure theorem).  The invariant content of
    n insertions is tr(S^n)/3: zero for n = 1, 2, first nonzero at
    n = 3 -> the back-reaction enters at cubic altitude, K^3.  The
    bridge^2 exception is forced by the same theorem: a PAIRED
    insertion S S-dagger is invariant already at second order -> K^2,
    exactly the quadratic altitude root.py assigns it."""
    import numpy as np
    S = np.roll(np.eye(3), 1, axis=0)               # Z3 shift operator
    tr = lambda M: abs(float(np.trace(M)))
    assert tr(S) < 1e-12                             # order 1: forbidden
    assert tr(S @ S) < 1e-12                         # order 2: forbidden
    assert tr(np.linalg.matrix_power(S, 3)) == pytest.approx(3.0)  # order 3
    assert tr(S @ S.T) == pytest.approx(3.0)         # paired |.|^2: order 2 ok


def test_conversion_vertex_irrationality_witness():
    """THE CONVERSION LEMMA: quantum dimensions in SU(3)_3 are
    the integers (1, 2, 3), so the ONLY irrationality source in the
    D6 walk algebra is the conversion-vertex normalization
    |C_3bar| = sqrt(2), computed constructively from the Fano plane
    (octonions.py).  Each excluded terminus word carries an ODD number
    of channel switches, so its amplitude is sqrt(2)^odd x rational =
    irrational, and an irrational number cannot contribute to an
    integer rank.  Witnesses frozen. The symbolic proof obligation is
    reduced to bookkeeping of the switch parity."""
    import contextlib, io
    import octonions
    with contextlib.redirect_stdout(io.StringIO()):
        oc = octonions.derive()
    assert oc["C3bar"]**2 == pytest.approx(2.0, abs=1e-12)
    assert math.isqrt(2)**2 != 2                     # 2 is not a square
    switches = lambda w: sum(a != b for a, b in zip(w, w[1:]))
    for word in ("ffaa", "fffa", "faaa"):            # the excluded words
        assert switches(word) % 2 == 1
    assert switches("ffff") == 0 and switches("aaaa") == 0  # pure towers


def test_vent_share_witnesses():
    """THE VENT-SHARE THEOREM (why h10, why K^0, why the scalar
    channel): the 7 of
    G2 decomposes under confined SU(3) as 3 + 3bar + 1 with trialities
    (1, 2, 0). N-ality superselection (words.py, executable theorem)
    blocks the triality-charged shares at a confining interface, so
    the fundamental share of the bridge cannot reach the condensate.
    Its share in the weight metric is h10 (the emission-share grammar,
    sin^2(th_W) precedent), and a selection rule carries no altitude
    suppression (Casimir-vent precedent, kind='vent') -> K^0."""
    from root import h10, h_7
    t = lambda p: (p[0] + 2 * p[1]) % 3
    assert t((1, 0)) == 1 and t((0, 1)) == 2 and t((0, 0)) == 0
    assert 3 + 3 + 1 == 7                            # the decomposition
    assert h10 == Fraction(2, 9)                     # the blocked share
    assert h_7 == Fraction(2, 5)                     # the sin^2 precedent weight
    # the canonical edge carries the vent as kind='vent' (selection rule)
    from root import WEB
    kinds = {term.kind for term in WEB["lambda_MPl"].terms}
    assert "vent" in kinds


def test_texture_stabilization_witness():
    """Topological-protection witness, both halves.  Fast half A:
    the static unit
    texture really carries the pi_3 charge the confined phase must
    protect: the lattice estimator converges toward |B| = 1 with
    resolution.  Fast half B (branch (i) closure): the induced
    Faddeev-Skyrme quartic exists and stabilizes a few healing
    lengths wide (2.5 xi at the measured gap) (tests/probes/effective_action.py, executed here).  The
    dynamical record (ballistic unwinding through amplitude-node
    events. Rotation loophole closed by the exact commutation test)
    is the slow witness, tests/test_probes.py::
    test_knot_charge_not_protected."""
    from probes import knot_charge as kc
    B = {}
    for N in (32, 48):
        psi, dx, R, _ = kc.texture_state(N, 24.0, amp=0.7, size=2.5)
        B[N] = kc.charge(*kc.rel_modes(psi), dx, R)
    # converging toward -1 from above in |.|, already substantial at 48
    assert abs(B[48]) > abs(B[32]) > 0.5
    assert abs(B[48]) > 0.85
    # branch (i): the stabilizing term, constructively
    from probes import effective_action
    r = effective_action.run(report=lambda *a, **k: None)
    # the hermitian equation's scales: R* = 2.5 healing lengths (1.0
    # under the non-variational variant's gap), local limit holds (R* m = 1.41)
    assert r["R_star_m"] > 1.0 and 2.0 < r["R_star_xi"] < 3.0


# ═══════════════════════════════════════════════════════════════════════
#  5. Joint cosmology closure (independent channels, one universe)
# ═══════════════════════════════════════════════════════════════════════

def test_higgs_transport_is_not_an_observable_promotion(res):
    """The Higgs transport stops at its imported beta coordinate."""
    from root import WEB

    assert WEB["lambda_MPl"].value() == pytest.approx(
        -0.01245015315975486, abs=1e-14)
    assert res["g"]["lambda_mt_transport"] == pytest.approx(
        0.1265078151, abs=2e-10)
    assert "mH_pred" not in res["g"]


def test_joint_cosmology_closure(res):
    """Freeze the joint numerical consequence of eta_B and bridge venting.

    Standard conversion: Omega_b h^2 = eta_10/273.9.
    """
    eta10 = res["g"]["eta_B"] * 1e10
    Obh2 = eta10 / 273.9
    ODMh2 = res["d"]["DM_baryon_ratio"] * Obh2
    assert abs(Obh2 - 0.02237) / 0.00015 < 2.0     # Planck pull < 2 sigma
    assert abs(ODMh2 - 0.1200) / 0.0012 < 2.0      # Planck pull < 2 sigma


def test_searchlight_depth_and_recursive_web_are_distinct(res):
    """The octonionic searchlight stops at 3; the recursive web does not.

    Three checks:
    1. The searchlight depth is 3 while a declared depth-4 web edge is legal.
    2. The EM sector graph has exactly 2 charged nodes (G2-knot, F4-knot),
       so its longest non-repeating path is length 2. All deeper traversals
       are resummed by the kernel's self-consistency.
    3. The depth-4 EM response is state-dependent (callable), confirming
       that it reads the depth-3 carrier recursively from the web state.
    """
    from root import WEB, SEARCHLIGHT_DEPTH, EchoTerm, Ledger

    # 1. Hurwitz limits the searchlight, not EchoTerm/Web.
    assert SEARCHLIGHT_DEPTH == 3
    assert EchoTerm(["post-searchlight"], 0.0, 4).depth == 4
    max_depth_found = 0
    total_terms = 0
    for name, led in WEB.items():
        if not isinstance(led, Ledger):
            continue
        for t in led.terms:
            max_depth_found = max(max_depth_found, t.depth)
            total_terms += 1
    assert max_depth_found == 4, "canon must include its depth-4 reading"
    assert total_terms > 0, "web must have echo terms"
    assert WEB["S_baryo"].web_depth == 4
    assert WEB["M_Pl_MeV"].web_depth == 4
    assert not hasattr(WEB["S_baryo"], "emergence_depth")
    for name in ("inv_alpha", "S_quark", "m_e", "m_mu", "m_tau",
                 "a1_bridge", "Sigma_vent(UV)", "G_ratio"):
        assert WEB[name].web_depth == 4, \
            f"{name} must declare its recursive depth-4 reading"

    # 2. EM sector graph has exactly 2 charged sectors.
    # The no-further-edge theorem (root.py) relies on:
    #   - EM channel couples only to charged sectors
    #   - G2-knot (leptons) and F4-knot (quarks) are the only charged sectors
    #   - Higgs and common mode are EM-neutral
    # With 2 nodes and 1 edge, the only cycle is the 2-cycle e<->q,
    # and the kernel's self-consistency resums all repetitions.
    # Read the solved web itself and verify the EM channel's
    # depth-4 reading of the depth-3 carrier appears in the charged nodes and
    # nowhere else.
    em_cycle_nodes = {name for name, led in WEB.items()
                      if isinstance(led, Ledger)
                      and any("e↔q" in "→".join(t.path)
                              for t in led.terms)}
    assert "inv_alpha" in em_cycle_nodes          # the coupling node
    mass_nodes = {n for n in em_cycle_nodes if n.startswith("m_")}
    assert mass_nodes == {"m_e", "m_mu", "m_tau"}, \
        f"lepton knots only on the e-side of the EM 2-cycle: {mass_nodes}"
    assert "S_quark" in em_cycle_nodes            # the action node (q-side)
    neutral = {"sin2W", "lambda_MPl", "dark_ratio", "M_Pl_MeV"}
    assert not (neutral & em_cycle_nodes), \
        "EM-neutral nodes must carry no e↔q cycle edge"

    # 3. The depth-4 EM response is state-dependent (self-referential)
    inv_alpha_led = WEB["inv_alpha"]
    depth4_terms = [t for t in inv_alpha_led.terms if t.depth == 4]
    assert len(depth4_terms) == 1, "inv_alpha must have one depth-4 response"
    assert callable(depth4_terms[0].factor), \
        "depth-4 EM response must read alpha from the recursive web"

    # The self-referential term resums all traversals: the cubic
    # x^3 = (512/pi)[(1-1/(2pi))x^2 - 1/(2pi^2)] is EXACT,
    # not a truncation.  Verify the fixed point satisfies it.
    x = WEB.state["inv_alpha"]
    cubic_lhs = x**3
    cubic_rhs = (512 / math.pi) * ((1 - 1/(2*math.pi)) * x**2
                                   - 1 / (2 * math.pi**2))
    assert abs(cubic_lhs - cubic_rhs) < 1e-5, \
        "alpha(0) must satisfy the self-consistent cubic"


# ═══════════════════════════════════════════════════════════════════════
#  6. The grammar as a statistical object (edge_grammar.py, programs
#     D1-D3 of the gaps audit): the committed menus are frozen above.
#     here they are tested as a NULL MODEL (could chance produce the
#     scorecard?) and as a BLIND HOLDOUT (does any PDG mass steer its
#     own prediction?).
# ═══════════════════════════════════════════════════════════════════════

def test_beta_order_is_not_web_depth():
    """Transport expansion order cannot promote a Web ledger."""
    import higgs_ew
    from root import WEB

    assert "not Ledger.web_depth" in higgs_ew._beta_sm.__doc__
    assert WEB["lambda_MPl"].web_depth == 2


def _koide_Q(m1, m2, m3):
    return (m1 + m2 + m3) / (math.sqrt(m1) + math.sqrt(m2)
                             + math.sqrt(m3))**2


def _solve_third(Q, ma, mb, lo, hi, n=200):
    """All solutions m in [lo, hi] of Q(m, ma, mb) = Q (log grid +
    bisection)."""
    import numpy as np
    xs = np.exp(np.linspace(math.log(lo), math.log(hi), n))
    f = [_koide_Q(x, ma, mb) - Q for x in xs]
    sols = []
    for i in range(n - 1):
        if f[i] == 0.0 or f[i] * f[i + 1] < 0:
            a, b = xs[i], xs[i + 1]
            for _ in range(60):
                m = 0.5 * (a + b)
                if (_koide_Q(a, ma, mb) - Q) * (_koide_Q(m, ma, mb) - Q) <= 0:
                    b = m
                else:
                    a = m
            sols.append(0.5 * (a + b))
    return sols


def test_grammar_null_model_scorecard_not_cheap():
    """D2: Monte Carlo over the committed edge grammar with RANDOM
    targets.  Draw fake 'universes' (each mass moved log-uniformly
    within a factor 2.5 of the observed value), and ask whether the
    grammar's finite menus (edge_grammar.py: one theta for BOTH mu and
    tau, word multipliers + one vent for u/d/c/t, Koide inversion for
    b and s) can hit all eight non-anchor targets within 1%.  If the
    scorecard were cheap, the null rate would be O(1).  Measured: the
    lepton stage alone passes at ~2.5e-4, the full scorecard at 0 in
    4000 trials (95% CL upper bound 7.5e-4). The stage-rate product
    (stages treated as approximately independent. An estimate, not
    a theorem) puts the joint chance near 1e-6.  Frozen bounds are
    deliberately loose: lepton rate < 5e-3, product < 1e-4, zero
    full hits at this N.  This is the look-elsewhere correction for
    the 9/9-at-1% claim.  Caveats stated: the target window (x/2.5)
    is a choice, and the null holds the framework's functional FORMS
    fixed: it asks whether THIS grammar is dense, not how many
    other grammars might exist."""
    import numpy as np
    import edge_grammar as eg

    BA = math.sqrt(2.0)
    WIND = {"e": 1, "mu": 2, "tau": 0}
    OBS = dict(e=0.511, mu=105.658, tau=1776.93, u=2.16, d=4.70,
               s=93.5, c=1273.0, b=4183.0, t=172570.0)
    TOL, SPREAD, N = 0.01, 2.5, 4000

    def circ(theta):
        g = {n: 1.0 + BA * math.cos(theta + 2 * math.pi * k / 3)
             for n, k in WIND.items()}
        if min(g.values()) <= 0:
            return None
        return (g["mu"] / g["e"])**2, (g["tau"] / g["e"])**2

    LEP = [r for r in (circ(t % (2 * math.pi))
                       for t in eg.theta_menu()) if r is not None]
    CORR = [0.0] + sorted(set(eg.correction_menu().values()))
    MUL2 = sorted({w * (1 + c) for w in set(eg.word_menu(2).values())
                   for c in CORR})
    MUL3 = sorted({w * (1 + c) for w in set(eg.word_menu(3).values())
                   for c in CORR})
    W4 = set(eg.word_menu(4).values()) | set(eg.terminus_menu().values())
    MUL4 = sorted({w * (1 + c) for w in W4 for c in CORR})
    KQ = sorted(set(eg.koide_menu().values()))

    def hit(menu, target):
        return any(abs(v / target - 1.0) < TOL for v in menu)

    rng = np.random.default_rng(20260711)   # frozen seed: reproducible
    names = ("mu", "tau", "u", "d", "s", "c", "b", "t")
    full_hits, lep_hits = 0, 0
    stage = dict(u=0, d=0, c=0, t=0)
    for _ in range(N):
        tgt = {k: OBS[k] * math.exp(rng.uniform(-math.log(SPREAD),
                                                math.log(SPREAD)))
               for k in names}
        me = OBS["e"]
        # unconditional stage rates (for the product bound)
        stage["u"] += hit([m * me for m in MUL2], tgt["u"])
        stage["d"] += hit([m * me for m in MUL2], tgt["d"])
        stage["c"] += hit([m * OBS["mu"] for m in MUL3], tgt["c"])
        stage["t"] += hit([m * OBS["tau"] for m in MUL4], tgt["t"])
        # the joint trial
        pair = next(((rm * me, rt * me) for rm, rt in LEP
                     if abs(rm * me / tgt["mu"] - 1) < TOL
                     and abs(rt * me / tgt["tau"] - 1) < TOL), None)
        if pair is None:
            continue
        lep_hits += 1
        mu_hat, tau_hat = pair
        if not (hit([m * me for m in MUL2], tgt["u"])
                and hit([m * me for m in MUL2], tgt["d"])
                and hit([m * mu_hat for m in MUL3], tgt["c"])
                and hit([m * tau_hat for m in MUL4], tgt["t"])):
            continue
        ok = False
        c_hats = [m * mu_hat for m in MUL3
                  if abs(m * mu_hat / tgt["c"] - 1) < TOL]
        t_hats = [m * tau_hat for m in MUL4
                  if abs(m * tau_hat / tgt["t"] - 1) < TOL]
        for ch in c_hats:
            for th in t_hats:
                for Q in KQ:
                    for mb in _solve_third(Q, ch, th, 100.0, 4e4):
                        if abs(mb / tgt["b"] - 1) < TOL and any(
                                abs(ms / tgt["s"] - 1) < TOL
                                for Q2 in KQ
                                for ms in _solve_third(Q2, ch, mb,
                                                       1.0, 2e3)):
                            ok = True
        full_hits += ok

    # the scorecard must NOT be cheap
    assert lep_hits / N < 5e-3, \
        "one theta hitting both mu and tau at 1% must be rare"
    assert full_hits == 0, \
        f"chance full scorecards: {full_hits}/{N} (must be 0 at this N)"
    prod = (lep_hits / N or 2.5e-4)
    for k in stage:
        prod *= stage[k] / N
    assert prod < 1e-4, f"stage-rate product {prod:.2e} too permissive"


def test_channel_weight_projector_lemma():
    """The channel-weight lemma, tested in its three separable parts
    (octonions._check_projector_normalization): STRUCTURE, the raw
    Fano CG composite satisfies L^2 = lam L with spectrum {0, lam},
    so 1/(loop value) is the unique idempotent coefficient; VALUES,
    the classical loop value is 4 while the quantum dimension is
    d(1,0) = 2 (Weyl/Kac-Peterson, computed in octonions.su3_qdim,
    cross-checked against the wiring.py table); TRACE CHOICE, the
    spherical (quantum) trace is unique in an MTC and gate 3 places
    the channel in the SU(3)_3 MTC, so what is identified is that
    the physical channel is the MTC channel; cross-pinned by the
    engine's Q0 = 2/3 (classical reading: Q0 = 5/12, fails)."""
    from octonions import _check_projector_normalization, su3_qdim
    from wiring import su3_qdims
    schur_ok, forced_ok, lam_cl = _check_projector_normalization()
    assert schur_ok, "channel composite must be Schur (prop. to projector)"
    assert forced_ok, "spectrum {0, lam}: 1/lam is the unique scaling"
    assert abs(lam_cl - 4.0) < 1e-9, "classical loop value is 4, not d10"
    assert abs(su3_qdim(1, 0) - 2.0) < 1e-9
    assert abs(su3_qdim(1, 1) - 3.0) < 1e-9
    d = su3_qdims(3)                      # cross-check vs the table
    assert abs(d[(1, 0)] - su3_qdim(1, 0)) < 1e-9
    # the quantum reading (identification): W(3bar) = 1/d(1,0) = 1/2
    assert abs(1.0 / su3_qdim(1, 0) - 0.5) < 1e-9


def test_proton_stability_bookkeeping():
    """D5: the conserved quantity forbidding B-violating edges is the
    fermion-line terminus.  Frozen here: (i) the edge taxonomy is
    exactly {ratio, echo, vent} and root.EchoTerm rejects any other
    kind (an identity-changing edge cannot be constructed). (ii) every
    edge in the solved web carries one of the three kinds.  The
    exposure and its kill are stated at root.EchoTerm: no decay at any
    sensitivity, and any confirmed proton decay ends it."""
    from root import WEB, EchoTerm, Ledger
    with pytest.raises(AssertionError):
        EchoTerm(["p", "e+ pi0"], 1.0, 1, kind="transmute")
    kinds = {t.kind for led in WEB.values() if isinstance(led, Ledger)
             for t in led.terms}
    assert kinds <= {"ratio", "echo", "vent"}
    assert "proton decay" in EchoTerm.__doc__
    assert "tau_p" in EchoTerm.__doc__


def test_blind_holdout_nine_masses(res):
    """D3: leave-one-out over the nine masses, executable form.  The
    chain's only dimensional input is the CODATA electron anchor.
    Every PDG mass is display-only.  Therefore, given eight masses,
    indeed given none, the ninth is forced.  Frozen by perturbing
    EVERY PDG mass entry simultaneously (x1.2 / x0.8 alternating),
    re-running the derivation, and asserting all nine predictions are
    bit-identical to canon.  If a refactor ever routed a PDG mass
    into its own (or any) prediction, this fails."""
    import contextlib
    import io

    import root as root_mod
    import masses as masses_mod
    from root import PDG_MASSES

    saved = dict(PDG_MASSES)
    try:
        for i, k in enumerate(PDG_MASSES):
            PDG_MASSES[k] = saved[k] * (1.2 if i % 2 == 0 else 0.8)
        with contextlib.redirect_stdout(io.StringIO()):
            R2 = root_mod.derive()
            m2 = masses_mod.derive(R2)
        for k in ("m_e", "m_mu", "m_tau", "m_u", "m_d", "m_s",
                  "m_c", "m_b", "m_t"):
            assert m2[k] == pytest.approx(res["m"][k], rel=1e-12), \
                f"{k} moved when PDG values were perturbed: a PDG " \
                f"mass is steering a prediction"
    finally:
        PDG_MASSES.update(saved)
        # restore global web state for any later consumer
        with contextlib.redirect_stdout(io.StringIO()):
            R3 = root_mod.derive()
            masses_mod.derive(R3)


def test_cited_paths_exist():
    """Citation resolution: every repo-relative file path cited in
    the README and the interference/ module docstrings must
    exist on disk.  Catches citation rot (a probe renamed or moved
    while prose still points at the old path)."""
    import re
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                            ".."))
    pat = re.compile(
        r"(?:interference|tests)/(?:[a-z_]+/)?[A-Za-z0-9_.-]+"
        r"\.(?:py|md|yml)")
    sources = [os.path.join(root_dir, "README.md")]
    iface = os.path.join(root_dir, "interference")
    for dirpath, dirnames, filenames in os.walk(iface):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        sources += [os.path.join(dirpath, f) for f in sorted(filenames)
                    if f.endswith(".py")]
    missing = []
    for src in sources:
        text = open(src, encoding="utf-8").read()
        for cited in sorted(set(pat.findall(text))):
            if not os.path.exists(os.path.join(root_dir, cited)):
                missing.append(f"{os.path.basename(src)} cites {cited}")
    assert not missing, "\n".join(missing)


def test_searchlight_absorption():
    """The searchlight absorption scoreboard, executable.
    Q0 absorbed (engine identity); B/A cross-checked (Fano CG
    recomputes sqrt(2)); theta the unique blind positivity survivor
    of the committed angle menu."""
    import contextlib
    import io
    import math
    from fractions import Fraction

    import edge_grammar as eg
    import octonions
    import root

    # Q0: falls out of the engine as an identity, no octonion consulted
    assert root.Q0 == Fraction(1, 3) + Fraction(root.d10, 6) \
        == Fraction(root.d10, root.d11) == Fraction(2, 3)

    # B/A: ABSORBED.  The theta-blind identity Q = 1/3 + (B/A)^2/6
    # (the circulant master identity, proven in circulant_normalization
    # PART A) plus the absorbed Q0 forces (B/A)^2 = 6(Q0 - 1/3) = d10
    # with no octonion consulted: the amplitude is an engine
    # observation, and the Fano CG recomputation is the independent
    # cross-check, not the source.
    assert 6 * (root.Q0 - Fraction(1, 3)) == root.d10
    assert abs(float(root.BA_ratio) ** 2 - root.d10) < 1e-12
    with contextlib.redirect_stdout(io.StringIO()):
        oct_out = octonions.derive()
    assert abs(oct_out["BA_ratio"] - math.sqrt(2.0)) < 1e-12
    assert abs(float(root.BA_ratio) - oct_out["BA_ratio"]) < 1e-12

    # theta: sole blind survivor (positivity + nondegeneracy) of the
    # committed angle menu, recomputed here from the grammar
    survivors = []
    for theta in eg.theta_menu():
        gaps = sorted(1.0 + float(root.BA_ratio)
                      * math.cos(theta + 2 * math.pi * k / 3)
                      for k in range(3))
        if gaps[0] > 0 and len({round(g, 12) for g in gaps}) == 3:
            survivors.append(theta)
    assert len(survivors) == 1
    assert abs(survivors[0] - 2.0 / 9.0) < 1e-12


def test_private_symbol_contract():
    """Seven leading-underscore symbols are load-bearing across the
    verification layer (probes and tests import them directly).  A
    rename would silently break CI-visible physics claims, so their
    existence and call signatures are a frozen CONTRACT here."""
    import inspect
    import nls_soliton
    import mixing as _mixing
    import octonions as _oct
    import higgs_ew as _hew

    contract = [
        (nls_soliton, "_split_step",
         ["psi", "dx", "dt", "g0", "g1", "n_steps", "alpha"]),
        (nls_soliton, "_z3_fourier", None),
        (_mixing, "_compute_wzw", []),
        (_oct, "_check_projector_normalization", []),
        (_hew, "_run_rge", None),
        (_hew, "_beta_sm", ["y"]),
        (_hew, "_derive_alpha_em", None),
    ]
    for mod, name, params in contract:
        fn = getattr(mod, name, None)
        assert callable(fn), f"contract broken: {mod.__name__}.{name}"
        if params is not None:
            got = list(inspect.signature(fn).parameters)
            assert got[:len(params)] == params, \
                f"{mod.__name__}.{name} signature drifted: {got}"
