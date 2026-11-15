"""Scorecard tests of the foundations branch: every number the engine prints for it is asserted."""

import math
from fractions import Fraction

import numpy as np
import pytest

from root import PDG_MASSES


def err_pct(predicted: float, expected: float) -> float:
    """Absolute percentage deviation from the expected value."""
    return abs((predicted - expected) / expected) * 100.0


# ═══════════════════════════════════════════════════════════════════════
#  Root constants  (frozen by value)
# ═══════════════════════════════════════════════════════════════════════


_ROOT_CONSTANTS = [
    ("hv_G2", 4), ("hv_F4", 9),         # dual Coxeter numbers
    ("C2_7", 2), ("C2_26", 6),          # Casimirs
    ("K", 6),                           # WZW altitude
    ("h10", 2 / 9), ("h11", 1 / 2),     # conformal weights
    ("sin2W", 3 / 13),                  # zero edge Weinberg base
    ("N_bridge", 182),                  # bridge mode count
]


# ═══════════════════════════════════════════════════════════════════════
#  Structural closure  (E_8 mode accounting)
# ═══════════════════════════════════════════════════════════════════════

def test_e8_mode_assignment_closure(res):
    """248 = 14 + 52 + 182.  Every E_8 mode has a physical role."""
    R = res["R"]
    assert R["dim_G2"] + R["dim_F4"] + R["N_bridge"] == 248
    assert R["dim_G2"] == 14
    assert R["dim_F4"] == 52
    assert R["N_bridge"] == 182


# ═══════════════════════════════════════════════════════════════════════
#  Polynomial identities  (the four integer reductions)
#  ──────────────────────────────────────────────────────────────────────
#  These check "does the named quantity equal the claimed polynomial in
#  (d_10, d_11, n_7, n_26)?". They carry the unified version's main
#  claim: 14+ named quantities reduce to these four integers.
# ═══════════════════════════════════════════════════════════════════════

def test_four_irreducible_integers(res):
    """(d_10, d_11, n_7, n_26) = (2, 3, 7, 26)."""
    R = res["R"]
    assert (R["d10"], R["d11"], R["n7"], R["n26"]) == (2, 3, 7, 26)


def test_sugawara_consistency_d11sq_minus_1_eq_4_d10(res):
    """The integer identity d_11^2 - 1 = 4 d_10 forces (d_10, d_11) = (2, 3).

    Among positive integers, this Diophantine constraint together with
    d_10 < d_11 and d_10, d_11 >= 2 (both quantum dimensions of nontrivial
    primaries) has the unique solution (d_10, d_11) = (2, 3).
    Equivalently: d_11^2 = 4 d_10 + 1.  The enumeration below runs the
    uniqueness claim over 2 <= d_10 < d_11 <= 200.
    """
    N = 200
    solutions = [(a, b) for b in range(3, N + 1) for a in range(2, b)
                 if b * b - 1 == 4 * a]
    assert solutions == [(2, 3)]
    R = res["R"]
    assert (R["d10"], R["d11"]) == solutions[0]


def test_vertex_current_trace_is_30(res):
    """The E8 positive root current index supplies the vertex factor."""
    R = res["R"]
    assert R["vertex_charge_trace"] == R["N_vertex"] == R["hv_E8"] == 30


def test_root_constants_frozen(res):
    """The defined constants of the algebraic root, frozen by value."""
    R = res["R"]
    for key, expected in _ROOT_CONSTANTS:
        assert math.isclose(float(R[key]), expected, rel_tol=1e-12), key


def test_Q0_is_d10_over_d11(res):
    """Koide value Q_0 = d_10 / d_11 = 2/3."""
    R = res["R"]
    assert math.isclose(R["Q0"], R["d10"] / R["d11"], rel_tol=1e-12)
    assert math.isclose(R["Q0"], 2 / 3, rel_tol=1e-12)


def test_h_bridge_is_one(res):
    """h_bridge = h(7,G_2) + h(26,F_4) = 2/5 + 3/5 = 1 (dimension one)."""
    R = res["R"]
    assert math.isclose(R["h_7"], 2 / 5, rel_tol=1e-12)
    assert math.isclose(R["h_26"], 3 / 5, rel_tol=1e-12)
    assert math.isclose(R["h_bridge"], 1.0, rel_tol=1e-12)
    # The Fraction identity: h(R) = C_2(R) / (1 + h_dual)
    assert math.isclose(R["h_7"], R["C2_7"] / (1 + R["hv_G2"]), rel_tol=1e-12)
    assert math.isclose(R["h_26"], R["C2_26"] / (1 + R["hv_F4"]), rel_tol=1e-12)


def test_c_coset_vanishes(res):
    """c_coset = c(E_8) - c(G_2) - c(F_4) = 8 - 14/5 - 26/5 = 0."""
    R = res["R"]
    assert math.isclose(R["c_coset"], 0.0, abs_tol=1e-12)


def test_linearity_from_closure():
    """L_x linear + |xy|=|x||y| on O. Sedenion zero divisor. Z3 rule."""
    import numpy as np
    import octonions as oc
    rng = np.random.default_rng(11)
    assert oc._check_linearity(rng)
    assert oc._check_hurwitz_gate()
    assert oc._check_z3_selection()


def test_sedenion_zero_divisor_census():
    """The '168 pairs' remark is computed, not quoted: full enumeration
    of (e_i ± e_j)(e_k ± e_l) = 0 on the unit coefficient sedenion
    lattice gives 336 ordered = 168 unordered zero divisor pairs."""
    import octonions as oc
    assert oc._check_hurwitz_gate(census=True)


# ═══════════════════════════════════════════════════════════════════════
#  Embedding uniqueness (E₈ coherence gate)
# ═══════════════════════════════════════════════════════════════════════

def test_embedding_uniqueness():
    """E₈(1) ⊃ G₂(1) × F₄(1) is the unique pair among the rank-8
    candidates that passes all six coherence gates. The result is
    computed from the scan, not asserted as literals. The E₇ × A₁
    contender passes (a)-(e) and fails only (f), trivial centres. It is
    frozen alongside. Gates (e) and (f) are framework closure
    requirements, and the scan shows that they decide the result."""
    import embedding_uniqueness as eu
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        result = eu.derive()
    assert result['n_survivors'] == 1
    assert result['survivor'] == ('E8', 'G2', 'F4')
    nm = {(p, frozenset((g, f))): failed
          for p, g, f, failed in result['near_misses']}
    assert nm == {('E8', frozenset(('E7', 'A1'))): '(f) trivial centres'}


# ═══════════════════════════════════════════════════════════════════════
#  Protected forgetting (G₂ harmonic)
# ═══════════════════════════════════════════════════════════════════════

def test_protected_forgetting():
    """The canonical Z₃ harmonic (B, θ) = (√2, 2/9): PvP = 0, total
    second moment Pv²P = 1, per mode power E[v²] = 1/2, Koide Q = 2/3.

    The excluded harmonic (1, -1/2, -1/2) is (B, θ) = (1, 0). It gives
    Q = 1/2 with two masses degenerate. θ = 0 is degenerate for every
    amplitude. Both facts are computed, not asserted.
    """
    import protected_forgetting as pf
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        result = pf.derive()
    assert abs(result['PvP']) < 1e-14
    # total vs per mode: Pv²P = 2·E[v²], and the web consumes E[v²]
    assert abs(result['Pv2P_total'] - 1.0) < 1e-12
    assert abs(result['E_v2_mode'] - 0.5) < 1e-12
    assert abs(result['Pv2P_total'] - 2 * result['E_v2_mode']) < 1e-12
    # Koide value, and the amplitude it pins
    assert abs(result['Q'] - 2.0 / 3.0) < 1e-12
    assert result['amplitude_pinned_by_koide'] is True
    # the excluded harmonic is degenerate at Q = 1/2
    assert abs(result['retired_Q'] - 0.5) < 1e-12
    assert result['theta0_degenerate'] is True
    # the short root direction gives three distinct masses
    assert result['short_root_three_distinct'] is True
    # canonical harmonic stays on the positive gap branch
    assert result['gap_min'] > 0


def test_bell_fork():
    """The local realist reading (outcomes as assigned +/-1 facts) is
    capped at 2. The framework's local nonrealist reading (the
    ledger) sits exactly at the Tsirelson bound."""
    from probes import bell
    local = bell.chsh_local_bound(n_strategies=2000, n_lambda=200)
    assert local["achieved"] == 2.0
    assert local["best_random"] <= 2.0 + 1e-12
    assert abs(bell.chsh_ledger() - 2.0 * math.sqrt(2.0)) < 1e-12


# ═══════════════════════════════════════════════════════════════════════
#  The derivation chain, asserted against the engine
#
#  These are the identities the chain stands on: the dimensional
#  transmutation actions and their two face readings, the depth-4
#  response weights, the per mode condensate reading behind the
#  alpha_s threshold, and the CKM lambda form.  Each is computed in
#  the module that owns it. Each is read back here so a refactor
#  cannot quietly move one.
# ═══════════════════════════════════════════════════════════════════════

def test_transmutation_actions_and_their_two_faces():
    """S_lepton = 2pi/(b0 alpha) on the F4 face. The h-dual metric
    transports it to the G2 face: b_face = b0 h(F4)/h(G2) = 24, which is
    also the embedding invariant h(E8) - C2(26). Hence S_baryo =
    (4/9) S_lepton = 2 pi^2, read at web depth 4."""
    import root

    assert root.b0_G2 == Fraction(32, 3)
    assert root.S_lepton == pytest.approx(
        2.0 * math.pi / (float(root.b0_G2) * root.alpha_G2_Pl), abs=1e-12)
    assert root.S_lepton == pytest.approx(
        9.0 * math.pi**2 / 2.0, abs=1e-12)
    assert root.b_baryo_face == root.hv_E8 - root.C2_26 == 24
    assert root.S_baryo == pytest.approx(
        (root.hv_G2 / root.hv_F4) * root.S_lepton, abs=1e-12)
    assert root.S_baryo == pytest.approx(2 * math.pi**2, abs=1e-12)
    assert root.WEB["S_baryo"].web_depth == 4
    assert root.alpha_G2_Pl == pytest.approx(1.0 / (24.0 * math.pi),
                                             abs=1e-15)
    assert root.delta_S == pytest.approx(15.0 / 512.0, abs=1e-15)
    # the G₂ lever arm ends on π/32: 24π − (16/3π)(9π²/2 − 6) = 32/π
    inv_mu_star = (1.0 / root.alpha_G2_Pl
                   - float(root.b0_G2) / (2.0 * math.pi)
                   * (root.S_lepton - float(root.C2_26)))
    assert inv_mu_star == pytest.approx(32.0 / math.pi, abs=1e-9)
    assert root.alpha_G2_WZW == pytest.approx(math.pi / 32.0, abs=1e-15)
    assert root._inv0 == 1.0 / root.alpha_EM
