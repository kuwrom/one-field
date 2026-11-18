"""Scorecard tests of the mixing branch: every number the engine prints for it is asserted."""

import math
from fractions import Fraction

import numpy as np
import pytest

from root import PDG_MASSES


def err_pct(predicted: float, expected: float) -> float:
    """Absolute percentage deviation from the expected value."""
    return abs((predicted - expected) / expected) * 100.0


# ═══════════════════════════════════════════════════════════════════════
#  Quark masses  (6 predictions)
# ═══════════════════════════════════════════════════════════════════════



# ═══════════════════════════════════════════════════════════════════════
#  CKM sector  (Wolfenstein + UT angles, summarised)
# ═══════════════════════════════════════════════════════════════════════

def test_ckm_diagonal_pull_score(res):
    """CKM correlated coordinates have a small descriptive pull score."""
    n = len(res["mx"]["ckm_obs"])
    score_per_coordinate = res["mx"]["diagonal_pull_score_ckm"] / n
    assert score_per_coordinate < 1.0


def test_ckm_max_pull(res):
    """No individual CKM observable deviates by more than 2 sigma."""
    assert res["mx"]["max_pull_ckm"] < 2.0


# ═══════════════════════════════════════════════════════════════════════
#  PMNS sector (3 mixing angles + convention dependent U_e3 phase)
# ═══════════════════════════════════════════════════════════════════════

def test_pmns_diagonal_pull_score(res):
    """PMNS marginal pulls have a small descriptive diagonal score."""
    score_per_angle = res["mx"]["diagonal_pull_score_pmns"] / 3
    assert score_per_angle < 1.0


def test_pmns_generation_block_is_structural(res):
    """The projected bulk fusion blocks have the stated TBM pattern."""
    assert res["mx"]["generation_orbits"] == (
        ((0, 1), (1, 0)), ((0, 2), (2, 0)), ((1, 2), (2, 1)))
    assert res["mx"]["democratic_eigenvalues"] == pytest.approx(
        (3.0, 0.0, 0.0), abs=1e-12)
    assert res["mx"]["fundamental_eigenvalues"] == pytest.approx(
        (-0.5, 1.0, 0.5), abs=1e-12)
    tbm_abs = np.array([
        [2/math.sqrt(6), 1/math.sqrt(3), 0],
        [1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)],
        [1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)],
    ])
    assert np.allclose(np.abs(res["mx"]["U_TBM"]), tbm_abs)


def test_pmns_fold_is_not_mislabeled_as_the_zc_nimrep(res):
    """Z_C has two boundaries. The three orbit projection is not closed."""
    assert res["mx"]["zc_boundary_states"] == 2
    assert np.array_equal(res["mx"]["zc_n10"], np.ones((2, 2), dtype=int))
    assert np.array_equal(
        res["mx"]["zc_n11"], np.array([[1, 2], [2, 1]], dtype=int))
    assert res["mx"]["zc_fusion_residual"] == 0.0
    assert res["mx"]["fundamental_generation_leakage"] > 1e-12
    assert res["mx"]["fold_fusion_residual"] > 1e-12


def test_ckm_d6_boundary_data_are_computed(res):
    """The pair fraction and braid magnitude come from the graph/cells."""
    assert res["mx"]["d6_boundary_pairs"] == ((0, 5), (1, 2), (3, 4))
    expected = np.array([[1.0, 2.0], [2.0, 1.0]])
    assert all(np.array_equal(block, expected)
               for block in res["mx"]["d6_adjoint_blocks"])
    assert res["mx"]["A"]**2 == pytest.approx(2/3, abs=1e-15)
    cell_hi, cell_lo = res["mx"]["d6_cell_magnitudes"]
    assert cell_hi**2 == pytest.approx((math.sqrt(3) + 1)/2, abs=1e-15)
    assert cell_lo**2 == pytest.approx((math.sqrt(3) - 1)/2, abs=1e-15)
    assert res["mx"]["d6_connection_magnitude"] == pytest.approx(
        1/math.sqrt(2), abs=1e-15)
    assert res["mx"]["d6_braid_residual"] < 1e-12
    assert sorted(res["mx"]["d6_depth4_carrier_dimensions"]) == [2, 3, 3]


def test_pmns_delta_invariant(res):
    """The invariant PMNS delta sits inside NuFIT 6.0 NO 3sigma (96-422)."""
    delta = res["mx"]["delta_invariant_deg"]
    assert 96 < delta < 422, f"delta = {delta:.2f} outside NuFIT 3sigma"


def test_ckm_lambda_is_the_base_plus_echo_form():
    """λ exactly reuses bridge²'s 73/72 factor: λ = h10 · 73/72."""
    import root

    assert 1 + root.h11 / root.K**2 == Fraction(73, 72)
    assert root.h10 * (1 + root.h11 / root.K**2) == Fraction(73, 324)
    assert abs(float(Fraction(73, 324)) - 0.2253086) < 1e-6


def test_lambda_echo_menu_unique_hit(res):
    """The committed backreaction h10(1 + h11/K^2) is the unique
    1 sigma hit over the declared mixed numerical menu
    {h10, h11, Q0, 1} x {K, K^2, K^3}.  The enumeration runs as code,
    and this test freezes the count."""
    menu = res["mx"]["lambda_menu"]
    assert menu["n_candidates"] == 12
    assert menu["n_within_1sigma"] == 1
    assert menu["committed"] == ("h11", 2)
    assert menu["committed_is_unique_hit"]
    assert menu["hit_keys"] == [("h11", 2)]
    assert abs(menu["committed_pull"]) < 0.5
    assert menu["next_best_abs_pull"] > 1.9


def test_modular_phase_candidate_for_eta_bar():
    """The computed modular phase is pi/9. Assigning it to eta_bar is IH."""
    import mixing
    import root

    wzw = mixing._compute_wzw()
    fund = wzw["IDX"][(1, 0)]
    phase = math.atan2(wzw["S_mat"][fund, fund].imag,
                       wzw["S_mat"][fund, fund].real)
    assert phase == pytest.approx(math.pi / root.d11**2, abs=1e-12)
