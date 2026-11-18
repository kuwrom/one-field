"""
test_g_path_trials.py: freeze the trials ledger of the G prediction.

What is frozen:
  1. each coarse menu on the G path (C_EFF, b0, loop order, alpha_G2)
     has every alternative off by an O(1) factor in M_Pl;
  2. the depth 2 coefficient carries M_Pl from -1140 to +22.2 ppm, on a
     menu spaced by alpha/2pi;
  3. G_F, through the M_Pl free ratio v_EW/m_e, admits exactly one of
     145 combinations (c2, w16) at N = 30 within 2 sigma, the canonical
     one, with c2 over the rival readings of frame_coefficients.py F2.
     The chance rate of such a match is computed for each scope of
     derived inputs;
  4. gate S, with no measured G, removes the reduced Planck boundary
     and bounds c.  The cutoff reading and the regulator at the
     boundary pass it, and the lane rule removes them.  Of the factors
     of general relativity, the gates as wired pass only k = 1 (G10).
     That exclusion rests on the wired threshold.
"""

import pytest

import g_path_trials


@pytest.fixture(scope="module")
def ledger():
    return g_path_trials.run(report=lambda *a: None)


def test_coarse_menus_are_order_one_apart(ledger):
    for name in ("C_EFF", "b0", "loop", "alpha_G2", "boundary"):
        alts = [abs(v) for k, v in ledger["menus"][name].items() if abs(v) > 100]
        assert len(alts) == len(ledger["menus"][name]) - 1
        assert min(alts) > 1e5


def test_the_boundary_is_the_identification_g_tests(ledger):
    b = ledger["menus"]["boundary"]
    assert b["full"] == pytest.approx(22.2, abs=0.1)
    assert (1 + 1e-6 * b["reduced"]) == pytest.approx(5.0132 * (1 + 22.2e-6), rel=1e-4)
    assert (1 + 1e-6 * b["cutoff"]) == pytest.approx(3**0.5 * (1 + 22.2e-6), rel=1e-4)
    assert ledger["menus"]["b0"]["21/2 (light G2 scalars)"] == pytest.approx(1.024e6, rel=2e-3)


def test_gate_s_bounds_the_boundary_and_g_selects(ledger):
    bg = ledger["boundary_gate"]
    assert bg["survives"] == {"full": True, "reduced": False, "cutoff": True, "regulator at M_b": True}
    assert [bg["c_max"][c] for c in (1.0, 0.0, -1.0)] == pytest.approx([4.112, 6.976, 9.841], abs=1e-3)
    assert max(bg["t0"]["reduced"].values()) == pytest.approx(0.392, abs=1e-3)
    assert min(bg["t0"]["cutoff"].values()) == pytest.approx(1.371, abs=1e-3)
    assert all(t == pytest.approx(1.0) for t in bg["t0"]["regulator at M_b"].values())
    assert bg["mpl_shift"]["regulator at M_b"] == pytest.approx(2.028, abs=1e-3)


def test_gr_menu_narrows_to_full_and_rationalized(ledger):
    bg = ledger["boundary_gate"]
    n = bg["n_conventions"]
    assert bg["gr_menu"]["1"] == n                      # passes under every convention
    assert bg["gr_menu"]["8pi"] == bg["gr_menu"]["16pi"] == bg["gr_menu"]["32pi"] == 0
    assert bg["gr_menu"]["4pi"] == 2                    # two of six conventions
    # the gates as wired (S at the nominal threshold, M's band) read no measured G
    assert bg["wired"]["1"] == (True, True)
    assert bg["wired"]["2pi"] == (False, True)          # fails S at the covariant contact
    assert bg["wired"]["4pi"] == (False, False)
    assert bg["wired"]["2 (off menu)"] == bg["wired"]["3 (off menu)"] == (True, True)
    assert bg["c_window"] == pytest.approx((0.656, 4.11), abs=2e-3)
    assert bg["gr_menu"]["2 (off menu)"] == bg["gr_menu"]["3 (off menu)"] == 4
    assert bg["x_half"] == pytest.approx(0.2674, abs=1e-4)


def test_depth2_coefficient_carries_the_ppm_step(ledger):
    assert ledger["depth1"] == pytest.approx(-1139.9, abs=0.5)
    assert ledger["G"]["G_ppm"] == pytest.approx(-54.6, abs=0.1)
    assert ledger["menus"]["c2"]["-1.0"] == pytest.approx(22.2, abs=0.1)
    assert ledger["spacing_ppm"] == pytest.approx(1161, abs=1)
    assert ledger["w4_span"] == pytest.approx(26.4, abs=0.2)


def test_gf_fixes_the_depth2_coefficient(ledger):
    gf = ledger["gf"]
    assert gf["n_combos"] == 145
    assert gf["hits"] == [((-1.0, 30, 16.0), 0.019)]
    assert gf["runner_up_ppm"] == pytest.approx(2.72, abs=0.02)
    assert gf["chance"]["N, w16 derived"] < 0.005
    assert 0.02 < gf["chance"]["N derived"] < 0.05
    assert set(gf["chance"]) == {"N, w16 derived", "N derived"}


def test_alpha_s_excludes_two_loop_g2_running(ledger):
    lo = ledger["loop"]
    assert lo["inv_alpha_star_one"] == pytest.approx(32 / 3.141592653589793, abs=1e-6)
    assert lo["alpha_MZ_one"] == pytest.approx(0.11839, abs=5e-5)
    assert (lo["alpha_MZ_two"] - 0.1180) / 0.0009 > 20


def test_alpha_s_needs_the_scheme_identification_at_mu_star(ledger):
    sc = ledger["scheme"]
    assert sc["k_needed"] == pytest.approx(0.8445, abs=1e-3)
    assert sc["Lambda_ratio"] == pytest.approx(2.327, abs=2e-3)
    assert (sc[0.0][1] - 0.1180) / 0.0009 < -15


def test_pi_over_32_is_the_msbar_coupling_to_one_percent(ledger):
    sb = ledger["band"]
    assert sb["band"][0] == pytest.approx(-0.0157, abs=5e-4)
    assert sb["band"][1] == pytest.approx(0.0101, abs=5e-4)
    assert (sb["alpha_s"]["V"] - 0.1180) / 0.0009 < -7
    assert (sb["alpha_s"]["MS"] - 0.1180) / 0.0009 > 30


def test_boundary_datum_has_one_scheme_free_reading(ledger):
    dr = ledger["datum_readings"]["mpl_shift"]
    assert dr["action"] == pytest.approx(1.0, abs=1e-9)
    assert dr["MS-bar"] == pytest.approx(0.2177, abs=1e-3)   # the G6 two loop alternative, 1/4.59
    assert dr["DR-bar"] == pytest.approx(0.2316, abs=1e-3)
    assert dr["V"] == pytest.approx(0.3030, abs=1e-3)        # analytic readings scatter with the scheme
