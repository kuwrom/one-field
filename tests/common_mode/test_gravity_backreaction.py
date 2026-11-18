"""Freezes tests/probes/gravity_backreaction.py: gravity's backreaction on the masses, on the running and on itself."""

import math

import pytest


def _silent(*args, **kwargs):
    pass


def test_no_road_back_and_the_running_holds_no_gravitational_term():
    """No web node reads a gravity node. The exchange coefficients run
    from 1.75e-45 (electron) to 2.0e-34 (top), 26 orders below the
    finest mass comparison. A gravitational term in the running at
    -3/pi moves S0 to 17.48 and M_Pl by twelve orders, at -1/pi to
    32.45 and five orders; 1/alpha(0) holds |a0| below 5e-11 with no
    measured G and G below 1e-6. The metric's own loop takes the
    stiffness from 25.83 to 23.17 or 18.17 and the gates pass. A
    self-echo of the common mode is held below 6e-4 and 3e-3."""
    from probes import gravity_backreaction

    r = gravity_backreaction.run(report=_silent)
    assert r["readers_of_gravity"]["G_ratio"] == [] and r["readers_of_gravity"]["Sigma_vent(UV)"] == ["G_ratio"]
    assert r["readers_of_ruler"] == ["m_e", "m_mu", "m_tau", "v_EW_GeV"]
    assert r["largest"] == "m_t" and r["coefficients"]["m_t"] == pytest.approx(2.0e-34, rel=0.02)
    assert r["orders_below"] == pytest.approx(26.0, abs=0.2)
    assert r["a0"]["Robinson and Wilczek"] == pytest.approx(-3 / math.pi)
    assert r["S_moved"]["Robinson and Wilczek"] == pytest.approx(17.48, abs=0.02)
    assert r["S_moved"]["Toms, proper time cutoff"] == pytest.approx(32.45, abs=0.02)
    assert r["orders"]["Robinson and Wilczek"] == pytest.approx(-11.7, abs=0.1)
    assert r["orders"]["Toms, proper time cutoff"] == pytest.approx(-5.2, abs=0.1)
    assert r["orders"]["dimensional regularization, momentum cutoff"] == 0.0
    assert r["dlnalpha_da0"] == pytest.approx(7.32, abs=0.05)
    assert r["a0_held_by_alpha"] == pytest.approx(5.5e-11, rel=0.1) and r["a0_held_by_G"] < 2e-6
    assert r["sigma_with_loop"][0.0] == pytest.approx(155 / 6)
    assert all(all(g) for g in r["gates_with_loop"].values())
    assert r["self_weight_bound"][1] < 1e-3 and r["self_weight_bound"][2] < 5e-3
