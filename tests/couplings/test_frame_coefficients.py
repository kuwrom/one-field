"""Freezes tests/probes/frame_coefficients.py: the ruler's three coefficients as objects the web holds."""

import pytest


def _silent(*args, **kwargs):
    pass


@pytest.fixture(scope="module")
def r():
    from probes import frame_coefficients
    return frame_coefficients.run(report=_silent)


def test_c_eff_is_the_closures_common_mode_read_whole(r):
    """P Delta^2 P = 2 and C_EFF = 1/2. The first order reading gives 1,
    which puts v_EW/m_e off by a factor of two. The three menu entries
    need Q = 3/4, 4/3, 3/2."""
    assert r["P_delta2"] == pytest.approx(2.0, abs=1e-12) and r["c_eff"] == pytest.approx(0.5, abs=1e-12)
    assert r["first_order"] == pytest.approx(1.0, abs=1e-12) and abs(r["gf_first"]) > 1e6
    assert r["menu_Q"] == pytest.approx({"4/9": 0.75, "1/4": 4 / 3, "2/9": 1.5}, abs=1e-12)
    assert r["menu_B2"]["4/9"] == pytest.approx(2.5, abs=1e-12)


def test_c2_is_the_bridge_vertex(r):
    """The response rule admits h_bridge alone, so c2 = -1. Every other
    weight, inserted by hand, fails G_F by more than 900 sigma."""
    assert r["admitted"] == ["h_bridge"] and r["c2"] == pytest.approx(-1.0, abs=1e-15)
    assert abs(r["c2_pulls"]["h_bridge"]) < 1
    assert min(abs(v) for k, v in r["c2_pulls"].items() if k != "h_bridge") > 900


def test_w4_is_the_share_of_the_unital_traversal(r):
    """w4 = 27/2. The traversal without its unit fails m_mu/m_e at 16
    sigma, the whole traversal at 440. The value as read puts G_F at
    +0.08 sigma and m_mu at +0.12 sigma."""
    p = r["w4_pulls"]
    assert r["w4"] == 13.5
    assert p["13 (no unit)"] == pytest.approx(16.5, abs=0.2)
    assert p["27 (the whole traversal)"] == pytest.approx(-441.8, abs=1)
    assert r["gf"] == pytest.approx(0.08, abs=0.01) and r["mu"] == pytest.approx(0.12, abs=0.01)
