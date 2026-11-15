"""Freezes tests/probes/common_mode_count.py: the common mode as the web's count, and the one reading that survives."""

import pytest


def _silent(*args, **kwargs):
    pass


@pytest.fixture(scope="module")
def out():
    from probes import common_mode_count
    return common_mode_count.run(report=_silent)


def test_c1_a_common_factor_is_no_record(out):
    """The depth two factor common to the three lepton rates moves the
    readings in M_b units (K1 to K3) by -1161 ppm.  It does not move the
    ratio readings K4, K5, K7, K8.  So K1, K2, K3 fail gate U."""
    u = out["gate_U"]
    assert all(u[k] * 1e6 == pytest.approx(-1161.4, abs=0.1) for k in ("K1", "K2", "K3"))
    assert all(abs(u[k]) < 1e-14 for k in ("K4", "K5", "K7", "K8"))


def test_c2_the_nodes_own_object(out):
    """Base 0 and two edges, 1/2 each on the searchlight's circulant,
    summing on the solved web to P v^2 P: the node's object."""
    from interference import root as r

    node = r.WEB["common_mode"]
    assert node.base == 0.0 and len(node.terms) == 2
    assert [round(t.value(dict(r.WEB.state, eq_loop=0.0)), 15) for t in node.terms] == [0.5, 0.5]
    assert r.COMMON_MODE == pytest.approx(out["K5"], abs=1e-14)


def test_c4_the_quark_triplets_are_not_the_population(out):
    """The crossing triplets' moments 3 Q - 1 = 1 + 6.9e-3 and 1 + 3.1e-3
    put G over 400 standard deviations high."""
    from interference import root as r
    from probes import common_mode_count as cc

    st = r.WEB.state
    assert 3 * st["Q_cbt"] - 2 == pytest.approx(6.944e-3, abs=1e-5)
    assert 3 * st["Q_scb"] - 2 == pytest.approx(3.086e-3, abs=1e-5)
    k8 = cc.readings(st)["K8"]
    assert cc.pulls(r.G_DEPTH_ONE * k8)["recommended"] > 400


def test_c5_the_survivor_and_the_measured_masses(out):
    """G M_b^2 = 3 Q_l - 1 = 1 - 7.50 L with Q_l 5.06 ppm below 2/3. The
    measured lepton masses give 1 + (-4.9 +- 11.3) L, and the engine
    lies -0.23 of the spread from it."""
    from interference import root as r

    L = r.WEB.state["eq_loop"]
    assert 3 * out["Q"] - 1 == pytest.approx(out["K5"], abs=1e-14)
    assert (out["K5"] - 1) / L == pytest.approx(-7.502, abs=1e-3)
    assert (out["Q"] / (2 / 3) - 1) * 1e6 == pytest.approx(-5.06, abs=0.01)
    assert (out["measured"] - 1) / L == pytest.approx(-4.9, abs=0.05)
    assert out["spread"] / L == pytest.approx(11.3, abs=0.05)
    assert out["pull_measured"] == pytest.approx(-0.23, abs=0.01)


def test_c6_the_test(out):
    """G = 6.673936e-11: -2.43 sigma of the recommended value, +0.17 with
    a spread added, -0.63 against the median."""
    from interference import root as r

    assert out["G"] == pytest.approx(r.G_PRED, rel=1e-14)
    assert out["G"] == pytest.approx(6.673936e-11, abs=2e-17)
    p = {k: round(v, 2) for k, v in out["pulls"].items()}
    assert p == {"recommended": -2.43, "spread added": 0.17, "median": -0.63}
