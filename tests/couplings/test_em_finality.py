"""Freezes tests/probes/em_finality.py: EM channel finality as a tested exclusion."""

import pytest


def _silent(*args, **kwargs):
    pass


@pytest.fixture(scope="module")
def r():
    from probes import em_finality
    return em_finality.run(report=_silent)


def test_the_sector_graph_is_bipartite(r):
    assert r["odd"] == [0, 0, 0] and r["odd_loop"] == [1, 4, 11]


def test_the_self_loop_fails_every_measurement(r):
    """n = 1 and 2 move 1/alpha by +5.1e-7 and +1.0e-6 and miss every
    value by 18 sigma or more. A double traversal is below 0.1 sigma."""
    t = r["table"]
    assert t["self-loop n = 1"][1] == pytest.approx(5.11e-7, rel=1e-2)
    assert min(abs(p) for p in t["self-loop n = 1"][2].values()) > 18
    assert min(abs(p) for p in t["self-loop n = 2"][2].values()) > 37
    assert abs(t["carrier twice"][1]) < 2.7e-9
    assert t["canon"][2]["Cs recoil"] == pytest.approx(0.13, abs=0.01)
    assert t["canon"][2]["Rb recoil"] == pytest.approx(-14.2, abs=0.1)


def test_what_a_rescue_needs(r):
    """The data require -0.01, +0.46, +0.50, +0.61. The self-loop supplies 2 and 4."""
    a = r["ask"]
    assert a["Cs recoil"] == pytest.approx(-0.014, abs=0.005)
    assert a["electron moment"] == pytest.approx(0.456, abs=0.005)
    assert a["CODATA 2022"] == pytest.approx(0.499, abs=0.005)
    assert a["Rb recoil"] == pytest.approx(0.613, abs=0.005)
    assert r["supplied"][1] == pytest.approx(2.0, abs=0.01) and r["supplied"][2] == pytest.approx(4.0, abs=0.02)


def test_the_one_landing_needs_a_matter_vertex_on_the_coupling(r):
    """Of seven conformal weights on the self-loop, only h10 = 2/9
    matches the data. Its pulls: electron moment -0.2, CODATA -0.7,
    Cs +4.3, Rb -3.9."""
    assert r["lands"] == ["2/9"]
    p = r["landing"]["2/9"]
    assert p["electron moment"] == pytest.approx(-0.2, abs=0.05)
    assert p["Cs recoil"] == pytest.approx(4.3, abs=0.1) and p["Rb recoil"] == pytest.approx(-3.9, abs=0.1)
