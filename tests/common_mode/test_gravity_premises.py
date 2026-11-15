"""Freezes tests/probes/gravity_premises.py: the amount's share and the vector contact."""

import pytest


def _silent(*args, **kwargs):
    pass


@pytest.fixture(scope="module")
def r():
    from probes import gravity_premises
    return gravity_premises.run(report=_silent)


def test_the_share_is_tested_by_alpha_without_G(r):
    """1/alpha_alg = 48 b0/pi. One real 7 running puts 1/alpha(0) at
    134.89, and the labelled bridge puts it at 81.36.  The amount comes
    out at +0.13 sigma from the Cs recoil value."""
    assert r["b0"] == pytest.approx(32 / 3, abs=1e-12) and r["base_over_b0"] == pytest.approx(48.0, abs=1e-9)
    s = r["shares"]
    assert s["amount, b0 = 32/3"][2] == pytest.approx(0.13, abs=0.01)
    assert s["one real 7 running"][1] == pytest.approx(134.8948, abs=1e-3)
    assert s["labelled, b0 = 19/3"][1] == pytest.approx(81.3646, abs=1e-3)


def test_the_contact_is_the_weitzenbock_term(r):
    """On the unit S^4 the vector less two ghosts has a1 = -8 for the
    Hodge Laplacian, +4 for the rough one. Without the Weitzenbock term
    b0 = -16/3 and S0 < 0."""
    assert r["a1_hodge"] == pytest.approx(-8.0, abs=1e-4)
    assert r["a1_rough"] == pytest.approx(4.0, abs=1e-4)
    assert r["b0_no_spin"] == pytest.approx(-16 / 3, abs=1e-12) and r["S_no_spin"] < 0
