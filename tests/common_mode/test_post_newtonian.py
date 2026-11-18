"""Freezes tests/probes/post_newtonian.py: gamma = beta = 1 and one G for static and radiated fields."""

import pytest

pytest.importorskip("sympy", reason="needs sympy (pip install -r requirements.txt)")


def _silent(*args, **kwargs):
    pass


def test_post_newtonian_gamma_beta_and_one_stiffness():
    """The vacuum equations of the Einstein action fix gamma = 1 at
    first order and beta = delta = 1 at second. The kernel's
    F0/F2 = -2 is gamma = 1, and its q^2 term is local only at that
    ratio (composite_pole.py K9). The kernel's stiffness is approached
    with equal and opposite slopes from spacelike and timelike
    momenta."""
    from probes import post_newtonian as pn

    r = pn.run(report=_silent)
    assert r["gamma"] == 1 and r["beta"] == 1 and r["delta"] == 1
    assert r["spacelike"] < 0 < r["timelike"] and abs(r["spacelike"] + r["timelike"]) < 2e-8
