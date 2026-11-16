"""Freezes tests/probes/coupling_normalization.py: the boundary factor k from the coupling variable."""

import math

import pytest


@pytest.fixture(scope="module")
def r():
    from probes import coupling_normalization
    return coupling_normalization.run(report=lambda *a: None)


def test_the_variable_and_the_factor(r):
    """The measured alpha selects the force law variable by 4 pi either
    way. In it the common mode gives k = 1, G at -2.43 sigma."""
    p = r["pulls"]
    assert abs(p["force law alpha"]) < 1
    assert abs(p["rationalized e^2 = 4 pi alpha"]) > 1e9 and abs(p["1/(4 pi alpha)"]) > 1e10
    assert r["g_pulls"]["1"] == pytest.approx(-2.43, abs=0.01)
    assert r["G"]["1"] / r["G"]["4pi"] == pytest.approx(4 * math.pi, rel=1e-12)
