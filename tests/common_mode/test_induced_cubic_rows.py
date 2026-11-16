"""Freezes tests/probes/induced_cubic_rows.py: the Weyl and vector rows' cubic term is Einstein's."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_weyl_and_vector_rows_induce_the_einstein_density_beyond_quadratic_order():
    """On backgrounds that are not small the squared Dirac operator and
    the vector's Hodge operator, with their spin connection and
    curvature terms, have c1 = (4/6 - 1) int sqrt g R to 1e-6 at both
    signs of the amplitude: the Dirac field weighs 1/3, four states at
    1/12, and the gauge field less its ghosts -2/3, two states at
    1/6 - 1/2."""
    from probes import induced_cubic_rows as icr

    r = icr.run(report=_silent)
    for rep in ("spinor", "vector"):
        assert r[rep]["c1"] == pytest.approx(-0.5337559, abs=2e-6)
        assert r[rep]["error"] < 1e-6 and r[rep]["odd_error"] < 1e-5
    assert r["dirac_weight"] == pytest.approx(1 / 3) and r["gauge_weight"] == pytest.approx(-2 / 3)
