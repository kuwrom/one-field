"""Freezes tests/probes/induced_cubic.py: the induced action's cubic term is Einstein's."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_induced_cubic_term_is_einsteins_with_the_kernels_stiffness():
    """On backgrounds whose components vary by a factor of five, the
    scalar row's heat trace has c1 = (1/6) int sqrt g R to 1e-6 at both
    signs of the amplitude. The part odd in the amplitude, the cubic
    and higher odd terms, is 11% of c1 and agrees with Einstein's to
    1e-5. Xi multiplies the whole density. And a row whose x-derivative
    is scaled by c^2 = 1.3 induces the density of its own metric,
    c1/c0 larger by 1.3."""
    from probes import induced_cubic

    r = induced_cubic.run(report=_silent)
    assert r["metric_spread"] > 5.0 and r["expansion_error"] < 1e-6
    assert r["c1"][1.0] == pytest.approx(0.2668780, abs=1e-6) and r["c1"][-1.0] == pytest.approx(0.2078311, abs=1e-6)
    assert r["odd_share"] == pytest.approx(0.1106, abs=1e-3) and r["odd_error"] < 1e-5
    assert r["xi_error"] < 1e-6
    assert r["frame_ratio"] == pytest.approx(1.3, abs=1e-5)
