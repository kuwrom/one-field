"""Freezes tests/probes/growth_form.py: the linear growth the count
reading gives on the expansion's form."""

import pytest

pytest.importorskip("sympy", reason="needs sympy (pip install -r requirements.txt)")


def _silent(*args, **kwargs):
    pass


def test_growth_form_one_scale_free_factor_and_lambda_does_not_perturb():
    """A region of the count minus its background obeys
    delta'' + 2 H delta' = 4 pi G rho delta with Lambda absent; D = a
    in the matter era. Under the constant vacuum term the growth index
    is 0.55.  The amount's share 1 - 1/(2 pi) gives the held exponent
    0.9006, and the released baryons reach 0.996 of its excess today."""
    from probes import growth_form as gf

    r = gf.run(report=_silent)
    assert r["growth_equation"] and not r["lambda_perturbs"]
    assert r["matter_era_error"] < 1e-8
    assert list(r["growth_index"].values()) == pytest.approx([0.5562, 0.5543, 0.5535], abs=2e-4)
    assert r["held_exponent"] == pytest.approx(0.90055, abs=1e-5)
    assert r["held_numeric"] == pytest.approx(r["held_exponent"], abs=1e-6)
    assert r["baryons_at_tenth"] == pytest.approx(0.9705, abs=1e-3)
    assert r["baryons_today"] == pytest.approx(0.9961, abs=1e-3)


def test_growth_form_what_would_change_it():
    """A sound speed at the CMB bound gives a 73 Mpc Jeans length. On
    boundary_whole.py B4's menu the remainder alone (dark ratio
    infinite, share 1) gives the held exponent 1 and the local share
    alone (ratio zero, share 0) gives 0."""
    from probes import growth_form as gf

    r = gf.run(report=_silent)
    assert r["jeans_Mpc"] == pytest.approx(72.7, abs=0.5)
    assert r["held_readings"]["remainder alone"] == 1.0 and r["held_readings"]["local alone"] == 0.0
    assert r["held_readings"]["both"] == pytest.approx(r["held_exponent"], abs=1e-12)
    assert gf.held_exponent(1.0) == 1.0 and gf.held_exponent(0.0) == 0.0
