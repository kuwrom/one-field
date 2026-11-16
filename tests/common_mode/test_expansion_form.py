"""Freezes tests/probes/expansion_form.py: the traceless equations on a
homogeneous isotropic metric give the rate equation with Lambda its
integration constant."""

import pytest

pytest.importorskip("sympy")


def _silent(*args, **kwargs):
    pass


def test_expansion_form_rate_equation_with_lambda_as_integration_constant():
    """Every component of the traceless equations is the same equation
    dH/dt - k/a^2 = -4 pi G (rho + p). With a conserved source
    3 H^2 + 3 k/a^2 - 8 pi G rho is a constant of the motion, fixed by
    the initial data and different for different initial rates.  Taken
    as a density it has w = -1. An amount with the CMB bound's w = 0.0024
    would move its ratio to the baryons by 5.2 percent since z = 1100."""
    from probes import expansion_form as ef

    r = ef.run(report=_silent)
    assert (r["components"], r["independent"]) == (4, 1) and r["constant"]
    assert r["lambdas"][1.0] == pytest.approx(0.6999, abs=1e-10)
    assert r["lambdas"][1.2] == pytest.approx(1.1399, abs=1e-10)
    assert r["drift"] < 1e-9 and r["form_error"] < 1e-9 and r["w_error"] < 1e-7
    assert r["ratio_drift_at_bound"] == pytest.approx(0.0517, abs=5e-4)


def test_expansion_form_a_source_that_is_not_conserved_leaves_no_constant():
    """With d rho/dt = -3 H (rho + p) + Q the same combination moves at
    -8 pi G Q: the constant exists only for a conserved source."""
    import sympy as sp
    from probes import expansion_form as ef

    rate, trace, Q, G = ef.constant_of_motion(ef.trace_free_components())
    assert sp.simplify(rate.subs(Q, 0)) == 0 and trace == 0
    assert sp.simplify(rate + 8 * sp.pi * G * Q) == 0 and sp.simplify(rate.subs(Q, 1)) != 0
