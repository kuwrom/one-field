"""Freezes tests/probes/strong_field.py: the static spherical vacuum
solutions of the traceless equations, their horizon, and the mass
above which the two derivative action holds there."""

import pytest

sp = pytest.importorskip("sympy", reason="needs sympy (pip install -r requirements.txt)")


def _silent(*args, **kwargs):
    pass


def test_strong_field_schwarzschild_de_sitter_with_two_integration_constants():
    """The equations force h = f and reduce to r^2 f'' - 2 f + 2 = 0,
    solved by f = 1 - 2 G M/r - Lambda r^2/3 with R = 4 Lambda. The
    horizon is regular for T = 1/(8 pi G M) at Lambda = 0, and
    dM = T dS holds with the entropy A/(4G) of the same G.  The
    horizon's curvature scale is 0.931 M_b^2/M.  It reaches the kernel's
    10 and 50 percent departures at 7.2 and 2.2 boundary masses."""
    from probes import horizon_amount as ha
    from probes import strong_field as sf

    r = sf.run(report=_silent)
    assert r["constants"] == 2 and r["scalar_curvature"] == "4 Lambda" and r["first_law"]
    assert r["temperature"] == "1/(8*pi*G*M)"
    assert r["horizon_scale"] == pytest.approx(0.9306, abs=1e-4)
    assert r["mass_floors"] == pytest.approx((7.21, 2.16), abs=0.01)
    assert r["solar_ratio"] < 1e-37
    assert ha.run(report=_silent)["q_departure"] == pytest.approx(list(sf.DEPARTURES), abs=1e-3)


def test_strong_field_other_static_metrics_leave_a_residual():
    """A linear term in f (the conformal gravity potential) and the
    screened Poisson potential with the healing length
    xi0 = l_Pl/d10 each fail the vacuum equations. The
    Schwarzschild de Sitter f does not."""
    from interference import root
    from probes import strong_field as sf

    G, M, Lam, gam = sp.symbols("G M Lambda gamma", positive=True)
    r = sf.r
    f = 1 - 2 * G * M / r
    xi0 = sp.sqrt(G) / root.d10                       # l_Pl/d10 with hbar = c = 1
    assert not sf.residual(f) and not sf.residual(f - Lam * r**2 / 3)
    assert sf.residual(f + gam * r)
    assert sf.residual(1 - 2 * G * M * sp.exp(-r / xi0) / r)
