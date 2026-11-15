"""Freezes tests/probes/quadrupole_radiation.py: a conserved slow source
radiates the quadrupole power through the kernel's pole, with Newton's G."""

import numpy as np
import pytest


def _silent(*args, **kwargs):
    pass


def test_quadrupole_radiation_two_polarizations_and_the_quadrupole_coefficient():
    """The kernel's F0/F2 = -2 makes the on shell exchange of a conserved
    source the squared transverse traceless part of its stress. Its
    angular integral keeps the traceless quadrupole alone. With the
    constant fixed by Newton's law the power is (G/5) <Q''' Q'''>, and
    32/5 G mu^2 d^4 W^6 for a circular binary."""
    from probes import quadrupole_radiation as qr

    r = qr.run(report=_silent)
    assert r["kernel_ratio"] == pytest.approx(-2.0, abs=3e-6)
    assert r["tt_error"] < 1e-12 and r["angular_error"] < 1e-12 and r["quadrupole_error"] < 1e-12
    assert abs(r["breathing"]) < 1e-12 and r["binary"] == pytest.approx(32 / 5, abs=1e-12)
    assert abs(r["einstein_pull"]) < 2.0


def test_quadrupole_radiation_a_source_not_conserved_reads_no_power():
    """A source that is not conserved gives the exchange both signs, so
    no power can be read. The power scales as the square of the source,
    and at F0/F2 = -2 a breathing sphere sends nothing."""
    from probes import quadrupole_radiation as qr

    r = qr.run(report=_silent)
    assert 0 < r["negative_share"] < 1
    grid = qr.directions()
    w, Qw = qr.binary()
    half = qr.power(w, Qw, 0.5, 1.0, grid)
    assert qr.power(w, 2 * Qw, 0.5, 1.0, grid) == pytest.approx(4 * half, rel=1e-12)
    assert qr.angular(np.eye(3, dtype=complex), qr.a_of(-2.0), grid) == pytest.approx(0.0, abs=1e-12)
