"""Freezes tests/probes/two_lumps.py: two lumps attract through the induced kernel, each weighing its energy."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_two_lumps_attract_with_their_whole_energy_as_weight():
    """The profile's own fields give int T_00 = E, a conserved stress
    and int T_kk = 0. Exchanging the one list's kernel, two lumps with
    binding +0.25% and -25% have V r/(G E1 E2) = -1 at three
    separations.  In units of the charges' rest mass it is -0.752.  The
    stress terms move it by less than 1e-6. The coefficient is the
    kernel's: F0/F2 = -2.2 gives -1.0303. The full form factors at
    eps = 0.01 and 0.1 move the potential by less than 1e-4."""
    from probes import two_lumps

    r = two_lumps.run(report=_silent)
    assert r["binding"] == pytest.approx((0.00246, -0.24964), abs=1e-4)
    assert r["sources"] < 1e-6 and r["conservation"] < 1e-4 and r["laue"] < 1e-6
    assert r["newton"] == pytest.approx([-1.0, -1.0, -1.0], abs=1e-6)
    assert r["charge_weight"] == pytest.approx(-0.7522, abs=1e-3)
    assert r["stress_share"] < 1e-6
    assert r["moved_kernel"] == pytest.approx(-2 * (2 / 3 - 1 / 6.6), abs=1e-6)
    assert max(r["reach"].values()) < 1e-4 and min(r["kernel_floor"]) > 0
