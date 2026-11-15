"""The dust gate at the classical level (tests/probes/dust_gate.py).

What is frozen: no name free constant 3-form exists on the carrier's
field space.  o(6) and u(3) have none.  su(3) has two, but they are
charged under the overall phase and odd under relabelling.  A name free
identically conserved current is a superpotential divergence.  The
nonexact current needs a winding profile.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "probes"))
import dust_gate as dg  # noqa: E402


@pytest.fixture(scope="module")
def res():
    return dg.run(report=lambda *a, **k: None)


def test_no_name_free_constant_three_form(res):
    d1 = res["d1"]
    assert (d1["o6"], d1["u3"], d1["su3"]) == (0, 0, 2)
    assert d1["phase_charge"] == pytest.approx(3.0)
    assert d1["odd_relabel"] == pytest.approx([-1.0, -1.0])


def test_conserved_currents_are_superpotentials(res):
    assert res["d2"]["divergence"] < 1e-5
    assert res["d2"]["superpotential"] < 1e-5


def test_non_exact_current_needs_a_winding_profile(res):
    assert abs(res["d3"]["winding_density"]) > 1e-3
    assert res["d3"]["no_winding"] == 0.0


def test_induced_route_needs_a_massive_cold_amount():
    import dust_gate_state as ds
    r = ds.run(report=lambda *a, **k: None)
    assert r["s1"]["w_rad"] == pytest.approx(1 / 3)
    assert r["s1"]["window_dlna"] == pytest.approx(0.0072, abs=1e-6)
    assert r["s2"]["v_max_cs"] == pytest.approx(3.098e-3, abs=1e-5)
    assert r["s3"]["per_scalar"] == pytest.approx(0.0268, abs=5e-4) and r["s3"]["total"] == pytest.approx(4.87, abs=0.05)
