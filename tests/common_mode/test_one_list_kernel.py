"""Freezes tests/probes/one_list_kernel.py: every row of a list as one kernel of the composite metric."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_one_list_kernel_carries_the_boundary_weights():
    """A heavy scalar row is the linearized Einstein-Hilbert operator
    with stiffness (1/6 - xi) E2(m^2 t0) per component, the weight
    boundary_content.py assumes.  The counted list's spin-2 form factor
    has no zero on the axis.  The labelled list's kernel has
    boundary_content's stiffness.  With the regulator at the boundary it
    is -0.175, with a ghost pole and a repulsive potential.  At the
    self-consistent t0 = 0.5885 it is 3.698 and attractive, with its
    spin-2 zero at 0.624 M_b, inside the boundary.  The carrier alone
    has no zero.  As m^2 t0 -> 0 it departs from the Einstein value at
    q^2 t0 = 0.26 (10 percent) and 3.3 (50 percent).  Heavy carriers
    depart at q^2 t0 of order one too, so the departure sits at the
    regulator scale for every mass."""
    from probes import one_list_kernel as olk

    r = olk.run(report=_silent)
    assert r["direct"] < 1e-10
    assert all(v["worst"] < 1e-5 and v["shape"] < 5e-5 for v in r["heavy"].values())
    assert r["heavy"][1.0]["E2"] == pytest.approx(0.148496, abs=1e-6)
    assert r["sigma_booked"] == pytest.approx(27.0, abs=1e-3)
    assert r["zero_booked"] is None
    assert r["regulator_booked"] == pytest.approx(0.4824, abs=1e-3)
    assert r["q_departure_booked"] == pytest.approx([0.1295, 0.4329], abs=1e-3)
    at_b, at_s = r["labelled_list"]["boundary"], r["labelled_list"]["self"]
    assert r["agree"] < 1e-4 and r["shape_labelled"] < 1e-3
    assert at_b["sigma"] == pytest.approx(-0.17486, abs=1e-4)
    assert at_b["sign"] == pytest.approx(-1.0, abs=1e-3) and at_b["negative"] and at_b["zero"] is None
    assert r["t0_self"] == pytest.approx(0.58853, abs=1e-4)
    assert at_s["sigma"] == pytest.approx(3.6978, abs=1e-3) and at_s["sign"] == pytest.approx(1.0, abs=1e-3)
    assert r["zero_labelled"] == pytest.approx(0.2295, abs=1e-3)
    assert r["q_zero_labelled"] == pytest.approx(0.6244, abs=1e-3)
    assert r["regulator_self"] == pytest.approx(1.3035, abs=1e-3)
    assert r["f0_zero_boundary"] == pytest.approx(0.0431, abs=5e-4)
    assert r["q_departure_labelled"] == pytest.approx([0.1366, 0.3797], abs=1e-3)
    dep = r["departure"]
    assert dep[0.0] == pytest.approx([0.2592, 3.3165], abs=1e-3)
    assert dep[0.01] == pytest.approx([0.3146, 3.5223], abs=1e-3)
    assert dep[0.1] == pytest.approx([0.4550, 4.2141], abs=1e-3)
    assert max(r["shift"][1e-4]) < 1e-2
    assert [v[0] for v in r["departure_heavy"].values()] == pytest.approx([0.737, 0.994, 1.068], abs=2e-3)
