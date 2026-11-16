"""Freezes tests/probes/gates_as_records.py: the six gates as statements about records, each with a failing case."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_each_gate_is_a_record_with_a_failing_case():
    """A lump's rate is its weight per added quantum.  A bias moves the
    rates of two kinds of lump by the same fraction, so their ratio
    records nothing.  A bias that depends on the count moves the ratio
    by 0.297 of itself. One symbol on the saddle, arrivals apart with a frame
    vector. The stiffness is 25.83 with the amount and -4.5 without.
    The grain is 4.11, 6.98, 9.84 boundary units, and the labelled
    bridge at the boundary keeps b0 and fails S and M. Two patterns
    and no breathing, and the labelled list's kernel has a zero at
    0.624 M_b, inside the boundary.  gamma = 1.  The stiffness depends
    on the rows' counts and not on their masses.  If the amount is read
    as massive below the boundary, its stiffness depends on that mass
    at order one."""
    from probes import gates_as_records

    r = gates_as_records.run(report=_silent)
    assert r["rates"] == pytest.approx((0.91006, 0.71644), abs=1e-4)
    assert r["rate_per_quantum"] == pytest.approx(r["rates"], rel=1e-3)
    assert r["shift"] == pytest.approx((1.0, 1.0), abs=1e-3) and abs(r["ratio_moves"]) < 1e-3
    assert r["shift_count"] == pytest.approx((1.0988, 1.3958), abs=1e-3)
    assert r["ratio_moves_count"] == pytest.approx(-0.297, abs=1e-3) and r["rows_frame"] > 1e-2
    assert r["arrivals_apart"] == 0.0 and min(r["arrivals_apart_frame"]) > 1e-3
    assert r["sigma"] == pytest.approx(25.8333, abs=1e-3) and r["sigma_bare"] == pytest.approx(-4.5, abs=1e-9)
    assert (r["bias_sign"], r["bias_sign_bare"]) == (-1, 1)
    assert r["grain"] == pytest.approx([4.11, 6.98, 9.84], abs=0.01) and r["readings_boundary"] == (True, False, False)
    assert r["patterns"] == 2 and r["breathing"] < 1e-12
    assert r["zero_labelled"] == pytest.approx(0.6244, abs=1e-3)
    assert r["gamma"] == 1.0 and r["stress_share"] < 1e-6
    assert r["mass_response_one"] == 0.0 and r["mass_response_massive"] < -0.5
