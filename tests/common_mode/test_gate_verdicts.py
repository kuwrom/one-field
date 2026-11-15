"""Freezes tests/audit/gate_verdicts.py: the engine's gate statuses equal the verdicts computed from their deciders."""

import os
import sys

import pytest

pytest.importorskip("sympy")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "audit"))
import gate_verdicts  # noqa: E402


@pytest.fixture(scope="module")
def computed():
    return gate_verdicts.run(report=lambda *a, **k: None)


def test_engine_statuses_equal_the_computed_verdicts(res, computed):
    gates = res["g"]["gravity_gates"]
    assert set(gates) == set(gate_verdicts.NAMES)
    for name in gate_verdicts.NAMES:
        assert gates[name]["status"] == computed["verdicts"][name] == "PASS", name


def test_each_decider_fails_on_its_failing_case(computed):
    """Every gate with a computed rival reading fails on it.  A frame
    vector in the vacuum fails C and U and leaves S and T.  The list
    without the amount fails S and M.  The amount carrying the label's
    coupling has a kernel zero inside the boundary and fails T.  The
    amount read as massive below the boundary passes S and M.  But its
    stiffness depends on that mass at order one, so it fails P."""
    framed = computed["framed"]
    assert framed["C one cone"] == framed["U universality"] == "FAIL"
    assert framed["S sign"] == framed["T tensor carrier"] == "PASS"
    fail = computed["failing"]
    assert set(fail) == set(gate_verdicts.NAMES)
    assert all(verdict == "FAIL" for _, verdict in fail.values())
    assert fail["S sign"][0] == fail["M magnitude"][0] == "the list without the amount"
    assert fail["C one cone"][0] == fail["U universality"][0] == "a frame vector in the vacuum"
    assert fail["T tensor carrier"][0] == "the amount with the label's coupling"
    assert fail["P no preferred location"][0] == "the amount read as massive below the boundary"
    amount = computed["evidence"]["amount"]
    assert amount["mass_response_one"] == 0.0 and amount["mass_response_massive"] == pytest.approx(-0.988, abs=1e-3)
    assert computed["evidence"]["amount"]["q_zero_label"] == pytest.approx(1.295, abs=2e-3)
