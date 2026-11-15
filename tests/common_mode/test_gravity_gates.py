"""
test_gravity_gates.py: the gravity lane gate verdicts are computed or witnessed.

What is frozen:
  1. the sign and magnitude gates are computed on the one list that the
     running of alpha_G2 reads (interference/boundary.py).  The bridge
     vents as amount and adds nothing to the running.  The stiffness
     does not depend on the proper time cutoff's shape or placement.
     With the regulator solved from the anchored G, S and M pass at
     every contact.  The bridge read with its searchlight label at the
     boundary's mass gives a repulsive stiffness there under the
     ledger's own cutoff.  The regulator literature leaves that sign
     unsettled.  The counted ledger passes both gates;
  2. every gate names its witness tests, and each named test exists,
     so a verdict cannot outlive the probe that supports it.  U, C and
     P name the tests that contain the case in which they fail;
  3. on the counted ledger the matter rows alone fail the sign gate, so
     the bridge makes it pass.
"""

import ast
import math
import os
from fractions import Fraction

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _test_names(path):
    tree = ast.parse(open(os.path.join(_ROOT, path)).read())
    return {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}


def test_sign_and_magnitude_are_computed_on_the_one_list(res):
    g = res["g"]
    one, gates = g["one_list"], g["gravity_gates"]
    assert gates["S sign"]["status"] == ("PASS" if all(v["inside"] for v in one.values()) else "FAIL")
    assert gates["M magnitude"]["status"] == (
        "PASS" if all(0.6 < v["t0"] < 15.0 for v in one.values()) else "FAIL")
    assert gates["S sign"]["status"] == "PASS" and gates["M magnitude"]["status"] == "PASS"
    assert [one[c]["sigma"] for c in (1.0, 0.0, -1.0)] == pytest.approx([25.8333, 43.8333, 61.8333], abs=1e-4)
    assert all(v["sigma_at_boundary"] == pytest.approx(v["sigma"], abs=1e-9) and v["inside"] for v in one.values())
    assert one[1.0]["t0"] == pytest.approx(4.1115, abs=1e-4) and one[1.0]["t0"] ** -0.5 < 1.0
    lab = g["labelled_bridge"]
    assert lab["sigma_at_boundary"] == pytest.approx(-0.175, abs=1e-3)
    assert lab["t0"] == pytest.approx(0.5885, abs=1e-4) and lab["w_need"] == pytest.approx(0.3702, abs=1e-4)
    assert g["booked_gates"] == {"S sign": "PASS", "M magnitude": "PASS"}


def test_sign_and_magnitude_hold_in_exact_arithmetic():
    """S and M on the one list without floating point: the weights are
    rationals, Sigma = 155/6, 263/6, 371/6 across the contacts, and pi
    lies between 333/106 and 355/113, so Sigma >= 2 pi and
    0.6 < Sigma/(2 pi) < 15 are inequalities between rationals. The
    xi band of the amount is certified the same way."""
    from interference import boundary, root
    lo, hi = Fraction(333, 106), Fraction(355, 113)
    assert lo < Fraction(math.pi) < hi
    weight = {"vector": lambda c: Fraction(1, 6) - Fraction(c, 2), "weyl": lambda c: Fraction(1, 12),
              "bridge": lambda c: Fraction(1, 6)}
    table = boundary.rows(root.d10, root.d11, root.n7, root.n26)
    live = [(n, kind) for _, n, kind, _, reading in table if reading in ("light", "amount")]
    sig = {c: sum(n * weight[kind](c) for n, kind in live) for c in (1, 0, -1)}
    assert sig == {1: Fraction(155, 6), 0: Fraction(263, 6), -1: Fraction(371, 6)}
    for c, s in sig.items():
        assert s >= 2 * hi                                     # S: Sigma >= 2 pi
        assert Fraction(3, 5) * 2 * hi < s < 15 * 2 * lo       # M: 0.6 < Sigma/(2 pi) < 15
        assert float(s) == pytest.approx(boundary.sigma(table, 0.0, float(c), 0.0), abs=1e-12)
    named = {c: s - Fraction(182, 6) for c, s in sig.items()}
    for xi in (Fraction(-178, 1000), Fraction(107, 1000)):     # the band's quoted edges, inside it
        s_xi = {c: named[c] + 182 * (Fraction(1, 6) - xi) for c in named}
        assert all(s >= 2 * hi and Fraction(3, 5) * 2 * hi < s < 15 * 2 * lo for s in s_xi.values())
    for xi in (Fraction(-179, 1000), Fraction(108, 1000)):     # one step outside, a gate fails
        s_xi = {c: named[c] + 182 * (Fraction(1, 6) - xi) for c in named}
        assert not all(s >= 2 * lo and s < 15 * 2 * hi for s in s_xi.values())


def test_the_running_and_the_stiffness_read_one_list(res):
    from interference import boundary, root
    table = boundary.rows(root.d10, root.d11, root.n7, root.n26)
    assert boundary.b0(table) == root.b0_G2
    sig = boundary.sigma(table, 1.0, 1.0, float(root.xi_bridge))
    assert sig == pytest.approx(res["g"]["one_list"][1.0]["sigma_at_boundary"], abs=1e-12)
    assert sig >= 2 * math.pi
    bridge = [row for row in table if row[0] == "bridge"][0]
    assert bridge[3] == 0 and bridge[4] == "amount"
    labelled = boundary.rows(root.d10, root.d11, root.n7, root.n26, bridge="light")
    assert boundary.b0(labelled) == boundary.b0(table) - Fraction(13, 3)


def test_every_gate_names_existing_witness_tests(res):
    gates = res["g"]["gravity_gates"]
    assert set(gates) == {"S sign", "U universality", "M magnitude", "C one cone",
                          "T tensor carrier", "P no preferred location"}
    names = {}
    for gate in gates.values():
        assert gate["status"] == "PASS" and gate["tier"] == "derived"
        assert len(gate["witness"]) >= 1
        for w in gate["witness"]:
            path, func = w.split("::")
            names.setdefault(path, _test_names(path))
            assert func in names[path], f"witness {w} does not exist"


def test_universality_cone_and_location_name_their_failing_case(res):
    """U and C are derived on the saddle's stabilizer and name the test
    in which a frame vector breaks them.  U and P name the two lump
    exchange test."""
    gates = res["g"]["gravity_gates"]
    frame = "tests/common_mode/test_saddle_stabilizer.py::test_saddle_stabilizer_a_frame_gives_each_row_its_own_speed"
    lumps = "tests/common_mode/test_two_lumps.py::test_two_lumps_attract_with_their_whole_energy_as_weight"
    assert frame in gates["U universality"]["witness"] and frame in gates["C one cone"]["witness"]
    assert lumps in gates["U universality"]["witness"] and lumps in gates["P no preferred location"]["witness"]


def test_matter_rows_alone_fail_the_sign_gate(res):
    g = res["g"]
    assert g["Sigma_total_UV"] > 0
    assert g["Sigma_total_UV"] - g["Sigma_bridge"] < 0
