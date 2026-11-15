"""Freezes tests/audit/gate_set.py: every closure gate, what it removes, and the premise it states."""

import pytest

import gate_set


def _silent(*a, **k):
    pass


@pytest.fixture(scope="module")
def out():
    return gate_set.run(report=_silent)


def test_each_stage_closes_on_one_survivor(out):
    g = out["gates"]
    assert g[("lattice", "even")]["survivors"] == [("E8", 8)]
    assert g[("pair", "emergent Z3")]["survivors"] == [("G2(1) x F4(1)",)]
    assert g[("level", "h_J = 1")]["survivors"] == [(3,)]
    assert g[("coherence", "(f) trivial centres")]["survivors"] == [("E8", "G2", "F4")]


def test_what_each_gate_removes(out):
    """Five gates select alone, six share the selection in pairs, and five
    are consequences: the survivor of the other gates already satisfies
    them.  Each selecting gate names the rival it removes."""
    g = out["gates"]
    status = {k: v["status"] for k, v in g.items()}
    assert {k for k, s in status.items() if s == "selects"} == {
        ("lattice", "unimodular"), ("lattice", "minimal_rank"), ("pair", "emergent Z3"),
        ("level", "h_J = 1"), ("coherence", "(f) trivial centres")}
    assert {k for k, s in status.items() if s == "shares"} == {
        ("lattice", "even"), ("lattice", "irreducible"), ("pair", "non-pointed"),
        ("level", "modular extension"), ("coherence", "(c) every mode assigned"),
        ("coherence", "(e) self-adjoint parent")}
    assert g[("lattice", "unimodular")]["rivals"] == [("A8", 8), ("D8", 8)]
    assert g[("lattice", "minimal_rank")]["rivals"] == [("D16+", 16)]
    assert g[("level", "h_J = 1")]["rivals"] == [(6,), (9,), (12,)]
    assert g[("coherence", "(f) trivial centres")]["rivals"] == [("E8", "A1", "E7")]


def test_every_gate_names_its_premise(out):
    """Every gate names its premise: R1 for the lattice axioms, the level
    and the bridge current, self-presence for the parent, and quantum
    numbers crystallising on breaking for the emergent Z3 and the centres."""
    assert out["unstated"] == []
    assert all(r["premise"] for r in out["gates"].values())
