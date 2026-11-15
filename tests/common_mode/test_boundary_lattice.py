"""Freezes tests/probes/boundary_lattice.py: the boundary reading c on the whole lattice of shares."""

import pytest


def _silent(*args, **kwargs):
    pass


@pytest.fixture(scope="module")
def r():
    from probes import boundary_lattice
    return boundary_lattice.run(report=_silent)


def test_c_is_a_lattice_bounded_below_by_one(r):
    """The 1023 subsets of the primaries give 216 readings and 162
    values of c. The least is 1, attained once, and the next is 36/35,
    1236 standard deviations of the recommended G away."""
    assert (r["n_subsets"], r["n_readings"], r["n_values"]) == (1023, 216, 162)
    assert r["c_min"] == 1.0
    assert r["nearest"] == pytest.approx(36 / 35, abs=1e-9)
    assert r["nearest_sigma"] == pytest.approx(-1236, abs=1)


def test_the_gates_and_symmetry_do_not_select_and_fusion_does(r):
    """With no measured G the wired gates pass 71 values and the six
    conventions of S pass four. Charge conjugation, the one nontrivial
    symmetry of the modular data, keeps 22 values and all four.
    The fusion closed sets give c = 1, 3, 12, 36.  Only the whole set
    contains the fundamental."""
    assert r["n_wired"] == 71
    assert r["six"] == pytest.approx([1.0, 36 / 35, 18 / 17, 12 / 11], abs=1e-9)
    assert len(r["sym_values"]) == 22
    assert r["sym_six"] == pytest.approx([1.0, 36 / 35, 18 / 17, 12 / 11], abs=1e-9)
    assert r["fusion_c"] == pytest.approx([1.0, 3.0, 12.0, 36.0], abs=1e-9)
    assert r["whole"] == pytest.approx(1.0, abs=1e-12)


def test_universality_is_the_input_that_excludes_c3(r):
    """The triality zero set (c = 3) reads only the (c,b,t) channel.
    Only the whole set also reads the leptons and (s,c,b)."""
    reads = r["reads"]
    assert len(reads[1.0]) == 3
    assert reads[3.0] == ["(c,b,t) (h11/K^3)"] and reads[12.0] == [] and reads[36.0] == []
