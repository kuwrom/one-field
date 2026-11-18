"""Freezes tests/probes/neutrino_mass_term.py: the mass term the two F4 singlets close on."""

import pytest


@pytest.fixture(scope="module")
def r():
    from probes import neutrino_mass_term
    return neutrino_mass_term.run(report=lambda *a: None)


def test_one_structure_passes_both_gates(r):
    """N1 to N4: of the eight mass structures, three conserve the lepton
    label, six give two massive neutrinos, and one does both: the
    singlet doublet term with the singlets carrying the lepton label.
    Neither gate decides alone."""
    assert r["counts"] == {"menu": 8, "label": 3, "data": 6, "both": 1}
    assert r["survivor"] == {"terms": "b", "label": 1, "massive": 2, "massless": 1, "kind": "Dirac"}


def test_the_seesaw_fails_for_every_label(r):
    """N2: the singlet mass needs label 0 and the singlet doublet term
    needs label 1, so no label admits both."""
    assert r["seesaw_labels"] == []


def test_the_data_do_not_separate_the_rows_today(r):
    """N5 to N7: the survivor has m_bb = 0. The seesaw row differs in
    m_bb alone and sits 7.6 times below the strongest limit. The right
    handed halves never share the early bath."""
    assert r["m_bb"] == 0.0
    assert r["m_bb_seesaw"] == pytest.approx((1.51, 3.70), abs=0.05)
    assert r["m_bb_margin"] == pytest.approx(7.6, abs=0.1)
    assert r["sum"] == pytest.approx(58.78, abs=0.05) and r["m_beta"] == pytest.approx(8.81, abs=0.05)
    assert r["yukawa"] == pytest.approx(2.88e-13, rel=0.01) and r["bath"] < 1e-9
