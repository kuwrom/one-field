"""
test_parameter_ledger.py: freeze the parameter audit's computed facts.

What is frozen: the impedance depth identity, the probes' depth at
kappa = 1 against the sextic that the committed sign's depth 3 needs,
the depth independence of the record area, the factor
between C_EFF and the record closure's xi0/l_Pl, and the equality of
two of the healing conventions on the clock.
"""

import pytest

import parameter_ledger


@pytest.fixture(scope="module")
def audit():
    return parameter_ledger.run(report=lambda *a: None)


def test_impedance_is_the_depth(audit):
    assert audit["P1"] < 1e-12


def test_probes_run_below_the_committed_sign_depth_three(audit):
    assert audit["depth_probe"] == pytest.approx(0.756, abs=1e-3)
    assert audit["kappa_for_range_depth_3"] == pytest.approx(125.9, abs=0.1)
    assert audit["kappa_for_o4_depth_3"] == pytest.approx(15.7, abs=0.1)


def test_gravity_does_not_see_the_depth(audit):
    assert audit["kappa_mu1"] == pytest.approx(audit["kappa_mu4"], rel=1e-6)


def test_c_eff_two_roles_split(audit):
    assert audit["xi0_closure"] / audit["C_EFF"] == pytest.approx(3.86, abs=0.05)
    assert audit["conv_match"] < 1e-12


def test_committed_sign_depth_bounded_only_by_classicality(audit):
    assert audit["depletion_max"] < 0.02
