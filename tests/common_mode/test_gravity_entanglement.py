"""Freezes tests/probes/gravity_entanglement.py: the framework's prediction for gravitationally induced entanglement."""

import pytest


@pytest.fixture(scope="module")
def r():
    from probes import gravity_entanglement
    return gravity_entanglement.run(report=lambda *a: None)


def test_the_masses_entangle_at_the_newtonian_phase(r):
    """Phi = -0.314, concurrence 0.156; G's two readings move it by
    under 1e-5; a separable mediator gives 0; Diosi-Penrose 6.6 ms."""
    assert r["Phi"] == pytest.approx(-0.314, abs=1e-3)
    assert r["concurrence"] == pytest.approx(0.156, abs=1e-3) and abs(r["C_shift"]) < 1e-5
    assert r["separable"] < 1e-12
    assert r["tau_dp"] == pytest.approx(6.58e-3, rel=1e-2)
