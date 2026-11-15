"""Freezes tests/probes/amplitude_gate.py: B/A against the angle menu, no measurement read."""

import pytest


@pytest.fixture(scope="module")
def r():
    from probes import amplitude_gate
    return amplitude_gate.run(report=lambda *a: None)


def test_the_window_holds_the_quantum_trace_and_not_the_classical(r):
    """The menu's unique survivor 2/9 exists for B in [1.30541, 1.47367].
    The quantum trace sqrt 2 gives one survivor, the classical 1/sqrt 2
    twelve."""
    assert r["window"] == pytest.approx((1.30541, 1.47367), abs=2e-5)
    assert r["n_quantum"] == 1 and r["n_classical"] == 12
    assert r["width"] == pytest.approx(0.059, abs=0.001) and r["held"] < 1e-9
