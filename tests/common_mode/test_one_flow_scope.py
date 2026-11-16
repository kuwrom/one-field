"""Freezes tests/probes/one_flow_scope.py: "coherence carries one flow" is the initial state."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_the_cold_count_folds_into_dispersion_without_pressure():
    """Before the first fold the count is one stream at every place. At
    s = 1 the single valued velocity gradient diverges. At s = 3 the
    count holds three streams over 0.51 of the period, carrying 0.87 of
    the mass, with a local velocity dispersion of 0.49 and no pressure.
    One velocity at every place would carry none."""
    from probes import one_flow_scope

    out = one_flow_scope.run(report=_silent)
    assert out["fold_fraction"] == pytest.approx(0.51, abs=0.01)
    assert out["fold_mass"] == pytest.approx(0.87, abs=0.01)
    assert out["dispersion"] == pytest.approx(0.49, abs=0.01)
    assert out["single_stream_dispersion"] < 1e-3
