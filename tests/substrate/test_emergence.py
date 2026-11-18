"""
test_emergence.py: freeze the emergent clock protocol and its verdict.

The probe (tests/probes/emergent_clock.py) is the substrate program's
first experiment, run forward: random fluctuations go in, the spectrum
is read out, and no target is used anywhere.  This file freezes the
STRUCTURE and the VERDICT.  It does not freeze any closeness to the
Brannen spectrum, which the probe is explicitly forbidden to assert.

The broadband verdict is frozen: the substrate only gas sign regime
reads broadband clocks, and the clock lives in the lump on the
self-bound sign (lump_formation.py).
"""

import pytest


def _silent(*args, **kwargs):
    pass


@pytest.mark.slow
def test_emergent_clock_protocol_and_verdict():
    from probes import emergent_clock
    out = emergent_clock.run(report=_silent)

    # linearised μ = (3g₀ + 2g₁)ρ₀ at ρ₀ = 1
    mu_lin = 3.0 + 2.0 * emergent_clock.G1_OVER_G0

    assert len(out["results"]) == len(emergent_clock.SEEDS)
    for r in out["results"]:
        # structural, target blind
        assert r["conservation"] < 1e-9
        assert r["growth"] > 3.0, \
            "knots do not grow from generic noise"
        assert all(0.0 < o < 2.0 * mu_lin for o in r["omegas"]), \
            "dominant line outside the physical band"
        assert len(r["peak_shares"]) == 3

    # The verdict freeze: substrate only 1D is broadband.
    assert out["broadband"] is True
    # Q is a readout; 1/3 <= Q < 1 holds identically for any three
    # positive lines and is not asserted.


@pytest.mark.slow
def test_emergent_clock_3d_rung():
    """The 3D rung of the protocol (one seed, several minutes).  Knots
    form from generic noise in 3D too.  The substrate only 3D clock
    verdict is frozen as broadband, matching 1D."""
    from probes import emergent_clock
    out = emergent_clock.run_3d(report=_silent)
    mu_lin = 3.0 + 2.0 * emergent_clock.G1_OVER_G0   # BdG μ
    assert out["growth"] > 3.0
    assert out["conservation"] < 1e-9
    assert all(0.0 < o < 2.0 * mu_lin for o in out["omegas"]), \
        "dominant line outside the physical band"
    assert out["broadband"] is True            # substrate only 3D
    # Q is a readout. Its 1/3..1 band holds identically, not asserted
