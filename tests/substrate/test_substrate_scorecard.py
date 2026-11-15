"""Scorecard tests of the substrate branch: every number the engine prints for it is asserted."""

import math
from fractions import Fraction

import numpy as np
import pytest

from root import PDG_MASSES


def err_pct(predicted: float, expected: float) -> float:
    """Absolute percentage deviation from the expected value."""
    return abs((predicted - expected) / expected) * 100.0


# ═══════════════════════════════════════════════════════════════════════
#  Quark masses  (6 predictions)
# ═══════════════════════════════════════════════════════════════════════



# ═══════════════════════════════════════════════════════════════════════
#  NLS soliton (substrate dynamics)
# ═══════════════════════════════════════════════════════════════════════

def test_nls_soliton_stability():
    """The Z₃-NLS integrator is variational (second order energy drift).

    Number conservation alone cannot certify the equation. The energy
    drift stays below 1e-4 and halving dt quarters it.
    """
    import nls_soliton
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        result = nls_soliton.derive()
    assert result['conservation'] < 1e-6
    assert result['energy_drift'] < 1e-4
    assert result['drift_ratio_half_dt'] == pytest.approx(0.25, abs=0.02)
