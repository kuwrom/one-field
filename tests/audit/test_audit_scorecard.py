"""Scorecard tests of the audit branch: every number the engine prints for it is asserted."""

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



def test_canonical_freeze(res):
    """The frozen canonical table (ledger closed)."""
    m, g, c = res["m"], res["g"], res["c"]
    import root as _r

    # exact (dimensionless web + anchor inversion)
    assert abs(_r.inv_alpha_phys - 137.0359990495613) < 1e-10
    assert abs(m["m_mu"] - 105.65837578191945) < 1e-9
    assert abs(m["m_tau"] - 1776.9092813590157) < 1e-7
    assert abs(m["v_EW_pred_GeV"] - 246.21964503014013) < 1e-9
    assert abs(_r.G_PRED - 6.673935680964234e-11) < 1e-22
    assert abs(_r.M_Pl_GeV - 1.220917145099806e+19) < 1e10

    # regulator diagnostic 6pi/Sigma (no G is read from the ledger)
    assert abs(g["G_ratio_UV"] - 0.698131701) < 1e-8
    assert abs(g["sin2W"] - 0.23128541) < 1e-7

    # EW transport chain (PDG 2024 imports)
    assert abs(g["M_W_derived"] - 80.3559) < 2e-3
    assert abs(g["M_Z_derived"] - 91.1885) < 2e-3

    # couplings + Higgs (declared imports, banded)
    # alpha_s matched at mu* = M_Pl e^{-(9pi^2/2-6)} = v_EW e^{15/512}
    # (the endpoint of the gauge lever arm)
    assert abs(c["mu_star_GeV"] - 253.534) < 1e-2
    assert abs(c["alpha_s_MZ_thresh"] - 0.118385) < 1e-5
    assert abs(g["lambda_mt_transport"] - 0.1265078) < 2e-7

    # words lemma bases (integer ranks, exact)
    import words as _w
    import io as _io, contextlib as _ctx
    with _ctx.redirect_stdout(_io.StringIO()):
        bases = _w.derive()["base"]
    assert bases == (4, 12, 97)


def test_frozen_web_nodes(res):
    """The graph nodes' fixed point values, frozen."""
    from root import WEB
    st = WEB.state
    assert abs(st["lambda_MPl"] - (-0.012450153)) < 1e-8  # ×(7/9): fundamental share vent
    assert abs(st["G_ratio"] - 0.698131701) < 1e-8
    assert abs(st["dark_ratio"] - (2.0 * math.pi - 1.0)) < 1e-12
    assert abs(st["v_EW_GeV"] - 246.219645) < 1e-5
    assert abs(st["M_Pl_MeV"] / 1e3 - 1.220917145e19) < 1e11


# ═══════════════════════════════════════════════════════════════════════
#  Attractor (Chladni) property of the web
# ═══════════════════════════════════════════════════════════════════════

def test_web_attractor(res):
    """Perturb every web node by up to ±50% and solve again: the whole web
    returns to its fixed point (machine precision, 5/5 trials).  The
    constants are stored nowhere.  They are the fixed point the
    recursion reaches, as a Chladni figure is the pattern a plate
    settles into."""
    import random
    from root import WEB

    WEB.solve()                      # ensure the state IS the fixed point
    baseline = dict(WEB.state)
    random.seed(3)
    try:
        for _ in range(5):
            WEB.state = {k: v * (1.0 + random.uniform(-0.5, 0.5))
                         for k, v in baseline.items()}
            WEB.solve(iters=2000)
            dev = max(abs(WEB.state[k] / baseline[k] - 1.0)
                      for k in baseline if baseline[k] != 0.0)
            assert dev < 1e-10, f"web did not return: max dev {dev:.2e}"
    finally:
        WEB.state = dict(baseline)
        WEB.solve()
