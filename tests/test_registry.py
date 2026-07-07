"""Registry <-> engine coupling: the ledger cannot drift from canon.

Every numerical value registered in interference/registry.py is
asserted against the solved engine.  If a canonical value moves
(deliberately or not), this file fails: updating the registry becomes
a forced, conscious act with a diff, never a forgotten edit.  This is
how the ledger stays honest: registered values are read back from the engine, never restated by hand.

Run with:
    pytest -q
"""

import math
import os
import sys
from fractions import Fraction

import pytest

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "interference")))

import registry


def _watch(i):
    return next(w for w in registry.WATCHES if w["id"] == i)


def _pred(i):
    return next(p for p in registry.PREDICTIONS if p["id"] == i)


def test_watches_match_engine(res):
    """Each registered watch value equals the engine's canon."""
    assert _watch(1)["value"] == pytest.approx(
        res["c"]["alpha_s_MZ_thresh"], abs=5e-5)
    assert _watch(2)["value"] == pytest.approx(
        res["c"]["inv_alpha_phys"], abs=5e-9)
    assert _watch(3)["value"] == pytest.approx(
        res["m"]["ms_over_mud"], abs=5e-3)
    assert _watch(4)["value"] == pytest.approx(
        res["g"]["sin2W"], abs=5e-6)
    assert _watch(5)["value"] == pytest.approx(res["m"]["m_b"], abs=0.05)
    assert _watch(6)["value"] == pytest.approx(
        res["d"]["DM_baryon_ratio"], abs=5e-4)


def test_predictions_match_engine(res):
    """Registered exposures equal the engine's outputs."""
    assert _pred("A")["value"] == pytest.approx(
        res["mx"]["delta_CP_deg"], abs=0.05)
    assert _pred("B")["value"] == pytest.approx(
        res["c"]["inv_alpha_phys"], abs=5e-9)
    assert _pred("E")["value"] == pytest.approx(2 * math.pi - 1, abs=5e-4)
    assert _pred("G")["value"]["mu_over_md"] == Fraction(38, 83)
    assert _pred("G")["value"]["Q_ellipse"] == pytest.approx(
        res["m"]["Q_ellipse"], abs=5e-3)
    assert _pred("F")["value"]["v_EW"] == pytest.approx(
        res["m"]["v_EW_pred_GeV"], abs=5e-5)
    assert _pred(9)["value"] == pytest.approx(2 * math.sqrt(2), abs=1e-12)


def test_promotion_matches_engine(res):
    """The promoted Higgs edge: registry record == canon, exactly."""
    from root import WEB, h10
    pr = registry.PROMOTIONS[0]
    assert pr["lambda_MPl"] == pytest.approx(
        WEB.state["lambda_MPl"], abs=1e-8)
    assert pr["mH"] == pytest.approx(res["g"]["mH_pred"], abs=5e-3)
    # the vent factor in the record is the one in the web
    undeflated = WEB.state["lambda_MPl"] / float(1 - h10)
    assert undeflated == pytest.approx(-0.016007340, abs=1e-8)


def test_derivation_programs_reference_real_witnesses():
    """Every program's witness test exists in this suite."""
    import test_coverage
    for prog in registry.DERIVATION_PROGRAMS:
        fn = prog["witness"].split("::")[-1]
        assert hasattr(test_coverage, fn), f"missing witness {fn}"


def test_registry_renders():
    text = registry.render()
    for key in ("WATCHES", "PROMOTIONS", "DERIVATION PROGRAMS",
                "(1 - h10)", "125.2965"):
        assert key in text
