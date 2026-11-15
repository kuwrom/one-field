"""Freezes tests/probes/gravity_resolution.py: gravity's two readings and
their resolutions.  The anchor is a reading at 22 ppm.  No input and
no deeper echo of its ladder can supply its residual.  The ledger's G
at the anchored regulator is an identity.  The ledger alone gives only
a sign and a size range of about two decades.  The stiffness nodes
are read at depth one and the common mode at depth four."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_gravity_resolution_anchor_reads_ledger_silhouettes(res):
    from probes import gravity_resolution as gr

    out = gr.run(report=_silent, build=False)
    lad = out["ladder_ppm"]
    assert [lad[1], lad[2], lad[4]] == pytest.approx([-1139.9, 21.5, 22.2], abs=0.1)
    at = out["attribution"]
    assert at["residual_over_sigma"] == pytest.approx(2.02, abs=0.02)
    assert at["delta_S"] == pytest.approx(2.22e-5, rel=1e-2)
    assert at["residual_over_next_echo"] > 1e4
    rd = out["readings"]
    assert rd["anchor_response"] == pytest.approx(1.0, abs=1e-6)
    assert rd["ledger_anchored_response"] == 0.0
    assert rd["ledger_placed_response"] == pytest.approx(-0.5, abs=1e-12)
    sh = out["silhouette"]
    assert sh["all_attractive"]
    assert (sh["sigma_min"], sh["sigma_max"]) == pytest.approx((109 / 6, 371 / 6), abs=1e-9)
    assert sh["decades"] == pytest.approx(1.93, abs=0.01)
    assert out["echoes"]["anchor_web_depth"] == 4 and out["echoes"]["common_mode_web_depth"] == 4
    assert all(n["max_term_depth"] <= 1 for n in out["echoes"]["nodes"].values())
