"""Freezes tests/probes/common_mode_web.py: the common mode as a node of the engine's web, solved by its kernel."""

import pytest


def _silent(*args, **kwargs):
    pass


def test_the_kernel_returns_the_common_mode_and_an_anchor_free_coupling():
    """The node has no base and two depth four edges, the powers of the
    lepton circulant's relative channels.  The kernel's fixed point is
    then 3 Q_l - 1 = 1 - 7.50 (alpha/2pi)^2, the engine's node.  The
    exchange between two electrons is 1.75171e-45, independent of the
    anchor and 54.6 ppm below the measured value.  One channel alone
    gives 1/2 on the searchlight's circulant, and a base of 0.1 gives
    1.1.  The formed lepton closure's second moment is 1 at depths one
    and two, and 1 - 1.0e-5 at depth four, the engine's count.  The
    probe then removes its nodes and restores the engine's common_mode,
    and the web is unchanged."""
    from probes import common_mode_web
    from interference import root

    names = set(root.WEB)
    r = common_mode_web.run(report=_silent)
    assert set(root.WEB) >= names and set(common_mode_web.NODES) & set(root.WEB) == {"common_mode"}
    assert root.WEB["common_mode"].web_depth == 4 and root.WEB.state["common_mode"] == root.COMMON_MODE
    assert r["common_mode"] == pytest.approx(root.COMMON_MODE, abs=1e-15) and r["vents"] == pytest.approx(1.0, abs=1e-12)
    assert r["depths"]["common_mode"] == (4, 4) and r["depths"]["common_mode(5)"] == (5, 5)
    assert r["alpha_G"] == pytest.approx(1.75171e-45, rel=1e-5)
    assert r["ppm"] == pytest.approx(-54.6, abs=0.1) and abs(r["anchor_moves_alpha"]) < 1e-12
    assert r["one_channel"] == pytest.approx(0.5, abs=1e-12) and r["with_base"] == pytest.approx(1.1, abs=1e-12)
    assert r["moment_by_depth"][1] == pytest.approx(1.0, abs=1e-12)
    assert r["moment_by_depth"][2] == pytest.approx(1.0, abs=1e-12)
    assert r["moment_by_depth"][4] - 1 == pytest.approx(-1.01e-5, abs=2e-7)
    assert r["g_ppm"] == pytest.approx(-54.6, abs=0.1) and r["g_ppm_formed"] == pytest.approx(-54.6, abs=0.1)
    assert r["formed"] == pytest.approx(root.COMMON_MODE, abs=1e-15) and abs(r["apart_ppm"]) < 1e-6
