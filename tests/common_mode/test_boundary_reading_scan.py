"""Freezes tests/audit/boundary_reading_scan.py: no term in the boundary
reading can be tested by G at the precision G has."""


def _silent(*args, **kwargs):
    pass


def test_boundary_reading_scan_grammar_lands_in_every_window_g_allows():
    """The three readings of the measurements each ask for a term
    beyond the engine's reading: eps = -54 +- 22, +7 +- 39 and
    -56 +- 9 ppm.  The
    engine's echo grammar puts 121 and 39 values in the two narrow
    windows, and similar counts in control windows.  To tell one term
    from its neighbours, G must be known to about 0.14 ppm. The engine's reading sits at
    -2.43, +0.17 and -6.22 standard deviations, and its depth four term,
    the count's -10.12 ppm, lies 0.003 ppm from the menu's
    -(15/2) (alpha/2pi)^2."""
    import boundary_reading_scan as b

    out = b.run(report=_silent)
    assert out["hits"]["recommended"] == 121 and out["hits"]["five most precise"] == 39
    assert out["controls"]["recommended"] == [302, 113] and out["controls"]["five most precise"] == [73, 41]
    assert abs(out["median_gap_ppm"] - 0.278) < 0.002
    p = out["engine_pulls"]
    assert abs(p["recommended"] + 2.43) < 0.02 and abs(p["added spread"] - 0.17) < 0.02
    assert abs(p["five most precise"] + 6.22) < 0.02
    assert abs(out["carrier_ppm"] + 10.120) < 1e-3
    assert out["menu_nearest"][0] == "-(15/2) (alpha/2pi)^2" and abs(out["menu_nearest"][1]) < 0.004
