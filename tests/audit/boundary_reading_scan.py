"""
boundary_reading_scan.py: what a term in the boundary reading would need from G to be tested.

THE QUESTION.  The anchor's G has one place where its residual can
sit: the boundary reading c in G = 1/(c M_b^2) (anchor_account.py A4).
G is the only observable that reads c.  The engine fixes c through the
common mode, c = 1/(3 Q_l - 1), which is the count of the lepton
circulant on the depth four web (root.py, THE COMMON MODE;
tests/probes/common_mode_count.py).  Write a further term as
G = G_engine/(1 + eps), so that eps = 0 is the engine's reading.
Three questions follow.  Which eps does each reading of the
measurements ask for?  How many terms of the engine's own echo grammar
fall inside that window by construction?  How precise would G have to
be before such a match meant anything?

PREMISES.  P1 the engine's G, 6.673936e-11 (root.py G_PRED).  P2
three readings of the sixteen inputs (g_measurements.py): the
recommended value,
6.67430(15). An added spread, 6.67389(26). And the five most precise
inputs alone, combined with their Birge ratio.  P3 the menu, fixed
before any comparison: eps = +-(p/q) x^n with p, q from 1 to 32
coprime, n = 1 or 2, and x one of the engine's three loop parameters,
alpha/2pi, alpha_G2(M_b)/2pi = 1/(48 pi^2) and alpha_G2(mu*)/2pi =
1/64.  Every depth-2 and depth-4 echo of the engine has this form
(root.py, THE ECHO LAW).  Its weights are the small rationals of the
engine's dictionary.

STATEMENTS (the numbers below check each one).

  B1 (the target depends on the reading).  The recommended value asks
     for eps = -54 ppm +- 22.  The added spread asks for +7 +- 39,
     which includes zero.  The five most precise inputs ask for
     -56 +- 9.  So a term chosen to close G is first a choice among
     readings of the data, and only then a statement about the
     boundary.
  B2 (the grammar matches everywhere).  The menu holds 7,710 distinct
     values, 2,978 of them between 1 and 300 ppm.  Within one standard
     deviation of the target, 121 values fall for the recommended
     value, 39 for the five most precise inputs and 2,615 for the
     added spread.  Windows of the same widths placed at -30 and
     -60 ppm hold 113 to 302 and 41 to 73 values.  A match inside such
     a window therefore measures the density of the menu, not the
     physics.  The added spread's window contains zero, where the menu
     is densest.
  B3 (what G would have to reach).  Between 30 and 60 ppm
     neighbouring menu values lie 0.28 ppm apart at the median.  To
     tell a term in c from its neighbours, G would have to be known to
     about 0.14 ppm.  That is 162 times below the recommended value's
     uncertainty and 65 times below that of the five most precise
     inputs taken together.
  B4 (the engine's reading).  eps = 0 is the engine's reading: the
     count, with no weight to choose.  It asks no choice of the
     measurements.  It sits at +0.2 standard deviations from the added
     spread, -2.4 from the recommended value and -6.2 from the five
     most precise inputs.  Its depth four term,
     3 Q_l - 1 = 1 - 7.50 (alpha/2pi)^2, is -10.1 ppm, or 0.45 of the
     recommended value's standard deviation.  The lepton closure fixes
     this term, not the menu.  The menu's nearest value,
     -(15/2) (alpha/2pi)^2, lies 0.003 ppm from it, which is a match of
     the kind B2 counts.  G does not resolve it.

RESULT.  No term in the boundary reading can be tested by G at any
precision G will reach soon.  The engine's own echo grammar puts a
candidate in every window G allows, and the windows themselves move
with the rule that combines the measurements.  The engine's reading,
the count with no weight to choose, is the only reading that needs no
choice.  The five most precise inputs taken alone exclude it.

Usage: python3 tests/audit/boundary_reading_scan.py     (a few seconds)
"""

import math
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "probes"))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402
import g_measurements as gm  # noqa: E402

G_ANCHOR = r.G_PRED * 1e11                  # P1, the engine's G, in 1e-11 m^3 kg^-1 s^-2
LOOPS = {"alpha/2pi": 1.0 / (2 * math.pi * r.inv_alpha_phys),
         "alpha_G2(M_b)/2pi": float(r.alpha_G2_Pl) / (2 * math.pi),
         "alpha_G2(mu*)/2pi": float(r.alpha_G2_WZW) / (2 * math.pi)}
P_MAX = 32


def menu():
    """Distinct eps values of P3, in ppm, with one label each."""
    out = {}
    for name, x in LOOPS.items():
        for n in (1, 2):
            for p in range(1, P_MAX + 1):
                for q in range(1, P_MAX + 1):
                    if math.gcd(p, q) != 1:
                        continue
                    for sign in (1, -1):
                        v = round(sign * 1e6 * p / q * x**n, 9)
                        out.setdefault(v, f"{'-' if sign < 0 else '+'}({Fraction(p, q)}) ({name})^{n}")
    return out


def readings():
    """P2: name -> (G, standard uncertainty), and the eps each asks for, in ppm."""
    rules = gm.combinations()
    precise = [row for row in gm.TABLE if row[2] / row[1] < gm.PRECISE]
    x = np.array([row[1] for row in precise])
    s = np.array([row[2] for row in precise])
    w = 1.0 / s**2
    mean = float((w * x).sum() / w.sum())
    birge = max(1.0, math.sqrt(float((w * (x - mean)**2).sum()) / (len(x) - 1)))
    data = {"recommended": gm.RECOMMENDED, "added spread": rules["added spread, chi-square = n - 1"],
            "five most precise": (mean, birge * float(w.sum())**-0.5)}
    eps = {k: (1e6 * (G_ANCHOR / v - 1.0), 1e6 * G_ANCHOR * u / v**2) for k, (v, u) in data.items()}
    return data, eps


def landings(values, centre, width):
    return int(np.sum(np.abs(values - centre) <= width))


def run(report=print):
    out = {}
    data, eps = readings()

    # B1: the target depends on the reading.
    assert abs(eps["recommended"][0] + 54.4) < 0.3 and abs(eps["recommended"][1] - 22.5) < 0.3
    assert abs(eps["added spread"][0]) < eps["added spread"][1]
    assert abs(eps["five most precise"][0] + 56.0) < 0.3 and abs(eps["five most precise"][1] - 9.0) < 0.1
    out.update(targets=eps)
    report("B1 targets (ppm): " + "; ".join(f"{k} {c:+.1f} +- {u:.1f}" for k, (c, u) in eps.items()))

    # B2: the grammar matches everywhere.
    m = menu()
    values = np.array(sorted(m))
    in_range = int(np.sum((np.abs(values) >= 1) & (np.abs(values) <= 300)))
    hits, controls = {}, {}
    for k, (c, u) in eps.items():
        hits[k] = landings(values, c, u)
        controls[k] = [landings(values, c0, u) for c0 in (-30.0, -60.0)]
    assert (len(m), in_range) == (7710, 2978) and (hits["recommended"], hits["five most precise"]) == (121, 39)
    assert all(min(controls[k]) > 0.3 * hits[k] for k in ("recommended", "five most precise"))
    nearest = {k: m[min(m, key=lambda v: abs(v - c))] for k, (c, _) in eps.items() if k != "added spread"}
    out.update(menu_size=len(m), in_range=in_range, hits=hits, controls=controls, nearest=nearest)
    report(f"B2 menu: {len(m)} distinct values, {in_range} between 1 and 300 ppm. Values within one standard "
           "deviation: " + "; ".join(f"{k} {hits[k]} (controls at -30, -60 ppm: {controls[k][0]}, {controls[k][1]})"
                                     for k in eps)
           + ". Nearest to the recommended target " + nearest["recommended"])

    # B3: what G would have to reach.
    near = np.sort(values[(values > -60) & (values < -30)])
    gap = float(np.median(np.diff(near)))
    needed = gap / 2                                    # a term must be told from its neighbours at about one sigma
    ratio = {k: u / needed for k, (_, u) in eps.items()}
    assert abs(gap - 0.278) < 0.002 and 150 < ratio["recommended"] < 170 and 60 < ratio["five most precise"] < 70
    out.update(median_gap_ppm=gap, needed_ppm=needed, ratio=ratio)
    report(f"B3 resolution: neighbouring menu values near 45 ppm lie {gap:.3f} ppm apart at the median. Telling a "
           f"term from its neighbours needs G to about {needed:.3f} ppm, " + ", ".join(
               f"{ratio[k]:.0f} times below the {k} uncertainty" for k in eps))

    # B4: the engine's reading, and its depth four term beside the menu.
    pulls = {k: (G_ANCHOR - v) / u for k, (v, u) in data.items()}
    carrier = 1e6 * (r.COMMON_MODE - 1.0)                              # the count, 3 Q_l - 1, minus one
    closest = min(m, key=lambda v: abs(v - carrier))
    assert abs(pulls["recommended"] + 2.43) < 0.01 and abs(pulls["added spread"] - 0.17) < 0.01
    assert abs(pulls["five most precise"] + 6.22) < 0.01
    assert abs(carrier / eps["recommended"][1] + 0.45) < 0.01
    assert m[closest] == "-(15/2) (alpha/2pi)^2" and abs(closest - carrier) < 0.004
    out.update(engine_pulls=pulls, carrier_ppm=carrier, menu_nearest=(m[closest], closest - carrier))
    report("B4 the engine's reading: " + "; ".join(f"{k} {p:+.2f} standard deviations" for k, p in pulls.items())
           + f". Its depth four term {carrier:+.2f} ppm, the menu's nearest {m[closest]} {closest - carrier:+.4f} ppm away")

    report("RESULT: the echo grammar puts a candidate term in every window G allows, and the windows move with the "
           "rule that combines the measurements. The engine's reading, a count, needs no choice, and the five most precise "
           "inputs alone exclude it.")
    return out


if __name__ == "__main__":
    run()
