"""EM channel finality as a tested exclusion: what a further edge on 1/alpha would have to be.

THE QUESTION. root.py reads 1/alpha(0) as the real root of
x^3 = (512/pi)[(1 - 1/(2 pi)) x^2 - 1/(2 pi^2)]. The two terms are the
bridge's depth one edge and the depth-4 reading of the closed carrier
with orientation count 2, and "no separate (alpha/2pi)^3 edge exists"
(EM CHANNEL FINALITY). The reading matches the Cs recoil value and
sits -7.8 and -14 standard deviations from the electron moment and the
Rb recoil. Is finality a computed exclusion? Which further edges does
the web's own structure allow, and what would an edge that rescues the
other two values have to be?

PREMISES. P1 the sector node rule: the EM channel couples the two
charged sectors, the G2-knot sector and the F4-knot sector, by one edge
(root.py). P2 the web's one EM self-loop is the lepton factor
lepton_EM = 1 - alpha/2pi at depth two (root.py). P3 the carrier rule:
every reader of the carrier carries one sign on its ledger, negative
on 1/alpha (tests/audit/carrier_ledger.py C3). P4 the three
determinations of 1/alpha(0) and their uncertainties, as tests.

STATEMENTS (the numbers below check each one).

  E1 (parity). The sector graph of P1 is bipartite, so its closed
     walks have even length. Every traversal of the carrier is then a
     power of L = (alpha/2pi)^2, and the kernel's fixed point resums
     them. An odd
     power, (alpha/2pi)^3, needs a self-loop, and the web's one EM
     self-loop is the lepton factor of P2.
  E2 (the self-loop composed). Inserted into the carrier the way
     1/alpha reads it, -2 L (1 - alpha/2pi)^n, the self-loop is a cubic
     term +2n in units of (512/pi)(alpha/2pi)^3. It moves 1/alpha(0) by
     +5.1e-7 at n = 1 and +1.0e-6 at n = 2, toward the Rb side and past
     it. Every measurement, the Rb recoil included, excludes both by 18
     to 79 standard deviations. A double traversal of the carrier, -4 L^2
     at the largest orientation weight, is below a tenth of a standard
     deviation of Cs, and the kernel resums it.
  E3 (what a rescue needs). In those units the measurements ask for a
     cubic coefficient of -0.01 +- 0.11 (Cs), +0.46 +- 0.06 (electron
     moment), +0.50 +- 0.08 (CODATA) and +0.61 +- 0.04 (Rb). The one odd
     source supplies 2n, an orientation count times a whole number of
     self-loops. 0 matches Cs, and 2 overshoots every other value by
     more than three times. A coefficient near 1/2 needs a vertex weight
     on the self-loop that the web does not register.

  E4 (the one match, and the rule it breaks). With a conformal weight
     w of the web on the self-loop, -2 L (1 - w alpha/2pi), the cubic
     coefficient is 2w. Of the seven weights {2/9, 2/5, 1/2, 5/9, 3/5,
     8/9, 1}, one matches: w = h10 = 2/9 gives 4/9, with the electron
     moment at -0.2 and CODATA at -0.7 standard deviations, Cs at +4.3
     and Rb at -3.9. The vertex composition rule gives the coupling
     node no matter vertex. Its weight is the orientation count alone
     (root.py). But h10 is a matter vertex, the lepton's phase, so this
     match needs the rule broken. A resolution of the recoil dispute
     toward the electron moment would therefore falsify the vertex
     composition rule, not move a coefficient.

RESULT. Finality is an exclusion with failing cases. The bipartite
sector graph leaves the lepton self-loop as the one source of an odd
term. Its coefficient is a whole multiple of two, and composing it
fails every measurement. The reading commits to the Cs recoil value.
A resolution of the recoil dispute on the Rb side would falsify it,
and the web has no further edge to supply.

Scope: the EM channel's sector graph and the web's registered vertices.
The three determinations all measure the on shell Thomson coupling,
which is the engine's coordinate.

Usage: python3 tests/probes/em_finality.py
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import root as r  # noqa: E402

ALPHA = {"Cs recoil": (137.035999046, 27e-9), "electron moment": (137.035999166, 15e-9),
         "CODATA 2022": (137.035999177, 21e-9), "Rb recoil": (137.035999206, 11e-9)}


def solve(n=0, twice=0.0, cubic=0.0):
    """1/alpha(0) with the carrier read as -2 L (1 - alpha/2pi)^n, a double traversal -twice L^2, and a
    cubic (alpha/2pi)^3, solved self-consistently."""
    inv0, x = r._inv0, r.inv_alpha_phys
    for _ in range(200):
        u = 1 / (2 * math.pi * x)
        x = inv0 * (1 - 1 / (2 * math.pi) - 2 * u**2 * (1 - u)**n - twice * u**4 + cubic * u**3)
    return x


def pulls(x):
    return {k: (x - v) / s for k, (v, s) in ALPHA.items()}


def run(report=print):
    out = {}

    # E1: parity of the sector graph.
    A = np.array([[0, 1], [1, 0]])
    loop = np.array([[1, 1], [1, 0]])
    odd = [int(np.trace(np.linalg.matrix_power(A, k))) for k in (1, 3, 5)]
    odd_loop = [int(np.trace(np.linalg.matrix_power(loop, k))) for k in (1, 3, 5)]
    assert odd == [0, 0, 0] and all(v > 0 for v in odd_loop)
    assert abs(solve() - r.inv_alpha_phys) < 1e-9
    out.update(odd=odd, odd_loop=odd_loop)
    report(f"E1 closed walks of odd length on the sector graph: {odd}. With the lepton self-loop {odd_loop}")

    # E2: the self-loop composed, and repeated traversals.
    base = solve()
    readings = {"canon": base, "self-loop n = 1": solve(n=1), "self-loop n = 2": solve(n=2),
                "carrier twice": solve(twice=4.0)}
    table = {k: (x, x - base, pulls(x)) for k, x in readings.items()}
    for k in ("self-loop n = 1", "self-loop n = 2"):
        assert min(abs(p) for p in table[k][2].values()) > 9
    assert abs(table["carrier twice"][1]) / ALPHA["Cs recoil"][1] < 0.1
    out.update(table=table)
    for k, (x, d, p) in table.items():
        report(f"E2 {k:<16s} 1/alpha = {x:.10f} ({d:+.2e}). " + ", ".join(f"{m} {v:+.1f}" for m, v in p.items()))

    # E3: the cubic coefficient each measurement asks for.
    u = 1 / (2 * math.pi * base)
    unit = r._inv0 * u**3
    ask = {k: (v - base) / unit for k, (v, _) in ALPHA.items()}
    ask_sigma = {k: s / unit for k, (_, s) in ALPHA.items()}
    supplied = {n: (solve(n=n) - base) / unit for n in (1, 2)}
    assert abs(supplied[1] - 2) < 0.01 and abs(supplied[2] - 4) < 0.02
    assert abs(ask["Cs recoil"]) < ask_sigma["Cs recoil"]
    assert all(0 < ask[k] < 2 / 3 for k in ("electron moment", "CODATA 2022", "Rb recoil"))
    for k in ALPHA:                                   # the asked coefficient reproduces each value
        assert abs(solve(cubic=ask[k]) - ALPHA[k][0]) < 0.05 * ALPHA[k][1]
    out.update(ask=ask, ask_sigma=ask_sigma, supplied=supplied)
    report("E3 cubic coefficient asked, in (512/pi)(alpha/2pi)^3: "
           + ", ".join(f"{k} {ask[k]:+.3f} +- {ask_sigma[k]:.3f}" for k in ALPHA))

    report(f"   the self-loop supplies {supplied[1]:.3f} (n = 1) and {supplied[2]:.3f} (n = 2)")

    # E4: a conformal weight on the self-loop.
    weights = {"2/9": 2 / 9, "2/5": 0.4, "1/2": 0.5, "5/9": 5 / 9, "3/5": 0.6, "8/9": 8 / 9, "1": 1.0}
    landing = {w: {k: (2 * v - ask[k]) / ask_sigma[k] for k in ALPHA} for w, v in weights.items()}
    lands = [w for w, p in landing.items() if min(abs(x) for x in p.values()) < 1]
    assert lands == ["2/9"]
    assert abs(landing["2/9"]["electron moment"] + 0.2) < 0.05 and abs(landing["2/9"]["Rb recoil"] + 3.9) < 0.1
    out.update(landing=landing, lands=lands)
    report("E4 a weight on the self-loop matches for " + ", ".join(lands) + ": "
           + ", ".join(f"{k} {v:+.1f}" for k, v in landing["2/9"].items()))
    report("RESULT: the sector graph is bipartite, the one odd source is the lepton self-loop, its coefficient is a "
           "whole multiple of two, and it fails every measurement.")
    return out


if __name__ == "__main__":
    run()
