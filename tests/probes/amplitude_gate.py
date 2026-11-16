"""The circulant amplitude B/A against the angle menu: a failing case that uses no measured value.

THE QUESTION. The lepton circulant is Delta_k = 1 + B cos(theta + 2 pi k/3)
with B = B/A = sqrt(d10) and theta = h10 = 2/9 (interference/masses.py).
B/A is the channel map's identification: the Fano CG ratio, the
orientation pair and the quantum Schur weight 1/d(fund) (the selection
table, interference/__init__.py). The angle has a finite menu with one
survivor (tests/audit/edge_grammar.py theta_menu). Does that survivor
depend on B, and does the classical value of the map fail the test?

PREMISES. P1 the committed angle menu: every nonzero primary weight of
SU(3)_3 and its multiples pi h and 2 pi h, 14 angles. P2 the blind
gates: every gap positive (sqrt m real) and three distinct masses.
P3 the two traces of the channel map: the quantum spherical trace gives
1/d(fund) = 1/2 and B = sqrt 2. The classical trace gives
B = 1/sqrt 2, Q = 5/12 (searchlight/octonions.py).

STATEMENTS (the numbers below check each one).

  A1 (the window). The menu has exactly one survivor of P2, and that
     survivor is h10 = 2/9, for B between 1.30541 and 1.47367 (Q between
     0.6173 and 0.6953). Below the window a second angle survives,
     above it none.
  A2 (the two traces). The quantum trace, B = sqrt 2, lies inside. The
     classical trace, B = 1/sqrt 2, leaves 12 of the 14 angles on the
     branch, and the angle is not selected.
  A3 (what the window does not do). It fixes B only to within six
     percent. Inside the window, the lepton mass ratios fix B to 5e-10.

RESULT. The angle menu's unique survivor constrains the amplitude. It
exists only on a window of B, and that window contains the quantum
trace and excludes the classical one. No measured value enters. B/A
stays one identification, the triality amplitude rule of the channel
map, and this probe is its failing case that does not use G.

Usage: python3 tests/probes/amplitude_gate.py
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import edge_grammar as eg  # noqa: E402
from interference import root as r  # noqa: E402

MU_E = (206.7682827, 0.0000046)


def survivors(B, menu):
    """Angles of the menu on the positive, nondegenerate branch at amplitude B."""
    out = []
    for th in menu:
        D = [1 + B * math.cos(th + 2 * math.pi * k / 3) for k in range(3)]
        if min(D) <= 0:
            continue
        m = sorted(d * d for d in D)
        if min(m[1] - m[0], m[2] - m[1]) < 1e-9 * m[2]:
            continue
        out.append(th)
    return out


def run(report=print):
    out = {}
    menu = eg.theta_menu()
    h10 = float(r.h10)
    assert len(menu) == 14

    # A1: the window.
    grid = np.linspace(1.2, 1.6, 40001)
    unique = [B for B in grid if survivors(B, menu) == [h10]]
    lo, hi = min(unique), max(unique)
    Q = lambda B: 1 / 3 + B * B / 6
    assert abs(lo - 1.30541) < 2e-5 and abs(hi - 1.47367) < 2e-5
    assert len(survivors(lo - 0.01, menu)) == 2 and survivors(hi + 0.01, menu) == []
    out.update(window=(lo, hi), Q_window=(Q(lo), Q(hi)))
    report(f"A1 unique survivor 2/9 for B in [{lo:.5f}, {hi:.5f}], Q in [{Q(lo):.4f}, {Q(hi):.4f}]")

    # A2: the two traces.
    quantum, classical = math.sqrt(r.d10), 1 / math.sqrt(2)
    n_q, n_c = len(survivors(quantum, menu)), len(survivors(classical, menu))
    assert lo < quantum < hi and abs(quantum - r.BA_ratio) < 1e-15
    assert n_q == 1 and n_c == 12
    out.update(n_quantum=n_q, n_classical=n_c)
    report(f"A2 quantum trace B = {quantum:.5f}: {n_q} survivor. Classical trace B = {classical:.5f}: {n_c} survivors")

    # A3: the window against the lepton ratio.
    gap = lambda B, k: 1 + B * math.cos(h10 + 2 * math.pi * k / 3)
    ratio = lambda B: (gap(B, 2) / gap(B, 1))**2
    slope = (ratio(quantum * (1 + 1e-7)) / ratio(quantum * (1 - 1e-7)) - 1) / 2e-7
    held = MU_E[1] / MU_E[0] / abs(slope)
    width = (hi - lo) / (2 * quantum)
    assert held < 1e-9 and 0.05 < width < 0.07
    out.update(held=held, width=width)
    report(f"A3 the window is +-{100 * width:.1f}% of B; m_mu/m_e fixes B to {held:.1e}")
    return out


if __name__ == "__main__":
    run()
