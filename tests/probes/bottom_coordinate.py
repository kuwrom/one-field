"""The bottom mass: the relation the engine uses and the coordinate it is compared in.

THE QUESTION. The engine computes m_b from the (c, b, t) Koide relation
with Q = Q0 + h11/K^3 = 289/432, using its own m_c and m_t. The result
is 4193.8 MeV, +1.5 standard deviations from the PDG m_b(m_b) =
4183(7). The triplet is compared in a mixed coordinate: m_c(m_c), m_b(m_b), and the
direct top mass. How much of the residual belongs to the relation, and
how much to the inputs and the coordinate?

PREMISES. P1 the relation of interference/masses.py, the negative free
branch bracketed between m_c and m_t. P2 reference values: PDG 2024
m_c(m_c) = 1273.0(4.6), m_t = 172.57(29) GeV direct, m_b(m_b) =
4183(7); FLAG 2024 m_b(m_b) = 4200(14) for Nf = 2+1+1 and 4171(20) for
Nf = 2+1 (arXiv:2411.04268, Eqs. 61, 63).

STATEMENTS (the numbers below check each one).

  B1 (the lever). d ln m_b / d ln m_t = +1.61 and d ln m_b / d ln m_c =
     -0.61. One GeV of top moves m_b by 39 MeV, so the 8 to 11 MeV
     residual equals a quarter of a GeV of top. That is the size of
     the spread between the top's direct and cross section masses.
  B2 (the relation fed with measured inputs). With the measured m_c and
     m_t the relation gives 4197.8 +- 14.6 MeV. With the spread of the
     inputs included, that is +0.9 standard deviations from 4183(7).
     This number tests the relation itself, independent of the
     engine's m_c and m_t.
  B3 (the references). The engine's 4193.8 is +1.5 sigma from PDG,
     -0.4 sigma from FLAG Nf = 2+1+1 and +1.1 sigma from FLAG Nf = 2+1.

RESULT. The residual on m_b is inside the spread of the coordinate
itself. The relation fed with measured inputs is at +0.9 sigma, and the
two lattice averages bracket the engine's value. The mixed coordinate
(MS-bar for c and b, direct for t) is a declared choice, and it moves
m_b by 39 MeV per GeV of top.

Usage: python3 tests/probes/bottom_coordinate.py
"""

import contextlib
import io
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import masses  # noqa: E402
from interference import root as r  # noqa: E402

M_C, M_T = (1273.0, 4.6), (172570.0, 290.0)
REFS = {"PDG m_b(m_b)": (4183.0, 7.0), "FLAG Nf=2+1+1": (4200.0, 14.0), "FLAG Nf=2+1": (4171.0, 20.0)}


def engine():
    """The engine's m_c, m_t and m_b in MeV, from the chain."""
    with contextlib.redirect_stdout(io.StringIO()):
        m = masses.derive(r.derive())
    return {k: m[k] for k in ("m_c", "m_t", "m_b")}


def bottom(m_c, m_t, Q=float(r.Q0 + r.h11 / r.K**3)):
    sc, st = math.sqrt(m_c), math.sqrt(m_t)
    a, b = Q - 1, 2 * Q * (sc + st)
    c = (Q - 1) * (m_c + m_t) + 2 * Q * sc * st
    for sign in (1, -1):
        x = (-b + sign * math.sqrt(b * b - 4 * a * c)) / (2 * a)
        if m_c < x * x < m_t:
            return x * x
    raise ValueError("no bracketed root")


def run(report=print):
    out = {}
    ENGINE = engine()
    assert abs(ENGINE["m_b"] - 4193.8) < 0.05
    assert abs(bottom(ENGINE["m_c"], ENGINE["m_t"]) - ENGINE["m_b"]) < 1e-6

    # B1: the lever.
    eps = 1e-6
    b0 = bottom(ENGINE["m_c"], ENGINE["m_t"])
    dt = (bottom(ENGINE["m_c"], ENGINE["m_t"] * (1 + eps)) / b0 - 1) / eps
    dc = (bottom(ENGINE["m_c"] * (1 + eps), ENGINE["m_t"]) / b0 - 1) / eps
    per_gev = dt * b0 / ENGINE["m_t"] * 1e3
    assert abs(dt - 1.608) < 0.005 and abs(dc + 0.608) < 0.005 and abs(per_gev - 39.1) < 0.3
    out.update(dlnt=dt, dlnc=dc, per_gev=per_gev)
    report(f"B1 d ln m_b/d ln m_t = {dt:+.3f}, d ln m_b/d ln m_c = {dc:+.3f}: {per_gev:.1f} MeV per GeV of top")

    # B2: the relation fed with measured inputs.
    fed = bottom(M_C[0], M_T[0])
    spread = math.hypot(bottom(M_C[0] + M_C[1], M_T[0]) - fed, bottom(M_C[0], M_T[0] + M_T[1]) - fed)
    pull = (fed - REFS["PDG m_b(m_b)"][0]) / math.hypot(REFS["PDG m_b(m_b)"][1], spread)
    assert abs(fed - 4197.8) < 0.1 and abs(spread - 14.6) < 0.1 and abs(pull - 0.91) < 0.01
    out.update(fed=fed, spread=spread, fed_pull=pull)
    report(f"B2 the relation on measured m_c, m_t: {fed:.1f} +- {spread:.1f} MeV, {pull:+.2f} sigma of PDG")

    # B3: the references.
    pulls = {k: (ENGINE["m_b"] - v) / s for k, (v, s) in REFS.items()}
    assert abs(pulls["PDG m_b(m_b)"] - 1.54) < 0.01 and abs(pulls["FLAG Nf=2+1+1"] + 0.44) < 0.01
    out.update(pulls=pulls)
    report("B3 the engine's 4193.8 MeV: " + ", ".join(f"{k} {p:+.2f}" for k, p in pulls.items()))
    return out


if __name__ == "__main__":
    run()
