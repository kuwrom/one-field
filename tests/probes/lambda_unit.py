"""Route B of notes/lambda_plan.md: is V0 = -rho_vac a tuning or a unit?

THE QUESTION. In the full form the metric is defined from the field,
g = h/U, and the regulator is a proper time measured in that metric.
The flat vacuum needs V0 = -rho_vac (composite_pole.py K2). If the
metric's own unit could absorb V0, flatness would be automatic and
Lambda would need no mechanism. If not, V0 is a constant of the action
that Lambda measures.

PREMISES. The carrier's real components, free at one loop, with the
proper time regulator of composite_pole.py, t0 and m fixed in physical
units (the regulator is physical). Leading order in the number of
components, as in composite_pole.py.

STATEMENTS (the numbers below check each one).

  B1 (the coordinate scale drops out). On the vacuum metric g = c delta,
     with the regulator measured in g, the gradient average is
     <d_mu phi d_nu phi> = (I1/4) g_mu nu, I1 the physical value, for
     every c (in coordinate units I1 becomes c^2 I1).  The gap equation
     g = <h>/V0 then gives V0 = I1/4 = -rho_vac, with no c left.  The
     unit of the metric is a coordinate choice and cannot absorb V0.
  B2 (the consequence of B1, not separate evidence). At fixed m^2 t0
     and content, Lambda G = 8 pi G^2 (V0 + rho_vac) is linear in V0 and
     vanishes at one point.  V0 is a constant of the action, not a unit,
     so the full form needs a tuning.
  B3 (K1 off the point is a cosmological term). The induced quadratic
     action at q = 0 plus V0 sqrt g is (V0 + rho_vac) times the
     second order part of sqrt det(1 + h). The first order part
     (V0 + rho_vac) tr h / 2 is a tadpole.  Off the point, flat space is
     not a solution: the vacuum is de Sitter or anti de Sitter with
     R = 4 Lambda.  The form h.h - (tr h)^2 / 2 is the cosmological term
     expanded about a nonsolution. It is not a graviton mass, and it has
     no ghost.  (A nonzero clock and ruler profile is a different object:
     Carone arXiv:1911.13166 Sec. III.)
  B4 (Weinberg). The one loop effective action of the auxiliary metric
     is a local functional, int sqrt g (V0 + rho_vac - R/(16 pi G) + ...),
     so the no go theorem applies to it as to any metric. The composite
     metric does not evade it.  An adjustment needs a variable of the
     state that the action does not fix.  In this theory the only
     candidates are the regulator's grain and the carrier mass, and both
     are fixed here.  Treating the grain as a variable of the state is
     route C.

RESULT. Route B closes: V0 = -rho_vac is a tuning, not a unit. One
byproduct: off the point the vacuum is de Sitter, with no graviton mass
and no ghost.

Usage: python3 tests/probes/lambda_unit.py     (a second)
"""

import math
import os
import sys

import numpy as np
from scipy.special import expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import composite_pole as cp  # noqa: E402

N = cp.N_COMPONENTS
M2, T0 = 1.0, 1e-2                     # carrier mass and regulator, physical units


def i1_coordinate(c, m2=M2, t0=T0):
    """I1 = 2 C16 J(2,0,0) of the rescaled field on g = c delta, in coordinate units.

    On g = c delta the operator is -d^2/c + m^2. With psi = sqrt(c) phi the
    field has coordinate mass c m^2 and coordinate proper time cutoff t0/c.
    <d phi d phi> = (1/c) <d psi d psi>.
    """
    J, _, _ = cp.proper_time(c * m2, t0 / c)
    return 2 * cp.C16 * float(J(2, 0, 0.0)) / c


def b1():
    J, rho, _ = cp.proper_time(M2, T0)
    i1_phys = 2 * cp.C16 * float(J(2, 0, 0.0))
    rows = {}
    for c in (0.1, 1.0, 7.0, 300.0):
        # <d_mu phi d_nu phi> = (I1_coord/4) delta should equal (I1_phys/4) g = (I1_phys/4) c delta
        rows[c] = i1_coordinate(c) / (i1_phys * c)
    v0_gap = i1_phys / 4                  # g = <h>/V0 with <h> = (I1/4) g
    return dict(ratio=rows, v0_gap=v0_gap, minus_rho=-rho)


def b2():
    _, rho, stiff = cp.proper_time(M2, T0)
    G = 1.0 / (16 * math.pi * N * stiff)          # composite_pole K5
    out = {}
    for f in (-1.0, -0.999, -1.001, -0.5, 0.0):
        V0 = f * N * rho                            # rho is per real component; f = -1 is flat
        out[f] = 8 * math.pi * G * G * (V0 + N * rho)
    return out


def b3(seed=5):
    rng = np.random.default_rng(seed)
    J, rho, _ = cp.proper_time(M2, T0)
    worst2 = 0.0
    V0 = 0.37 * rho                                 # off the flat point
    for _ in range(5):
        h = rng.normal(size=(4, 4))
        h = 0.5 * (h + h.T)
        # sqrt det(1 + h) = 1 + tr h / 2 + [(tr h)^2 / 8 - tr h^2 / 4] + O(h^3)
        second = np.trace(h)**2 / 8 - np.trace(h @ h) / 4
        e = 1e-3
        dets = [math.sqrt(np.linalg.det(np.eye(4) + s * e * h)) for s in (-1, 0, 1)]
        assert abs((dets[2] - 2 * dets[1] + dets[0]) / (2 * e * e) / second - 1) < 1e-4
        assert abs((dets[2] - dets[0]) / (2 * e) / (0.5 * np.trace(h)) - 1) < 1e-4
        total2 = V0 * second + cp.gamma2(h, np.zeros(4), J)
        worst2 = max(worst2, abs(total2 / ((V0 + rho) * second) - 1))
    return dict(second_order=worst2, tadpole_coefficient=(V0 + rho) / 2)


def run(report=print):
    out = dict(b1=b1(), b2=b2(), b3=b3())
    r1 = out["b1"]
    for c, ratio in r1["ratio"].items():
        assert abs(ratio - 1) < 1e-12, c
    assert abs(r1["v0_gap"] / r1["minus_rho"] - 1) < 1e-12
    report("B1 the coordinate scale drops out: <d phi d phi> / ((I1/4) g) = "
           + ", ".join(f"{v:.12f} (c = {c:g})" for c, v in r1["ratio"].items())
           + f". The gap fixes V0 = I1/4 = -rho_vac = {r1['v0_gap']:.6e}, with c gone.")
    lg = out["b2"]
    assert lg[-1.0] == 0.0 and lg[-0.999] * lg[-1.001] < 0
    assert abs(lg[0.0] / lg[-0.5] - 2.0) < 1e-12
    report("B2 Lambda G at V0 = f (-rho_vac), fixed m^2 t0: "
           + ", ".join(f"f = {-f:g}: {v:+.3e}" for f, v in lg.items())
           + ".  Linear in V0, zero at one point: a tuning, not a unit.")
    r3 = out["b3"]
    assert r3["second_order"] < 1e-6
    report(f"B3 off the point, induced q = 0 term + V0 sqrt g = (V0 + rho_vac) x second order of "
           f"sqrt det(1 + h) to {r3['second_order']:.1e}, with the tadpole "
           f"(V0 + rho_vac) tr h / 2: a cosmological term about a nonsolution, no graviton mass.")
    report("RESULT: route B closes. V0 = -rho_vac is a tuning. Off it the vacuum is de Sitter, "
           "with no mass and no ghost. An adjustment needs a variable of the state: route C.")
    return out


if __name__ == "__main__":
    run()
