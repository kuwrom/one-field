"""Two lumps and the induced kernel: the force between them.

THE QUESTION. lump_in_empty_field.py L3 reads a lump's weight from the
flat space integral of its stress, with the metric's response as a
premise. Here the response is the kernel: two lumps of different
binding exchange the one list's induced kernel (one_list_kernel.py),
each sourcing it with its own energy momentum tensor, the tensor
T = 2 dS/dh through which a row couples (saddle_stabilizer.py Z8). What is
the potential between them, and what does each lump weigh in it?

PREMISES. P1 the lump of lump_in_empty_field.py, Phi = sigma(r)
exp(-i Omega t), m = c = 1, f = sqrt 2 sigma, potential
V = a f^4 + kappa f^6 / 3. P2 its energy momentum tensor, static and
spherical:
    T_00 = (Omega^2 + 1) f^2/2 + f'^2/2 + V,
    T_rr = mu f^2 + f'^2/2 - V,   T_tt = mu f^2 - f'^2/2 - V (tangential).
P3 the one list's kernel at its anchored regulator, K = F2 x (spin 2)
+ F0 x (spin 0), and the exchange between conserved sources
    W(q) = -(1/8) [S2(T1, T2)/F2 + S0(T1, T2)/F0]
(composite_pole.py K5), with the sources' time components continued to
the kernel's Euclidean signature. P4 the carrier's mass in units of the
regulator, eps = m sqrt(t0): the lump's size against the kernel's.

STATEMENTS (the numbers below check each one).

  W1 (the sources). The profile's own fields give int T_00 = E, the
     energy of lump_in_empty_field.py L2 with binding and clock energy
     included. The stress is conserved, T_rr' + 2 (T_rr - T_tt)/r = 0,
     so its transform is transverse, T_ij(q) = (delta_ij - q_i q_j/q^2)
     s(q), and s(0) = int T_kk / 2 = 0 (von Laue).
  W2 (the potential). V(r) = (1/2 pi^2) int q^2 dq j0(q r) W(q). For
     two lumps whose binding energies are +0.25 and -25 percent of
     their rest mass, V r / (G E1 E2) = -1 at separations of 60, 120
     and 240 carrier lengths: Newton's law, attractive, with each
     lump's total energy as its weight. In units of the charges' rest
     mass, m Q for each, the same potential is -0.752: the binding
     energy gravitates.
  W3 (the stress does not weigh). The cross terms of energy density
     and stress, the pressure's share of the active mass, change the
     potential by less than 1e-6 at those separations: a static lump's
     stress has no monopole.
  W4 (the kernel is read). With the spin-0 form factor moved by ten
     percent, F0/F2 = -2.2, the same sources give -1.0303, the value
     2 (2/3 - 1/(6 x 1.1)): the coefficient is the kernel's, and the
     -1 of W2 is its F0/F2 = -2.
  W5 (the regulator's reach). With the kernel's full form factors, at
     eps = 0.01 and 0.1, the potential at those separations moves by
     less than 1e-4. The regulator changes the force only where the
     lumps overlap the regulator scale. F2 and F0 have no zero on the
     lumps' momenta.

RESULT. Two lumps attract with the potential -G E1 E2 / r, E the
lump's whole energy, binding included, and the coefficient is the
induced kernel's own. No value of G is read here: the potential is
quoted in units of G = 1/(16 pi x stiffness).

Scope: static sources at linear order in the kernel, spherical lumps on
the clock orientation, separations beyond the lumps' overlap.

Usage: python3 tests/probes/two_lumps.py     (about ten seconds)
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r_  # noqa: E402
import lump_in_empty_field as lf  # noqa: E402
import one_list_kernel as olk  # noqa: E402

FRACS = (0.3, 0.85)                          # mu / mu0 of the two lumps
SEPARATIONS = (60.0, 120.0, 240.0)           # carrier lengths 1/m
EPS = (0.01, 0.1)                            # m sqrt(t0)
Q = np.arange(0.0, 14.0, 0.002)              # momentum grid, carrier units
DR = 0.004


def source(frac):
    """A lump's energy density and stress on a radial grid, and their transforms on Q."""
    d = lf.lump(frac * lf.MU0)
    rel = lf.lift(d)
    mu, k = d["mu"], math.sqrt(-2 * d["mu"])
    r = np.arange(DR / 2, d["rc"] + 45.0 / k, DR)
    f = d["profile"](r)
    fp = np.gradient(f, r, edge_order=2)
    V = lf.A * f**4 + lf.KAPPA * f**6 / 3
    rho = (rel["omega"]**2 + 1) * f * f / 2 + fp * fp / 2 + V
    Trr = mu * f * f + fp * fp / 2 - V
    Ttt = mu * f * f - fp * fp / 2 - V
    w = 4 * math.pi * r * r * DR
    div = np.gradient(Trr, r, edge_order=2) + 2 * (Trr - Ttt) / r
    scale = np.abs(np.gradient(Trr, r, edge_order=2)).max()

    def transform(g):
        out = np.empty_like(Q)
        for i in range(0, Q.size, 500):                   # j0(q r), in blocks
            out[i:i + 500] = np.sinc(np.outer(Q[i:i + 500], r) / math.pi) @ (w * g)
        return out
    return dict(E=rel["E"], mQ=rel["omega"] * d["N"], binding=rel["binding"], E_fields=float(w @ rho),
                laue=float(w @ (Trr + 2 * Ttt)) / float(w @ np.abs(Trr + 2 * Ttt)),
                conservation=float(np.abs(div[5:-5]).max() / scale),
                rho=transform(rho), s=0.5 * transform(Trr + 2 * Ttt), radius=d["rc"])


def kernel_ratios():
    """The one list's F2 and F0 over their Einstein values, as functions of y = q^2 t0."""
    table = B.rows(r_.d10, r_.d11, r_.n7, r_.n26)
    t0 = B.anchored(table, 1.0, float(r_.xi_bridge))[0] / (2 * math.pi)
    rows = olk.kernel_rows(table)
    sig = olk.stiffness(rows, t0)[0]
    ly = np.linspace(math.log(1e-7), math.log(30.0), 70)
    e2, e0 = [], []
    for y in np.exp(ly):
        F2, F0 = olk.form_factors(rows, t0, y)[:2]
        ein = sig * olk.UNIT * y / (4 * t0 * t0)
        e2.append(F2 / ein)
        e0.append(F0 / (-2 * ein))

    def at(y):
        ly_ = np.log(np.maximum(y, 1e-7))
        return np.interp(ly_, ly, e2), np.interp(ly_, ly, e0)
    return at, (min(e2), min(e0))


def potential(s1, s2, sep, ratios=None, eps=0.0, f0_shift=1.0, stress=True, weights=None):
    """V r / (G w1 w2): the exchange of P3 summed over q, in units of the weights' Newton value."""
    rho1, rho2 = s1["rho"], s2["rho"]
    p1, p2 = (s1["s"], s2["s"]) if stress else (0.0 * rho1, 0.0 * rho2)
    e2, e0 = ratios(eps * eps * Q * Q) if ratios and eps else (1.0, 1.0)
    A = rho1 * rho2 + 2 * p1 * p2                         # tr(theta T1 theta T2)
    b1, b2 = -rho1 + 2 * p1, -rho2 + 2 * p2               # tr(theta T): time component continued
    integrand = (A - b1 * b2 / 3) / e2 - (b1 * b2 / 3) / (2 * e0 * f0_shift)
    w1, w2 = weights or (s1["E"], s2["E"])
    V = -(4 / math.pi) * np.trapezoid(np.sinc(Q * sep / math.pi) * integrand, Q)
    return V * sep / (w1 * w2)


def run(report=print):
    out = {}
    a, b = source(FRACS[0]), source(FRACS[1])

    # W1: the sources.
    for s in (a, b):
        assert abs(s["E_fields"] / s["E"] - 1) < 1e-6 and abs(s["laue"]) < 1e-6 and s["conservation"] < 1e-4
        assert abs(s["rho"][0] / s["E"] - 1) < 1e-6 and abs(s["s"][0]) < 1e-6 * s["E"]
        assert abs(s["rho"][-1]) < 1e-9 * s["E"]             # the grid holds the whole transform
    assert a["binding"] > 0 and b["binding"] < -0.2 and min(SEPARATIONS) > a["radius"] + b["radius"]
    out.update(binding=(a["binding"], b["binding"]), energy=(a["E"], b["E"]),
               sources=max(abs(s["E_fields"] / s["E"] - 1) for s in (a, b)),
               conservation=max(s["conservation"] for s in (a, b)), laue=max(abs(s["laue"]) for s in (a, b)))
    report(f"W1 sources: binding {100 * a['binding']:+.2f}% and {100 * b['binding']:+.2f}%. Int T_00 = E to "
           f"{out['sources']:.0e}, stress conserved to {out['conservation']:.0e}, int T_kk = 0 to {out['laue']:.0e}")

    # W2: Newton's law with the whole energy as weight.
    newton = [potential(a, b, x) for x in SEPARATIONS]
    charge = [potential(a, b, x, weights=(a["mQ"], b["mQ"])) for x in SEPARATIONS]
    expect_charge = -(1 + a["binding"]) * (1 + b["binding"])
    assert max(abs(v + 1) for v in newton) < 1e-6
    assert max(abs(v - expect_charge) for v in charge) < 1e-6 and abs(expect_charge + 1) > 0.2
    out.update(newton=newton, charge_weight=charge[0])
    report("W2 potential: V r/(G E1 E2) = " + ", ".join(f"{v:+.7f}" for v in newton)
           + f" at r = {SEPARATIONS}. In units of the charges' m Q it is {charge[0]:+.4f}")

    # W3: the stress does not weigh.
    bare = [potential(a, b, x, stress=False) for x in SEPARATIONS]
    out.update(stress_share=max(abs(u - v) for u, v in zip(newton, bare)))
    assert out["stress_share"] < 1e-6
    report(f"W3 stress: its cross terms move the potential by {out['stress_share']:.0e}")

    # W4: the coefficient is the kernel's.
    moved = potential(a, b, SEPARATIONS[0], f0_shift=1.1)
    assert abs(moved + 2 * (2 / 3 - 1 / (6 * 1.1))) < 1e-6 and abs(moved + 1) > 0.02
    out.update(moved_kernel=moved)
    report(f"W4 kernel: with F0/F2 = -2.2 the same sources give {moved:+.4f}")

    # W5: the regulator's reach.
    ratios, floor = kernel_ratios()
    assert min(floor) > 0                                    # no zero of F2 or F0 on the lumps' momenta
    reach = {e: max(abs(potential(a, b, x, ratios=ratios, eps=e) + 1) for x in SEPARATIONS) for e in EPS}
    assert max(reach.values()) < 1e-4
    out.update(reach=reach, kernel_floor=floor)
    report("W5 regulator: with the full form factors the potential moves by "
           + ", ".join(f"{v:.0e} at eps = {e}" for e, v in reach.items())
           + f"; F2 and F0 stay above {floor[0]:.2f} and {floor[1]:.2f} of their Einstein values to q^2 t0 = 30")
    report("RESULT: two lumps attract with -G E1 E2 / r, each weighing its whole energy, binding included; "
           "the coefficient is the induced kernel's F0/F2 = -2.")
    return out


if __name__ == "__main__":
    run()
