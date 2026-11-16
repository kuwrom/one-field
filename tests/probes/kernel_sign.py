"""The sign of the one list's spin-2 form factor, in ball arithmetic.

THE QUESTION. one_list_kernel.py evaluates the summed kernel with a
fixed set of Feynman parameter nodes and scans its spin-2 form factor
F2 for zeros on a grid. A grid proves no inequality, and the nodes lose
the endpoints at large momentum. In closed form, with rigorous error
bounds, where is F2 positive?

PREMISES. P1 the one list's massless rows under the proper time cutoff
(boundary.py; ledger_kernel.py's heat trace form factors): N_v = 18
vector fields with their ghosts, N_D = 22.5 Dirac fields, N_s = 182
real scalars at coupling xi. P2 on a transverse traceless perturbation
the kernel of ledger_kernel.py reduces, with y = q^2 t0,
Z = (y/4)(1 - v^2) and P_k = int_0^1 E_k(Z) dv, to
    F2 = (C16 / (2 t0^2)) Phi(y),
    Phi = -N_v G_v + N_D G_D - N_s [xi y/4 + (P_3 - 1/2)/4],
    G_v = y/4 + (P_3 - 1/2)/2 + y^2 P_1/8 + y (P_2 - 1)/2,
    G_D = y/4 + (P_3 - 1/2) + y (P_2 - 1)/4.
The terms linear in y cancel in G_v and G_D. With
B = (y/2) int_0^1 E_1(Z) Z/(1 + v) dv the cancellation is done by hand,
    G_v = B - (y/4) E_2(y/4) + (y/2) P_2 - 1/4 + P_3/2,
    G_D = P_3 - 1/2 + (y/4) P_2,
and for large y the same integrals are taken over z = Z with the weight
(1 - 4z/y)^(-1/2), where y enters only through that weight. P3 Arb
(python flint): every integral is a rigorous enclosure, y may be an
interval, and an inequality holds for every point of it.

STATEMENTS (the numbers below check each one).

  S1 (the closed form is the kernel). At y = 1e-3 to 10 Phi equals
     one_list_kernel.py's F2 to 1e-6. The two forms of P2 agree to 1e-9
     at y = 301, 1e3 and 5e4.
  S2 (the nodes' reach). At y = 5e4 the node kernel gives a negative
     F2 where Phi = 8.1241: the zero that kernel shows near y = 3.4e4
     is the nodes' loss of the endpoints, not the form factor's.
  S3 (the sign, xi = 0). Phi > 0 for every y in [1e-3, 1e12]: the axis
     is covered by intervals on each of which the enclosure of Phi is
     positive, its least lower bound 0.0034 at the lower end.
  S4 (the limit). As y grows Phi tends to (N_s - 2 N_D - 4 N_v)/8, the
     one list's 65/8 = 8.125. The enclosure at y in [1e9, 1e12] contains
     it. The limit is positive when the scalars outnumber
     2 N_D + 4 N_v = 117: with 110 scalars Phi is negative at y = 1e6.
  S5 (a coupling xi). xi takes N_s xi y/4 from Phi, a term that grows
     without bound, so every xi > 0 puts a zero on the axis: for the
     label's 1/(48 pi) it lies in 6.5 < y < 6.7 (q = 1.295 M_b at that
     list's regulator), sign certified on both sides.

RESULT. The one list's spin-2 form factor has the Einstein sign on the
whole Euclidean axis from 1e-3 to 1e12 in q^2 t0 and tends to the
positive constant 65/8 C16/(2 t0^2).  With the amount minimally
coupled the kernel has no zero. The amount's count sits above the 117
the limit needs, and any positive xi on the amount brings a zero in.

Scope: Euclidean momenta, massless rows, the covariant contact (two
ghosts per vector). Below 1e-3 Phi is Sigma y/4 to first order.

Usage: python3 tests/probes/kernel_sign.py     (about a minute; needs python flint)
"""

import math
import os
import sys

from flint import acb, arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402

ctx.prec = 100
DELTA = arb(2)**-60                          # the cut at the integrand's nonanalytic endpoint
TOL = arb(2)**-50
NAN = acb(float("nan"))
UNIT = arb(0.5, 0.5)                         # the interval [0, 1]
SWITCH = 300                                 # y above which the z-form is used
ONE = B.rows(r.d10, r.d11, r.n7, r.n26)
COUNTS = (sum(n for _, n, k, _, rd in ONE if k == "vector" and rd == "light") // 2,
          arb(sum(n for _, n, k, _, rd in ONE if k == "weyl" and rd == "light")) / 4,
          sum(n for _, n, k, _, rd in ONE if k == "bridge"))
XI = arb(1) / (48 * arb.pi())                # the label's coupling


def _combine(Bv, A2, A3y, edge, y, counts, xi):
    """Phi from B, A2 = (y/2) P2, A3y = P3/2 and edge = (y/4) E2(y/4)."""
    nv, nd, ns = counts
    Gv = Bv - edge + A2 - arb(1) / 4 + A3y
    Gd = 2 * A3y - arb(1) / 2 + A2 / 2
    return -nv * Gv + nd * Gd - arb(ns) * (xi * y / 4 + (2 * A3y - arb(1) / 2) / 4)


def _pieces_v(y):
    """B, (y/2) P2, P3/2 from integrals over v in [0, 1]."""
    def I(k, extra=False):
        def f(v, analytic):
            Z = (1 - v * v) * y / 4
            if analytic and not (Z.real > 0):
                return NAN
            return Z.expint(k) * (Z / (1 + v) if extra else 1)
        return acb.integral(f, 0, 1 - DELTA, rel_tol=TOL, abs_tol=TOL).real
    # beyond the cut: E_2 <= 1, E_3 <= 1/2, Z E_1(Z) <= 1
    return ((y / 2) * (I(1, True) + DELTA * UNIT), (y / 2) * (I(2) + DELTA * UNIT),
            (I(3) + DELTA * UNIT / 2) / 2)


def _pieces_z(y):
    """The same three from integrals over z in [0, y_low/8] with weight w = (1 - 4z/y)^(-1/2)."""
    c = y.lower() / 8

    def I(k, kind):
        def f(z, analytic):
            u = 1 - 4 * z / y
            if analytic and not (z.real > 0 and u.real > 0):
                return NAN
            w = u.rsqrt()
            return z.expint(k) * (w if kind == "A" else z * w * w / (w + 1))
        return acb.integral(f, DELTA, c, rel_tol=TOL, abs_tol=TOL).real
    Ec = {k: acb(c).expint(k).real for k in (1, 2, 3)}
    # below the cut the integrands are at most 2, 1 and 1. Beyond c they are bounded by E_k(c)
    A2 = I(2, "A") + DELTA * 2 * UNIT + Ec[2] * (y / 2) * UNIT
    A3 = I(3, "A") + DELTA * UNIT + Ec[3] * (y / 2) * UNIT
    return I(1, "B") + DELTA * UNIT + Ec[1] * (y * y / 8) * UNIT, A2, A3 / y


def phi(y, counts=COUNTS, xi=0, form=None):
    """Enclosure of Phi on the ball y."""
    y = arb(y)
    use_z = (y.lower() > SWITCH) if form is None else form == "z"
    Bv, A2, A3y = _pieces_z(y) if use_z else _pieces_v(y)
    edge = (y / 4) * acb(y / 4).expint(2).real
    return _combine(Bv, A2, A3y, edge, y, counts, arb(xi))


def certify(lo, hi, stats=None, depth=0, **kw):
    """Prove Phi > 0 on [lo, hi] by bisection in log y. Stats holds the intervals used and the least bound."""
    stats = {"intervals": 0, "least": math.inf} if stats is None else stats
    p = phi(arb(lo).union(arb(hi)), **kw)
    stats["intervals"] += 1
    if p.is_finite() and p > 0:
        stats["least"] = min(stats["least"], float(p.lower()))
        return stats
    assert depth < 40, (lo, hi)
    mid = math.sqrt(lo * hi)
    certify(lo, mid, stats, depth + 1, **kw)
    return certify(mid, hi, stats, depth + 1, **kw)


def node_kernel(y):
    """one_list_kernel.py's F2 in Phi's units, at the list's anchored regulator."""
    import one_list_kernel as olk
    t0 = B.anchored(ONE, 1.0, float(r.xi_bridge))[0] / (2 * math.pi)
    F2 = olk.form_factors(olk.kernel_rows(ONE), t0, y)[0]
    return F2 * t0 * t0 * 32 * math.pi**2


def run(report=print, full=True):
    out = {}
    assert COUNTS[0] == 18 and COUNTS[2] == 182 and abs(float(COUNTS[1]) - 22.5) < 1e-12

    # S1, S2: the closed form against the node kernel, and the two forms against each other.
    low = max(abs(node_kernel(y) / float(phi(y).mid()) - 1) for y in (1e-3, 1e-2, 0.1, 1.0, 10.0))
    forms = max(float(abs(phi(y, form="v") - phi(y, form="z")).upper()) for y in (301.0, 1e3, 5e4))
    far_node, far = node_kernel(5e4), phi(5e4)
    assert low < 1e-6 and forms < 1e-8 and far_node < 0 < far and abs(float(far.mid()) - 8.1241) < 1e-4
    out.update(node_agreement=low, form_agreement=forms, node_at_5e4=far_node, phi_at_5e4=float(far.mid()))
    report(f"S1 closed form: equals the node kernel to {low:.0e} for y <= 10. The two forms agree to {forms:.0e}.  "
           f"S2 at y = 5e4 the node kernel gives {far_node:+.2f} where Phi = {float(far.mid()):.4f}")

    # S3: the sign on the axis.
    spans = ((1e-3, 1e3), (1e3, 1e12)) if full else ((1e-3, 1e-1), (1e3, 1e12))
    cert = [certify(a, b) for a, b in spans]
    out.update(certified=spans, intervals=sum(c["intervals"] for c in cert), least=min(c["least"] for c in cert))
    assert out["least"] > 0
    report("S3 sign: Phi > 0 on " + " and ".join(f"[{a:g}, {b:g}]" for a, b in spans)
           + f", {out['intervals']} intervals, least lower bound {out['least']:.4f}")

    # S4: the limit and the count it needs.
    nv, nd, ns = COUNTS
    limit = (ns - 2 * float(nd) - 4 * nv) / 8
    tail = phi(arb(1e9).union(arb(1e12)))
    need = 2 * float(nd) + 4 * nv
    few = phi(1e6, counts=(nv, nd, 110))
    assert limit == 65 / 8 and tail.contains(arb(65) / 8) and float(tail.rad()) < 1e-3 and need == 117 and few < 0
    out.update(limit=limit, scalars_needed=need, phi_110=float(few.mid()))
    report(f"S4 limit: Phi -> (N_s - 2 N_D - 4 N_v)/8 = {limit}. Positive for N_s > {need:.0f}. With 110 scalars "
           f"Phi(1e6) = {float(few.mid()):+.3f}")

    # S5: the label's coupling on the amount.
    below, above = phi(arb(6.5), xi=XI), phi(arb(6.7), xi=XI)
    assert below > 0 > above
    out.update(label_zero=(6.5, 6.7))
    report(f"S5 with xi = 1/(48 pi): Phi(6.5) = {float(below.mid()):+.4f}, Phi(6.7) = {float(above.mid()):+.4f}: "
           "a zero between")
    report("RESULT: minimally coupled, the one list's spin-2 form factor has no zero from q^2 t0 = 1e-3 to 1e12 "
           "and tends to 65/8. A positive xi on the amount brings one in.")
    return out


if __name__ == "__main__":
    run()
