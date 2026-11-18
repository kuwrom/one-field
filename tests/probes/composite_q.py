"""Lambda route C: a composite vacuum variable q.

THE MECHANISM. Klinkhamer and Volovik (arXiv:0711.3170) give the
gravitating vacuum energy as eps - q d eps/dq (Eq. 28), zero in a
self-sustained equilibrium vacuum (Eqs. 30, 31). q must take the form
q g^{mu nu} in the vacuum, be conserved, carry no global charge, and
have a nonlinear eps(q) with nonnegative compressibility.  In their
minimal volume example (Sec. V E) "it is the dimensionless vacuum
variable q^ which adjusts itself to zero pressure", q = q^ / l^4 against
a fixed microscopic 4-volume l^4.

C1 (candidates the carrier offers).
   - the U(1) charge density: a global charge, frame dependent, zero in
     the empty field (excluded by KV: "any global charge ... must be
     zero in the Lorentz invariant vacuum");
   - dJ ^ dJ: zero in expectation (composite_volume.py A4);
   - <h> = (I1/4) g: the right form, but h is not a conserved tensor
     and the field's energy momentum tensor vanishes identically;
   - the regulator's grain, q = 1/t0^2, a minimal 4-volume: Lorentz
     invariant (a proper time), of the form q g, no charge.  The one
     survivor.
C2 (the grain is a unit, not a variable). In the empty field the grain
   is the only scale, and the anchor reads it in the boundary's units,
   t0 M_b^2 = Sigma/(2 pi).  Rescaling the grain with the boundary
   (t0 -> t0/lambda, M_b -> M_b sqrt lambda, every mass with it) leaves
   every dimensionless observable fixed: G M_b^2, Lambda G, m^2 t0.  The
   numbers check only dimensional analysis, an identity by
   construction. The content is the premise behind them: the empty
   field carries no second scale (no medium, every scale read off the
   anchor).  Under that premise a change of the grain is a change of
   units.  KV's mechanism needs q to move against a fixed microscopic
   scale, and the empty field has none.  So the dimensionless q^ has
   no realization here.
C3 (what a q-reading of the grain would leave). Per real component the
   regulated vacuum energy is eps(q) = -q E3(m^2 t0)/(32 pi^2), q = t0^-2,
   and eps - q d eps/dq = m^2 E2(m^2 t0)/(64 pi^2 t0).  The quartic
   class cancels identically.  The dimension two moment remains. It is
   the moment of the Higgs ratio (regulator_free_ratios.py R1, R3).  It
   vanishes for content massless at the regulator.  This is
   conditional on a q that, by C2, the theory does not have.  It is
   recorded as the form such a premise would take.

RESULT. Route C is excluded within the premises. It needs a dimensionless
density of records per grain as a variable of the state, and the state
carries none.

Usage: python3 tests/probes/composite_q.py     (a second)
"""

import math
import os
import sys

from scipy.special import expn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bootstrap  # noqa: E402,F401
import composite_pole as cp  # noqa: E402

N = cp.N_COMPONENTS


def observables(t0, m2, V0):
    """Dimensionless observables of the carrier sector at regulator t0, mass m^2, constant V0."""
    _, rho, stiff = cp.proper_time(m2, t0)
    G = 1.0 / (16 * math.pi * N * stiff)
    return dict(G_over_t0=G / t0, LambdaG=8 * math.pi * G * G * (V0 + N * rho), m2t0=m2 * t0)


def c2():
    base = observables(1e-2, 1.0, -3.0 * N * cp.proper_time(1.0, 1e-2)[1])
    out = {}
    for lam in (0.1, 3.0, 50.0):
        # grain t0 -> t0/lam. Every mass^2 -> lam m^2. Every energy density -> lam^2
        V0 = -3.0 * N * cp.proper_time(1.0, 1e-2)[1] * lam * lam
        o = observables(1e-2 / lam, lam * 1.0, V0)
        out[lam] = max(abs(o[k] / base[k] - 1) for k in base)
    return out


def eps_tilde_numeric(m2, t0, d=1e-6):
    q = t0**-2
    eps = lambda qq: cp.proper_time(m2, qq**-0.5)[1]
    deps = (eps(q * (1 + d)) - eps(q * (1 - d))) / (2 * q * d)
    return eps(q) - q * deps


def c3():
    rows = {}
    for m2, t0 in ((0.0, 1.0), (1e-3, 1.0), (0.3, 2.0), (2.0, 0.5)):
        num = eps_tilde_numeric(m2, t0)
        closed = m2 * float(expn(2, m2 * t0)) / (64 * math.pi**2 * t0)
        quartic = abs(cp.proper_time(m2, t0)[1])
        rows[(m2, t0)] = dict(numeric=num, closed=closed, quartic=quartic)
    return rows


def run(report=print):
    out = dict(c2=c2(), c3=c3())
    for lam, dev in out["c2"].items():
        assert dev < 1e-12, lam
    report("C2 rescaling the grain with the boundary and every mass: largest change of "
           "G/t0, Lambda G, m^2 t0 = " + ", ".join(f"{d:.1e} (x{l:g})" for l, d in out["c2"].items())
           + " (dimensional analysis, an identity).  With no second scale in the empty field the "
           "grain is a unit: nothing for q to move against.")
    for k, r in out["c3"].items():
        assert abs(r["numeric"] - r["closed"]) < 1e-6 * r["quartic"] + 1e-14, k
    assert abs(out["c3"][(0.0, 1.0)]["numeric"]) < 1e-12
    report("C3 eps - q deps/dq against m^2 E2(m^2 t0)/(64 pi^2 t0): "
           + "; ".join(f"m^2 t0 = {k[0] * k[1]:g}: {r['numeric']:+.4e} vs {r['closed']:+.4e} "
                       f"(quartic class {r['quartic']:.2e})" for k, r in out["c3"].items())
           + ".  Massless content leaves zero. Massive content leaves the dimension two moment.")
    report("RESULT: route C is excluded within the premises. The grain is a "
           "unit, and a dimensionless record density per grain would be a new premise.")
    return out


if __name__ == "__main__":
    run()
