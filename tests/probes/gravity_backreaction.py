"""Gravity's backreaction: on the masses, on the running of the boundary datum, and on itself.

THE QUESTION. Every closure vents into the common mode. Does the
common mode act back: on the quark and lepton rates, on the running
that carries the boundary datum alpha_G2 = 1/(24 pi) down the lever
arm, and on its own strength? Where the engine contains such a term,
the probe reports it.  Where it contains none, the probe computes what
the term would do and which reading excludes it.

PREMISES. P1 the engine's web after the full chain. P2 the exchange
coefficient of a closure of rate Omega through the common mode,
(Omega / M_b)^2 (common_mode_web.py K2). P3 a gravitational term in
the running of a gauge coupling, d ln g / d ln E = a0 G E^2 beside the
gauge term. Robinson and Wilczek (arXiv:hep-th/0509050 Eqs. 1, 19)
find "a0 = -3/pi" with M_P the full Planck mass, so a0 = -3/pi in
units of Newton's G. Pietrykowski (arXiv:hep-th/0606208) finds zero in
the harmonic gauge: the correction is "specific to this gauge, and
[is] devoid of any physical meaning". Ebert, Plefka and Rodigast
(arXiv:0710.1002): "all gravitational quadratic divergencies cancel in
cut off regularization and are trivially absent in dimensional
regularization". Toms (arXiv:1010.0793 Eq. 24) finds -1/pi with a
proper time cutoff, the lane's own regulator. Anber, Donoghue and
El-Houssieny (arXiv:1011.3229): "the quadratic cutoff dependence in
gravitational corrections does not signal the appearance of a running
coupling in physical processes but totally drops out as an unphysical
artifact". P4 the lever
arm: 1/alpha_G2(mu*) = 24 pi - (b0/2 pi)(S0 - C2(26)) = 32/pi, the
base of 1/alpha(0) (root.py). P5 the composite metric's own loop on
the one list (horizon_amount.py, ledger_kernel.py R5).

STATEMENTS (the numbers below check each one).

  R1 (the web has no path back). No node of the web reads a gravity
     node.  Only its own diagnostic ratio reads the venting sum. The ruler M_Pl_MeV is read by the three leptons and
     v_EW, and through them by every quark.
  R2 (on the masses). A closure's exchange coefficient is its rate
     over the boundary's, squared: 1.75e-45 for the electron, 2.0e-34
     for the top, the largest of the nine. The finest mass comparison
     the engine makes is 2.2e-8, twenty five orders above. A bias
     common to every clock is a factor on the ruler and not a record
     (gates_as_records.py U).
  R3 (on the running: what a gravitational term would do). With the
     term of P3 the datum the logarithmic running starts from is
     1/alpha_eff = 24 pi e^(a0) - (b0/2 pi) I(a0),
     I = int_0^1 (e^(a0 x) - 1)/(2 x) dx, in place of 24 pi. At
     a0 = -3/pi the exponent S0 falls from 44.41 to 17.48 and M_Pl by
     twelve orders. At the proper time cutoff's -1/pi, to 32.45 and by
     five orders. At zero nothing moves.
  R4 (two readers exclude the term). The lever arm's endpoint moves
     with the datum, so 1/alpha(0) depends on a0 without any measured
     G.  It bounds |a0| below 5e-11 at two standard deviations. G bounds
     it below 1e-6. In the engine's reading the datum is the value the
     logarithmic running starts from, read whole at the boundary.  Both
     readers agree with a0 = 0 to those parts.
  R5 (on itself, in the ledger). The metric's own loop adds -8/3 to
     the one list's stiffness on shell, 25.83 to 23.17 at the
     covariant contact, and -23/3 in the de Donder gauge. Gates S and
     M pass at every contact with either. The anchored G does not move
     with it: an identity (gravity_resolution.py R3).
  R6 (on itself, in the strength). The common mode reads the closed
     charged carrier at depth four (root.py, THE COMMON MODE) and
     contains no self-echo. A self-echo with the common mode's own loop factor,
     its coefficient over 2 pi, would move G by 16 percent per unit
     weight at first order and 2.5 percent at second. The recommended
     G bounds those weights below 6e-4 and 3e-3.

RESULT. The web has no path from the common mode back to a rate, and
none is needed to twenty five orders. A gravitational term in the
running at the size first proposed would move the ruler by twelve
orders, and at the proper time cutoff's size by five.  1/alpha(0)
excludes it to 5e-11 without G.  That bound agrees with the zero that
dimensional regularization and the momentum cutoff give.
Gravity's own loop is counted in the ledger, and the gates pass with
it.  G alone bounds its self-echo in the strength below a part in a
thousand of its natural size.

Scope: one loop running with a power law term of the stated form. The
self-echo is bounded by G alone, below a part in a thousand of its
natural size.

Usage: python3 tests/probes/gravity_backreaction.py     (about 30 s, runs the chain)
"""

import math
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "audit"))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402
import boundary_reading_scan as scan  # noqa: E402
import common_mode_web as cw  # noqa: E402

A0 = {"Robinson and Wilczek": -3.0 / math.pi, "Toms, proper time cutoff": -1.0 / math.pi,
      "dimensional regularization, momentum cutoff": 0.0}     # P3, in units of Newton's G
A0_RW = A0["Robinson and Wilczek"]
GRAVITY_NODES = ("a1_bridge", "Sigma_vent(UV)", "Sigma_vent(broken)", "G_ratio")
MASSES = ("m_e", "m_mu", "m_tau", "m_u", "m_d", "m_s", "m_c", "m_b", "m_t")
FINEST = 2.2e-8                              # relative uncertainty of m_mu/m_e, the finest mass comparison
INV_ALPHA_SIGMA = 2.0e-10                    # relative uncertainty of the caesium recoil 1/alpha
G_SIGMA = r.G_REL_SIGMA                      # relative uncertainty of the recommended G
METRIC_LOOP = (-8.0 / 3.0, -23.0 / 3.0)      # on shell, de Donder (ledger_kernel.py R5)
BAND = (0.6, 15.0)


def readers(web):
    """For every node, the nodes that read it, by the web's own read tracker."""
    out = {}
    for name, led in web.items():
        if not isinstance(led, r.Ledger):
            continue
        tracker = web._ReadTracker(web.state)
        try:
            led.value_at(tracker)
        except web._SOFT:
            pass
        for key in tracker.reads:
            out.setdefault(key, set()).add(name)
    return out


def datum(a0, x_b=1.0):
    """1/alpha_eff: the value the logarithmic running starts from when d ln g/d ln E carries a0 G E^2.

    y = 1/alpha obeys dy/dt = b0/2pi - 2 a0 x y with x = G mu^2 and dx/dt = 2x, so
    d(y e^(a0 x))/dt = (b0/2pi) e^(a0 x). Below the boundary y = y_eff - (b0/2pi) ln(M_b/mu).
    """
    b = float(r.b0_G2) / (2 * math.pi)
    tail = quad(lambda x: math.expm1(a0 * x) / (2 * x) if x > 0 else a0 / 2, 0.0, x_b)[0]
    return (1.0 / r.alpha_G2_Pl) * math.exp(a0 * x_b) - b * tail


def inv_alpha_of(base):
    """The web's cubic for 1/alpha(0) with the algebraic base moved: x^3 = I[(1 - 1/2pi) x^2 - 1/(2 pi^2)]."""
    f = lambda x: x**3 - base * ((1 - 1 / (2 * math.pi)) * x * x - 1 / (2 * math.pi**2))
    return brentq(f, 0.5 * base, base)


def moved(a0):
    """Relative shifts of S0, of M_Pl (in e-folds) and of 1/alpha(0) under the term a0."""
    b = float(r.b0_G2) / (2 * math.pi)
    dy = datum(a0) - 1.0 / r.alpha_G2_Pl
    S = r.S_lepton + dy / b                                # the same low scale is reached from the moved datum
    endpoint = 1.0 / r.alpha_G2_WZW + dy                   # the lever arm's length is held
    base = float(r.singh_ratio) * endpoint
    return dict(S=S, efolds=S - r.S_lepton,
                inv_alpha=inv_alpha_of(base) / inv_alpha_of(float(r.singh_ratio) / r.alpha_G2_WZW) - 1)


def run(report=print, build=True):
    out = {}
    if build and "Sigma_vent(UV)" not in r.WEB:
        cw._chain()
    web = r.WEB

    # R1: no path back.
    read_by = readers(web)
    back = {n: sorted(read_by.get(n, ())) for n in GRAVITY_NODES}
    ruler = sorted(read_by["M_Pl_MeV"])
    assert back["G_ratio"] == [] and back["Sigma_vent(UV)"] == ["G_ratio"]
    assert all(set(v) <= set(GRAVITY_NODES) for v in back.values())
    assert ruler == ["m_e", "m_mu", "m_tau", "v_EW_GeV"]
    out.update(readers_of_gravity=back, readers_of_ruler=ruler)
    report(f"R1 readers of the gravity nodes: {back}. Of the ruler: {ruler}")

    # R2: on the masses.
    coeff = {n: (web.state[n] / web.state["M_Pl_MeV"]) ** 2 for n in MASSES}
    largest = max(coeff, key=coeff.get)
    orders = math.log10(FINEST / coeff[largest])
    assert largest == "m_t" and abs(coeff["m_e"] / 1.75173e-45 - 1) < 1e-4 and orders > 25
    out.update(coefficients=coeff, largest=largest, orders_below=orders)
    report(f"R2 exchange coefficients: electron {coeff['m_e']:.3e}, top {coeff['m_t']:.2e}, the largest; "
           f"{orders:.1f} orders below the finest mass comparison")

    # R3: a gravitational term in the running.
    assert abs(datum(0.0) - 24 * math.pi) < 1e-12
    terms = {name: moved(a) for name, a in A0.items()}
    rw = terms["Robinson and Wilczek"]
    out.update(a0=A0, S_moved={k: v["S"] for k, v in terms.items()},
               orders={k: v["efolds"] / math.log(10) for k, v in terms.items()})
    assert rw["efolds"] / math.log(10) < -8 and terms["Toms, proper time cutoff"]["efolds"] / math.log(10) < -3
    assert terms["dimensional regularization, momentum cutoff"]["efolds"] == 0.0
    report("R3 a gravitational term in the running: " + "; ".join(
        f"{name}, a0 = {A0[name]:+.3f}: 1/alpha_eff = {datum(A0[name]):.2f}, S0 = {v['S']:.2f}, M_Pl by "
        f"10^{v['efolds'] / math.log(10):+.1f}" for name, v in terms.items())
        + f" (the engine: {24 * math.pi:.2f}, {r.S_lepton:.2f})")

    # R4: the two readers.
    eps = 1e-9
    d_alpha = (moved(eps)["inv_alpha"] - moved(-eps)["inv_alpha"]) / (2 * eps)
    d_S = (moved(eps)["efolds"] - moved(-eps)["efolds"]) / (2 * eps)
    held_alpha = 2 * INV_ALPHA_SIGMA / abs(d_alpha)
    g_room = abs(1e-6 * scan.readings()[1]["recommended"][0]) + 2 * G_SIGMA   # the residual and two sigma
    held_G = g_room / (2 * abs(d_S))
    assert held_alpha < 1e-10 and held_G < 2e-6 and abs(A0_RW) / held_alpha > 1e10
    out.update(dlnalpha_da0=d_alpha, efolds_per_a0=d_S, a0_held_by_alpha=held_alpha, a0_held_by_G=held_G)
    report(f"R4 1/alpha(0) moves by {d_alpha:+.2f} a0 and bounds |a0| < {held_alpha:.0e}; G bounds |a0| < "
           f"{held_G:.0e}. The proposed value is {abs(A0_RW) / held_alpha:.0e} times the first")

    # R5: the metric's own loop, in the ledger.
    one, xi = B.rows(r.d10, r.d11, r.n7, r.n26), float(r.xi_bridge)
    sig = {loop: [B.sigma(one, 0.0, c, xi) + loop for c in B.CONTACTS] for loop in (0.0, *METRIC_LOOP)}
    gates = {}
    for loop in METRIC_LOOP:
        S = all(B.sigma(one, 1.0, c, xi) + loop >= 2 * math.pi for c in B.CONTACTS)
        t0 = [B.anchored(one, c, xi, extra=loop)[0] / (2 * math.pi) for c in B.CONTACTS]
        gates[loop] = (S, all(BAND[0] < t < BAND[1] for t in t0))
    assert abs(sig[0.0][0] - 155 / 6) < 1e-9 and all(all(g) for g in gates.values())
    out.update(sigma_with_loop={k: v[0] for k, v in sig.items()}, gates_with_loop=gates)
    report(f"R5 the metric's own loop: stiffness {sig[0.0][0]:.2f} to {sig[METRIC_LOOP[0]][0]:.2f} on shell and "
           f"{sig[METRIC_LOOP[1]][0]:.2f} in the de Donder gauge at the covariant contact; S and M pass with each")

    # R6: a self-echo in the strength.
    own = 1.0 / (2 * math.pi)                              # the common mode's coefficient, 1, over 2 pi
    bound = {n: g_room / own**n for n in (1, 2)}
    assert bound[1] < 1e-3 and bound[2] < 5e-3
    out.update(self_loop_factor=own, self_weight_bound=bound)
    report(f"R6 a self-echo would move G by {100 * own:.0f}% per unit weight at first order and "
           f"{100 * own * own:.1f}% at second; G holds the weights below {bound[1]:.0e} and {bound[2]:.0e}")

    report("RESULT: no node reads the common mode back, and none is owed to twenty five orders; a gravitational "
           "term in the running is excluded by 1/alpha(0) to 5e-11 without G. Gravity's own loop is counted in "
           "the ledger, where the gates pass, and its self-echo in the strength is bounded by G alone.")
    return out


if __name__ == "__main__":
    run()
