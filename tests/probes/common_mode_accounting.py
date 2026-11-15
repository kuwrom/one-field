"""The common mode as the account of the closures, every object read at depth four.

THE QUESTION. The background is no separate object: "if you account for
every interference, the effect of the objects that close is seen in the
background as bias". Read that way the common mode carries no vertex
weight of its own. Its depth four content is what the closures do at
depth four. Which object of a closure does it read, over which closures,
and with which sign? The criteria are fixed before G is read: (1) every
closure the carrier moves is accounted for, each wave once; (2) the
object and the sign come from one classification; (3) the anchor's
inversion is not counted twice.

PREMISES. P1 the solved web after the full chain. Every response is
d ln(x)/dL at the carrier L = (alpha/2pi)^2, measured by scaling the
carrier and solving again (tests/audit/carrier_ledger.py). P2 a
closure's amplitudes psi_k = +-sqrt(m_k), read in M_b units. Its Z3
parts are the q = 0 power |mean psi|^2, the harmonic power (the q = 1, 2
channels) and their sum, the mean rate. P3 a live base is "the same
wave reached through another node" (root.py, Ledger). P4 the common
mode's depth one edges are the powers |v_q|^2 = E[v^2] of the relative
channels q = 1, 2 of the forgotten harmonic, normalized by Pv^2P = 1
(root.py, THE COMMON MODE).

STATEMENTS (the numbers below check each one).

  A1 (the closures). The carrier moves the lepton circulant through its
     three rates. The quark triplets move only as followers: m_c is
     (217/18) m_mu, m_t follows m_tau, m_b and m_s are Koide constraints
     on them. No quark rate reads the action's 16 (only v_EW does).
     Their responses, -32.3 for (c,b,t) and -37.3 for (s,c,b), are the
     lepton response seen through ratios. By P3 they are the same wave,
     so they are counted once. Counting them again gives 96.6 units.
  A2 (the objects). On the depth four lepton closure the q = 0 power
     falls by 27.00 L, the mean rate by 30.76 L, the harmonic power by
     34.51 L. The normalized harmonic power falls by 7.50 L: Pv^2P =
     1 - 1.01e-5. That is the object the common mode's own edges read
     (P4). The 27 is the q = 0 part alone. A transfer of what the
     closure sheds moves its whole power, 30.76. The node's own edges
     read 7.50 with their own sign. No classification gives the q = 0
     part alone.
  A3 (the anchor). M_b is inverted from the measured m_e with its
     depth four term in it, and moves by +0.55 L. The inversion fixes
     the unit and puts the shed power nowhere. Every response above is
     read in that moving unit, so nothing is counted twice.
  A4 (the test). With the node's edges read at depth four, common_mode
     = Pv^2P = 1 - 7.50 L and G = 6.673936e-11. That is -2.43
     standard deviations from the recommended value, +0.18 with a
     spread added to each input, and -0.63 against the median. The
     transfer of the whole shed power gives -0.13, +1.50, +0.87. The
     q = 0 part gives -0.36, +1.37, +0.72. The rule's weight 1 gives
     -1.92, +0.47, -0.29. The ranking reverses with the rule used to
     combine the measurements.

RESULT. The +27 has no classification behind it. It is the q = 0 part
of a shed whose whole is 30.76. One classification reads the common
mode's own edges on the depth four closure. That gives the normalized
harmonic power, with the sign of the edges: Pv^2P = 1 - 7.50 L. It
needs no transfer principle and no edge of its own to the carrier. The
transfer reading needs the principle that the shed power reappears in
the background. The depth two factor and the quark depth one terms
break that principle, unless losses through a named channel stay in
the background. G does not decide between the two readings.

Usage: python3 tests/probes/common_mode_accounting.py     (about 20 s)
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
sys.path.insert(0, os.path.join(os.path.dirname(_here), "audit"))
import _bootstrap  # noqa: E402,F401
import carrier_ledger  # noqa: E402
import g_measurements as gm  # noqa: E402
from interference import root as r  # noqa: E402

CLOSURES = {"leptons": (("m_e", 1), ("m_mu", 1), ("m_tau", 1)),
            "(c,b,t)": (("m_c", 1), ("m_b", 1), ("m_t", 1)),
            "(s,c,b)": (("m_s", -1), ("m_c", 1), ("m_b", 1))}
RULES = gm.comparison_rules()


def scaled(web, eps=1e-3):
    """The web's state, and the state with the carrier scaled by 1 + eps. The web is restored."""
    term = web["eq_loop"].terms[0]
    factor, before = term.factor, dict(web.state)
    term.factor = lambda s: (1.0 + eps) * factor(s)
    try:
        web.solve()
        after = dict(web.state)
    finally:
        term.factor = factor
        web.solve()
    return before, after, eps * before["eq_loop"]


def parts(state, closure):
    """q = 0 power, harmonic power and mean rate of a closure, in M_b units."""
    unit = state["M_Pl_MeV"]
    amp = sum(s * math.sqrt(state[k] / unit) for k, s in closure) / 3
    rate = sum(state[k] / unit for k, _ in closure) / 3
    return {"q0": amp * amp, "harmonic": rate - amp * amp, "rate": rate, "normalized": (rate - amp * amp) / (amp * amp)}


def pulls(weight, L):
    G = r.G_DEPTH_ONE * (1 + weight * L)
    return G, {k: (G - v) / s for k, (v, s) in RULES.items()}


def run(report=print):
    web = carrier_ledger.chain()
    before, after, dL = scaled(web)
    L = before["eq_loop"]
    resp = {}
    for name, closure in CLOSURES.items():
        p0, p1 = parts(before, closure), parts(after, closure)
        resp[name] = {k: (p1[k] / p0[k] - 1) / dL for k in p0}
    out = {"response": resp}

    # A1: the closures, each wave once.
    lep = resp["leptons"]
    counted_again = sum(resp[n]["q0"] for n in CLOSURES)
    assert abs(resp["(c,b,t)"]["q0"] + 32.30) < 0.02 and abs(resp["(s,c,b)"]["q0"] + 37.32) < 0.02
    assert abs(counted_again + 96.62) < 0.05
    rd = carrier_ledger.readers(web)
    assert not ({"m_c", "m_b", "m_t", "m_s", "m_u", "m_d"} & set(rd))
    report("A1 quark triplets follow the leptons: " + ", ".join(f"{n} {resp[n]['q0']:+.2f}" for n in CLOSURES)
           + f". Counted again {counted_again:+.1f}")

    # A2: the objects on the depth four lepton closure.
    assert abs(lep["q0"] + 27.00) < 0.01 and abs(lep["rate"] + 30.76) < 0.01
    assert abs(lep["harmonic"] + 34.51) < 0.01 and abs(lep["normalized"] + 7.50) < 0.01
    moment = parts(before, CLOSURES["leptons"])["normalized"]
    assert abs(moment - 1 + 1.012e-5) < 1e-8
    report("A2 the lepton closure per unit L: " + ", ".join(f"{k} {v:+.2f}" for k, v in lep.items())
           + f"; Pv^2P = 1 {moment - 1:+.3e}")

    # A3: the anchor moves the unit and places nothing.
    unit = (after["M_Pl_MeV"] / before["M_Pl_MeV"] - 1) / dL
    assert abs(unit - 0.547) < 0.005 and abs(after["m_e"] / before["m_e"] - 1) < 1e-12
    report(f"A3 M_b moves by {unit:+.3f} L with m_e held. Every response is read in that unit")

    # A4: G under each reading.
    readings = {"own edges at depth four": lep["normalized"], "rule's weight 1": 1.0,
                "q = 0 part of the shed": -lep["q0"], "whole shed (transfer)": -lep["rate"]}
    tests = {k: pulls(w, L) for k, w in readings.items()}
    own = tests["own edges at depth four"]
    assert abs(own[0] - 6.673936e-11) < 2e-17 and abs(own[1]["recommended"] + 2.43) < 0.01
    assert abs(tests["whole shed (transfer)"][1]["recommended"] + 0.13) < 0.01
    assert abs(tests["q = 0 part of the shed"][1]["recommended"] + 0.36) < 0.01
    out.update(readings=readings, tests=tests, counted_again=counted_again, unit=unit)
    for k, (G, p) in tests.items():
        report(f"A4 {k:26s} w {readings[k]:+7.2f}  G = {G:.6e}  "
               + ", ".join(f"{n} {v:+.2f}" for n, v in p.items()))
    return out


if __name__ == "__main__":
    run()
