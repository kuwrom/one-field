"""
carrier_ledger.py: every node of the web against the closed carrier.

THE PICTURE.  There is one field.  Couplings emerge and structures
close.  The background is not a separate object.  It is the field's
common mode, the component that carries no label.  Once every
interference is accounted for, the closed objects appear in the common
mode as a bias, and that bias, measured, is gravity.  The searchlight
resolves the closed loop at depth three, and the web has one node
there, the e<->q carrier L = (alpha/2pi)^2.  Every number compared with
a measurement is the web's reading at depth four.  So the account of a
node is its account with that carrier.

THE QUESTION.  For each node of the solved web: does it read the
carrier?  With what weight and sign?  Where does the weight come from,
and what range of weight does the measurement of that node allow?  The
common mode carries no weight, so for it the question is how its count
responds to the carrier.  For a node that reads no carrier: does it
follow a reader, do the data exclude a weight, or can no measurement
resolve one?

GENERATED from root.WEB after the full chain.  The readers are the
nodes with terms whose parent is eq_loop.  The response of every node
is measured by scaling the carrier, solving the web again and then
restoring it.

STATEMENTS (each frozen by tests/audit/test_carrier_ledger.py).

  C1 (the readers).  Six nodes read the carrier: the coupling 1/alpha
     with weight 2, the quark action with 16, the three lepton rates
     with (27/2) Delta_k, and the common mode through the closure.  The
     common mode's count moves by -7.50 L through the three lepton
     terms (tests/probes/common_mode_count.py C3).  The ruler reads it
     through the electron's term.
  C2 (who follows).  Of the 27 other nodes, six read the carrier and
     ten follow a reader through the web's own reads.  The followers are
     sin2W and the lepton factor through the coupling, the ruler through
     the electron, the six quarks through their leptons and each other,
     and v_EW through the action and the ruler.  Eleven read none.  When
     the carrier is scaled, every reader and follower moves except the
     anchor m_e and the two quarks that read it by ratio.
  C3 (the signs, clause by clause).  The vertex composition rule
     gives each weighted reader its sign along with its weight.  The
     terms of the coupling, the action and the rates are negative on
     their ledgers.  1/alpha falls, so alpha rises.  S_quark falls.  The
     three lepton rates fall (m_e's own term, m_mu, m_tau).  The common
     mode's slope comes from the closure.  The carrier lowers the tau
     rate most, so the ratio of the harmonic power to the q = 0 power
     falls, and the count falls with it.
  C4 (the weight is the node's own for the leptons; the common mode
     carries none).  A lepton, a knot, reads the whole traversal
     dim J3(O) = d11^3 = 27 with the one channel's share C_EFF = 1/2,
     times its gap.  The common mode's two edges are the powers of the
     relative channels on the depth four rates.  On the searchlight's
     circulant they are 1/2 and 1/2, summing to 1.  The carrier reaches
     them only through the lepton terms.  The common mode has no
     response of its own, so no vertex weight is its reading
     (tests/probes/common_mode_count.py C2).  The coupling's 2 and the
     action's 16 follow their own clauses of the rule: the orientation
     count, and the charge trace times the Casimir.
  C5 (what the data allow).  The measurement that reads each weight
     allows 2.0000 +- 0.0001 (Cs recoil; 1.9993 by Rb recoil and
     1.9995 by the electron moment), 16.0 +- 0.2 (G_F), 13.50 +- 0.03
     (m_mu/m_e), and 33 +- 17 for the common mode's slope (G).  Under
     the recommended value, the count's -7.50 sits 2.43 standard
     deviations below that.  This distance is set by the rule that
     combines the measurements (tests/probes/common_mode_count.py C6).
  C6 (the nodes that read none).  A unit weight on the lepton action
     moves 1/alpha by over a thousand standard deviations through the
     lever arm.  The data exclude it, and the action's reading sits on
     S_quark.  A unit weight on sin2W or on the dark ratio is below a
     tenth of a standard deviation of the measurement that reads it.
     The quark references are a thousand times coarser than L.  So no
     measurement resolves a weight there.  The venting ledger is a
     count of rows.
  C7 (the readings that are not web nodes).  The mixing matrices are
     read from finite maps, not ledgers.  Every depth-4 term is either
     sector blind or diagonal in its sector's mass basis.  Neither kind
     moves an eigenvector, so V = U_u^dag U_d is the same at depth four,
     to 1e-12.  The one live path is the lepton masses inside the PMNS
     map, which moves each PMNS pull by less than 1e-3
     (tests/probes/mixing_depth4.py B2, B3, frozen in
     tests/mixing/test_mixing_probes.py).  eta_B reads J_lep from that
     map and the G2 face action, and the data keep the action off the
     carrier (C6).  The dark ratio is a ratio of the two shares of the
     bridge response.  It has no depth one edge, so the rule gives it
     weight zero.  Each of these is therefore compared at its depth
     four reading.
  C8 (each reader's sign, by its own data).  The opposite sign on
     1/alpha's term misses the Cs recoil value by 3.3e4 standard
     deviations.  On the quark action it misses G_F by 168, and on the
     leptons it misses m_mu/m_e by about 900.  Each tests its own clause
     of the rule (C3).  None is a test of the common mode, whose slope
     comes from the closure.

RESULT.  Six nodes read the carrier and ten follow it.  Five readers
carry the weight and sign of their own clause of the vertex
composition rule, and their own data hold their three weights to a
percent or better.  The common mode carries no weight.  It is a count
of the lepton closure.  The carrier moves it by -7.50 L through the
three lepton terms, fixed with no measured G.

Usage: python3 tests/audit/carrier_ledger.py     (about 20 s, runs the chain)
"""

import contextlib
import io
import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_here), "probes"))
import _bootstrap  # noqa: E402,F401
from interference import (couplings, dark_sector, gravity,  # noqa: E402
                          masses, mixing, root)

CARRIER = "eq_loop"
# Reference rows read by this audit alone: value, one standard deviation.
ALPHA = {"Cs recoil": (137.035999046, 27e-9), "Rb recoil": (137.035999206, 11e-9),
         "electron moment": (137.035999166, 15e-9)}
GF = (1.1663788e-5, 0.0000006e-5)            # GeV^-2
MU_E = (206.7682827, 0.0000046)              # m_mu / m_e
G_REL_SIGMA = 0.00015 / 6.67430
SIN2W_SIGMA, DARK_SIGMA = 0.00004, 0.0645


def chain():
    with contextlib.redirect_stdout(io.StringIO()):
        R = root.derive()
        m = masses.derive(R)
        mx = mixing.derive(R, m)
        gravity_out = gravity.derive(R, m, mx, couplings.derive(R, m))
        dark_sector.derive(R, gravity_out)
    return root.WEB


def readers(web):
    """node -> slope, per unit carrier, of the node's terms whose parent is the carrier, with the ledger's sign.

    For a weighted reader the slope is its weight. For the common mode, whose edges are a count of the
    depth four lepton rates, it is the count's slope through the closure."""
    L, off = web.state[CARRIER], dict(web.state, **{CARRIER: 0.0})
    out = {}
    for n, led in web.items():
        terms = [t for t in led.terms if t.parent == CARRIER]
        if terms:
            out[n] = sum(t.value(web.state) - t.value(off) for t in terms) / L
    return out


def response(web, eps=1e-3):
    """d ln(node) / d ln(carrier), in units of L, for every node. The web is restored."""
    term = web[CARRIER].terms[0]
    factor, before = term.factor, dict(web.state)
    term.factor = lambda s: (1.0 + eps) * factor(s)
    try:
        web.solve()
        after = dict(web.state)
    finally:
        term.factor = factor
        web.solve()
    assert all(abs(web.state[k] - v) <= 1e-12 * abs(v) for k, v in before.items())
    L = before[CARRIER]
    return {k: (after[k] / v - 1.0) / (eps * L) for k, v in before.items() if k != CARRIER and v}


def followers(web, read):
    """Nodes that reach a reader through the reads the web's own tracker records."""
    reads = {}
    for n, led in web.items():
        tracker = web._ReadTracker(web.state)
        led.value_at(tracker)
        reads[n] = tracker.reads - {n}
    reach = set(read)
    while True:
        more = {n for n, rs in reads.items() if n not in reach and n != CARRIER and rs & reach}
        if not more:
            return sorted(more | (reach - set(read)))
        reach |= more


def run(report=print):
    web = chain()
    st, L = web.state, web.state[CARRIER]
    d10, d11 = root.d10, root.d11

    # C1: the readers.
    read = readers(web)
    gaps = {n: read[n] / -(d11**3 / d10) for n in ("m_e", "m_mu", "m_tau")}
    inv0 = 2**9 / math.pi
    assert set(read) == {"inv_alpha", "S_quark", "m_e", "m_mu", "m_tau", "common_mode"}
    assert abs(read["inv_alpha"] / inv0 + 2) < 1e-12 and abs(read["S_quark"] + 16) < 1e-12
    count_slope = (st["common_mode"] - 1.0) / L                       # the count 3 Q_l - 1 = 1 + slope L
    assert abs(read["common_mode"] - count_slope) < 1e-8 and abs(read["common_mode"] + 7.502) < 0.001
    report("C1 readers of the carrier: " + ", ".join(f"{n} {w:+.4g}" for n, w in read.items()))

    # C2: who moves.
    resp = response(web)
    moved = sorted(n for n, v in resp.items() if abs(v) > 1e-6)
    follow = followers(web, read)
    still = sorted(set(resp) - set(read) - set(follow))
    held = sorted((set(read) | set(follow)) - set(moved))
    assert abs(resp["m_e"]) < 1e-6 and resp["M_Pl_MeV"] > 0.5      # the anchor holds. The ruler reads it
    assert set(moved) <= set(read) | set(follow) and held == ["m_d", "m_e", "m_u"]
    assert (len(read), len(follow), len(still)) == (6, 10, 11)
    report(f"C2 {len(read)} readers, {len(follow)} followers ({', '.join(follow)}), {len(still)} read none "
           f"({', '.join(still)}). Held by the anchor: {', '.join(held)}")

    # C3: the signs, clause by clause. Coupling, action and rates are negative. The common mode is a count.
    signs = {"1/alpha": resp["inv_alpha"], "S_quark": resp["S_quark"],
             "m_mu": resp["m_mu"], "m_tau": resp["m_tau"], "m_e (its own term)": read["m_e"]}
    assert all(v < 0 for v in signs.values()) and read["common_mode"] < 0
    assert abs(resp["common_mode"] - read["common_mode"]) < 0.01
    report("C3 the term on each ledger: " + ", ".join(f"{k} {v:+.3g}" for k, v in signs.items())
           + f". Common mode (the closure's count) {read['common_mode']:+.3g}")

    # C4: the common mode's and the leptons' weights from their own edges.
    cm = web["common_mode"]
    own = sum(t.value(dict(st, **{CARRIER: 0.0})) for t in cm.terms)   # the searchlight's circulant
    traversal = d11**3
    assert abs(own - 1.0) < 1e-14 and cm.base == 0.0 and len(cm.terms) == 2
    assert all(t.parent == CARRIER and t.depth == 4 for t in cm.terms)
    assert all(abs(read[n] + root.C_EFF * traversal * g) < 1e-9 for n, g in gaps.items())
    report(f"C4 common mode: base 0, two edges summing to {own:.0f} on the searchlight's circulant, no weight of "
           f"its own. Lepton: C_EFF {root.C_EFF} x {traversal} x gap")

    # C5: what the data allow for each weight.
    allow = {}
    for name, (value, sigma) in ALPHA.items():
        allow["1/alpha, " + name] = (2 + (st["inv_alpha"] - value) / (inv0 * L), sigma / (inv0 * L))
    v_meas = (math.sqrt(2) * GF[0]) ** -0.5
    allow["S_quark, G_F"] = (16 + math.log(v_meas / st["v_EW_GeV"]) / L, GF[1] / (2 * GF[0]) / L)
    lever = (gaps["m_mu"] - gaps["m_e"]) * L
    allow["leptons, m_mu/m_e"] = (d11**3 / d10 - (MU_E[0] / (st["m_mu"] / st["m_e"]) - 1) / lever,
                                  MU_E[1] / MU_E[0] / lever)
    allow["common mode, G"] = ((root.G_CODATA * root.COMMON_MODE / root.G_PRED - 1) / L, G_REL_SIGMA / L)
    assert abs(allow["1/alpha, Cs recoil"][0] - 2) < allow["1/alpha, Cs recoil"][1]
    assert abs(allow["S_quark, G_F"][0] - 16) < allow["S_quark, G_F"][1]
    assert abs(allow["leptons, m_mu/m_e"][0] - 13.5) < allow["leptons, m_mu/m_e"][1]
    assert abs((read["common_mode"] - allow["common mode, G"][0]) / allow["common mode, G"][1] + 2.43) < 0.01
    report("C5 the weight each measurement allows: "
           + "; ".join(f"{k} {v:.4f} +- {s:.4f}" for k, (v, s) in allow.items()))

    # C8: the opposite sign on each reader that is not G.
    x_flip = st["inv_alpha"]
    for _ in range(200):
        u = 1 / (2 * math.pi * x_flip)
        x_flip = inv0 * (1 - 1 / (2 * math.pi) + 2 * u * u)
    flip = {"1/alpha, Cs recoil": (x_flip - ALPHA["Cs recoil"][0]) / ALPHA["Cs recoil"][1],
            "S_quark, G_F": -32 * L / (GF[1] / (2 * GF[0])),
            "leptons, m_mu/m_e": (-(d11**3 / d10) - allow["leptons, m_mu/m_e"][0]) / allow["leptons, m_mu/m_e"][1]}
    assert abs(flip["1/alpha, Cs recoil"]) > 3e4 and abs(flip["S_quark, G_F"]) > 150
    assert abs(flip["leptons, m_mu/m_e"]) > 800
    report("C8 the opposite sign, in standard deviations of each reader's own data: "
           + ", ".join(f"{k} {v:+.3g}" for k, v in flip.items()))

    # C6: a unit weight on the nodes that read none.
    lever_arm = float(root.b0_G2) / (2 * math.pi) * float(root.singh_ratio) * st["inv_alpha"] / inv0
    unit = {"S_lepton (1/alpha, Cs recoil)": lever_arm * L / ALPHA["Cs recoil"][1],
            "sin2W": L / SIN2W_SIGMA,
            "dark_ratio": st["dark_ratio"] * L / DARK_SIGMA}
    assert unit["S_lepton (1/alpha, Cs recoil)"] > 1000 and max(unit["sin2W"], unit["dark_ratio"]) < 0.1
    report("C6 a unit weight, in standard deviations of its reader: "
           + ", ".join(f"{k} {v:.3g}" for k, v in unit.items()))

    report("RESULT: six nodes read the carrier and ten follow. Five carry the weight and sign of their own clause "
           "of the vertex composition rule. The common mode, a count of the lepton closure, carries none. The "
           "carrier moves it by -7.50 L through the three lepton terms.")
    return dict(readers=read, response=resp, moved=moved, followers=follow, still=still, signs=signs,
                own=own, traversal=traversal, allow=allow, unit=unit, flip=flip)


if __name__ == "__main__":
    run()
