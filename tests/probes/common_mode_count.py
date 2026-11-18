"""The common mode as the web's count: one reading survives the framework's own criteria.

THE QUESTION. The background is no separate object: "if you account for
every interference, the effect of the objects that close is seen in the
background as bias". Read as a count rather than a search, the common
mode carries no weight of its own. It is a Z3 count of the closure it
reads, made on the depth four web. Nine readings of that count are on
the table. Which survive the criteria the framework already holds, with
G read only at the end?

PREMISES. P1 the solved web after the full chain. Responses d ln(x)/dL
at the carrier L = (alpha/2pi)^2 by scaling the carrier and solving
again (common_mode_accounting.py). P2 gate U, derived with its failing
case: "a common factor is no record" (gates_as_records.py). P3 the
common mode's depth one definition: base 0, "the vacuum is the empty
field, with no amplitude of its own", and two edges, the relative
channels' powers |v_q|^2 = E[v^2], q = 1, 2 (root.py, THE COMMON MODE),
whose sum is P v^2 P = H/q0, the harmonic power over the q = 0 power of
the lepton circulant, and G M_b^2 = P v^2 P (boundary_whole.py B5,
c = 2/B^2). P4 every object is read at depth four. P5 a live base is
"the same wave reached through another node" (root.py, Ledger), and
the common mode reads no names (interference/__init__.py, THE HORIZON). The engine's common mode is this
count (root.py, THE COMMON MODE). The probe checks the node against it.

THE CANDIDATES. K1 the closure's whole shed power in M_b units (+30.76),
K2 its q = 0 part (+27.00), K3 its harmonic part (+34.51), K4 the
closure's mean rate over its q = 0 power, C_EFF (1 + P v^2 P), K5
P v^2 P on the depth four closure, K6 the vertex rule's carrier edge,
weight 1, K7 K5 with K6 added, K8 P v^2 P summed over every Koide
triplet, K9 the shares of the whole, which sum to one at every depth.

STATEMENTS (the numbers below check each one).

  C1 (gate U). A closure's records are the ratios of its rates: "a
     common factor is no record". A factor common to the three lepton
     rates, the depth two (1 - alpha/2 pi), moves K1 to K3 by -1161 ppm
     and leaves K4, K5, K7, K8 unchanged. The readings in M_b units read
     the closure's scale, which the anchor fixes, not its records. K1,
     K2, K3 fail.
  C2 (the node's own object). The node's base is 0 and its edges are
     the two relative powers, so its object is P v^2 P. K4 carries the
     q = 0 power as a half, and the empty field has no q = 0 power. K6
     reads the carrier through a vertex weight, and a node with no
     response of its own has no vertex weight. K9 does not move when
     the harmonic does (raising m_tau by 1e-6 moves P v^2 P and leaves
     K9 unchanged), so it reads nothing. K4, K6, K9 fail.
  C3 (the carrier once). The carrier reaches the depth four closure
     through the three lepton terms and moves K5 by -7.50 L. K7 reads
     it a second time through K6. K7 fails.
  C4 (the population). The node's edges are the channels of the lepton
     circulant. The quark triplets reach their rates through the
     leptons by ratio and through named WZW vents (Q_cbt, Q_scb). K8
     counts the same wave again and reads names. Their moments,
     3 Q - 1 = 1 + 6.9e-3 and 1 + 3.1e-3, put G over 400 standard
     deviations high. K8 fails.
  C5 (the survivor). K5 alone passes: G M_b^2 = P v^2 P = 3 Q_l - 1 on
     the depth four closure, 1 - 7.50 L, with Q_l 5.06 ppm below 2/3.
     The measured lepton masses give 1 + (-4.9 +- 11.3) L. The tau mass
     sets the spread, and the engine sits at -0.23 of it.
  C6 (the test). G = 6.673936e-11. That is -2.43 standard deviations
     from the recommended value, +0.17 with a spread added to each
     input, and -0.63 against the median (g_measurements.py, the same
     sixteen inputs).

RESULT. Counted rather than searched for, the common mode is the
second moment of the lepton circulant's relative part over its common
part, read on the depth four web: G M_b^2 = 3 Q_l - 1. Gate U removes
every reading in absolute units. The node's own definition removes all
the others but one. The carrier enters once, through the closure. The
reading has no weight and no edge to choose. Anything beyond it
belongs to the next depth.

Usage: python3 tests/probes/common_mode_count.py     (about 20 s)
"""

import math
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
sys.path.insert(0, os.path.join(os.path.dirname(_here), "audit"))
import _bootstrap  # noqa: E402,F401
import carrier_ledger  # noqa: E402
from common_mode_accounting import CLOSURES, scaled  # noqa: E402
import g_measurements as gm  # noqa: E402
from interference import root as r  # noqa: E402

LEPTONS = CLOSURES["leptons"]
MEASURED = {"m_e": 0.51099895069, "m_mu": 105.6583755, "m_tau": 1776.93}
TAU_SIGMA = 0.09
QUARK_Q = ("Q_cbt", "Q_scb")            # the crossing triplets' Koide ledgers: P v^2 P = 3 Q - 1
RULES = gm.comparison_rules()


def count(masses, closure=LEPTONS, unit=1.0):
    """Z3 count of a closure: q = 0 power, harmonic power, mean rate, and P v^2 P = H/q0."""
    amp = sum(s * math.sqrt(masses[k] / unit) for k, s in closure) / 3
    rate = sum(masses[k] / unit for k, _ in closure) / 3
    q0 = amp * amp
    return {"q0": q0, "H": rate - q0, "rate": rate, "Pv2P": (rate - q0) / q0}


def readings(state):
    """The nine candidates as values of G M_b^2 on a state (K1 to K3 in M_b units, as responses later)."""
    c = count(state)
    return {"K4": 0.5 * (1 + c["Pv2P"]), "K5": c["Pv2P"], "K7": c["Pv2P"] + state["eq_loop"],
            "K8": c["Pv2P"] + sum(3 * state[q] - 2 for q in QUARK_Q), "K9": 1.0}


def pulls(G):
    return {k: (G - v) / s for k, (v, s) in RULES.items()}


def run(report=print):
    web = carrier_ledger.chain()
    st, L = web.state, web.state["eq_loop"]
    out = {}

    # C1: gate U, a factor common to the three rates.
    f = st["lepton_EM"]
    common = dict(st, **{k: st[k] * f for k, _ in LEPTONS})
    absolute = {k: count(common, unit=st["M_Pl_MeV"])[x] / count(st, unit=st["M_Pl_MeV"])[x] - 1
                for k, x in (("K1", "rate"), ("K2", "q0"), ("K3", "H"))}
    ratio = {k: readings(common)[k] / readings(st)[k] - 1 for k in ("K4", "K5", "K7", "K8")}
    assert all(abs(v * 1e6 + 1161.4) < 0.1 for v in absolute.values())
    assert all(abs(v) < 1e-14 for v in ratio.values())
    out["gate_U"] = {**absolute, **ratio}
    report("C1 a common factor on the three rates: " + ", ".join(f"{k} {v * 1e6:+.1f} ppm" for k, v in out["gate_U"].items()))

    # C2: the node's own object; K9 reads nothing.
    eps = 1e-6
    moved = dict(st, m_tau=st["m_tau"] * (1 + eps))
    k5_moves = readings(moved)["K5"] / readings(st)["K5"] - 1
    k9_moves = readings(moved)["K9"] - readings(st)["K9"]
    assert abs(k5_moves) > 0.1 * eps and k9_moves == 0.0
    own = [round(t.value(dict(st, eq_loop=0.0)), 15) for t in web["common_mode"].terms]
    assert own == [0.5, 0.5] and web["common_mode"].base == 0.0
    assert abs(sum(t.value(st) for t in web["common_mode"].terms) - readings(st)["K5"]) < 1e-14
    report(f"C2 the node: base {web['common_mode'].base}, edges {own} on the searchlight's circulant. Its object P v^2 P; "
           f"K4 adds a q = 0 half, K6 a vertex weight, K9 moves by {k9_moves} when m_tau moves by {eps:.0e} and K5 by {k5_moves:+.1e}")

    # C3: the carrier once.
    before, after, dL = scaled(web)
    k5_resp = (readings(after)["K5"] / readings(before)["K5"] - 1) / dL
    assert abs(k5_resp + 7.50) < 0.01
    report(f"C3 the carrier moves K5 by {k5_resp:+.2f} L through the leptons; K7 adds it again: {k5_resp + 1:+.2f} L")

    # C4: the population.
    quark = {q: 3 * st[q] - 2 for q in QUARK_Q}
    k8_pull = pulls(r.G_DEPTH_ONE * readings(st)["K8"])["recommended"]
    assert abs(quark["Q_cbt"] - 6.944e-3) < 1e-5 and abs(quark["Q_scb"] - 3.086e-3) < 1e-5 and k8_pull > 400
    report("C4 quark moments " + ", ".join(f"{n} {v:+.3e}" for n, v in quark.items()) + f"; K8 at {k8_pull:+.0f} sigma")

    # C5: the survivor, and the measured lepton masses.
    k5 = readings(st)["K5"]
    Q = sum(st[k] for k, _ in LEPTONS) / sum(math.sqrt(st[k]) for k, _ in LEPTONS) ** 2
    meas = count(MEASURED)["Pv2P"]
    spread = abs(count(dict(MEASURED, m_tau=MEASURED["m_tau"] + TAU_SIGMA))["Pv2P"] - meas)
    pull_meas = ((k5 - 1) - (meas - 1)) / spread
    assert abs((k5 - 1) / L + 7.502) < 0.001 and abs(3 * Q - 1 - k5) < 1e-14
    assert abs((Q / (2 / 3) - 1) * 1e6 + 5.06) < 0.01 and abs(pull_meas + 0.23) < 0.01
    out.update(K5=k5, Q=Q, measured=meas, spread=spread, pull_measured=pull_meas)
    report(f"C5 K5 = 3 Q_l - 1 = 1 {(k5 - 1) / L:+.3f} L (Q_l {(Q / (2 / 3) - 1) * 1e6:+.2f} ppm of 2/3); "
           f"measured masses 1 {(meas - 1) / L:+.1f} +- {spread / L:.1f} L: {pull_meas:+.2f} sigma")

    # C6: the test, read last.
    G = r.G_DEPTH_ONE * k5
    p = pulls(G)
    assert abs(G / r.G_PRED - 1) < 1e-14 and abs(G - 6.673936e-11) < 2e-17 and abs(p["recommended"] + 2.43) < 0.01
    out.update(G=G, pulls=p)
    report(f"C6 G = {G:.6e}: " + ", ".join(f"{k} {v:+.2f}" for k, v in p.items()))
    return out


if __name__ == "__main__":
    run()
