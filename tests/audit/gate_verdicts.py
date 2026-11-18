"""
gate_verdicts.py: the gravity lane's six verdicts, computed from the probes that decide them.

The engine records each gate's status together with the tests that
witness it (interference/gravity.py).  This module runs the deciding
probes and reduces each gate to a boolean on named numbers, so each
status is checked against a computation.
tests/common_mode/test_gate_verdicts.py fails when the engine says PASS
and the decider does not. It also shows each decider failing on a rival
reading that the repository computes:
  C, U  a frame vector in the vacuum;
  S, M  the list without the amount;
  T     the amount carrying the label's coupling xi = 1/(48 pi), whose
        kernel has a zero inside the boundary at 1.268 M_b
        (tests/probes/horizon_amount.py q_zero_label);
  P     the amount treated as massive below the boundary. It passes S
        and M, but its stiffness depends on that mass at order one
        (horizon_amount.py A7).
The gates stated as claims about records are in
tests/probes/gates_as_records.py.

  S  the one list's stiffness reaches 2 pi at every contact
  M  its anchored regulator is Planckian, 0.6 < t0 < 15, at every contact
  C  the saddle's symmetry leaves no vector and one rank-2 tensor, and
     the scalar, vector and Weyl rows each have one symbol on its cone
  U  a row's weight is its energy, and two lumps exchange the kernel
     with their whole energies as weights
  T  the kernel has the Einstein form and sign on its axis, and the
     induced action's cubic term is Einstein's
  P  the stress has no monopole, gamma = beta = 1, and the stiffness
     depends on the rows' counts, not on their masses

Usage: python3 tests/audit/gate_verdicts.py     (about half a minute)
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "probes"))
import _bootstrap  # noqa: E402,F401
from interference import boundary as B  # noqa: E402
from interference import root as r  # noqa: E402

NAMES = ("S sign", "U universality", "M magnitude", "C one cone", "T tensor carrier", "P no preferred location")


def _silent(*args, **kwargs):
    pass


def evidence():
    """Run the deciding probes once. Returns their outputs by probe name."""
    import horizon_amount
    import induced_cubic
    import post_newtonian
    import saddle_stabilizer
    import two_lumps
    one = B.rows(r.d10, r.d11, r.n7, r.n26)
    xi = float(r.xi_bridge)
    roots = {c: B.anchored(one, c, xi) for c in B.CONTACTS}
    return dict(
        list=dict(inside=[B.inside(one, c, xi) for c in B.CONTACTS],
                  t0=[(roots[c][0] / (2 * math.pi)) if roots[c] else float("nan") for c in B.CONTACTS]),
        saddle=saddle_stabilizer.run(report=_silent), lumps=two_lumps.run(report=_silent),
        amount=horizon_amount.run(report=_silent), cubic=induced_cubic.run(report=_silent),
        pn=post_newtonian.run(report=_silent))


def verdicts(ev):
    """Each gate as a boolean on the evidence."""
    s, l, a, c, pn = ev["saddle"], ev["lumps"], ev["amount"], ev["cubic"], ev["pn"]
    counts = s["invariant_counts"]
    one_symbol = (s["stabilizer_dim"] == 6 and counts["rank 1"] == 0 and counts["rank 2"] == 1
                  and s["weyl_invariants"] == 1 and s["vector_det"]["null"] < 1e-12)
    out = {
        "S sign": all(ev["list"]["inside"]),
        "M magnitude": all(0.6 < t < 15.0 for t in ev["list"]["t0"]),
        "C one cone": one_symbol,
        "U universality": (one_symbol and s["weight_over_energy"] < 1e-12 and s["weyl_vector_error"] < 1e-12
                           and max(abs(v + 1) for v in l["newton"]) < 1e-6 and l["sources"] < 1e-6),
        "T tensor carrier": (a["zero"] is None and a["kernel_limit"] > 0 and abs(a["kernel_shape"] + 2) < 1e-3
                             and abs(a["static_sign"] - 1) < 1e-3 and c["expansion_error"] < 1e-6
                             and c["odd_error"] < 1e-5),
        "P no preferred location": (l["stress_share"] < 1e-6 and l["laue"] < 1e-6 and pn["gamma"] == 1
                                    and pn["beta"] == 1 and abs(a["mass_response_one"]) < 1e-6),
    }
    return {name: "PASS" if out[name] else "FAIL" for name in NAMES}


def framed(ev):
    """The evidence with the saddle's numbers replaced by those of a vacuum that holds a frame vector."""
    s = dict(ev["saddle"])
    s.update(stabilizer_dim=s["frame_stabilizer_dim"], weyl_invariants=s["frame_weyl_invariants"],
             invariant_counts=dict(s["invariant_counts"], **{"rank 1": s["frame_vectors"],
                                                             "rank 2": s["frame_symmetric_rank2"]}),
             weight_over_energy=1.0 - min(s["frame_weight"]))
    return dict(ev, saddle=s)


def failing(ev):
    """Each gate's verdict on the rival reading that fails it: gate -> (the case, the verdict)."""
    import gates_as_records as gr
    import horizon_amount
    bare = gr.grain(gr.with_amount(None))
    q_label = horizon_amount.zero_q(r.n7 * r.n26, float(r.xi_bridge))
    assert abs(q_label - ev["amount"]["q_zero_label"]) < 1e-12
    labelled = dict(ev["amount"], zero=q_label)
    cases = {
        "S sign": ("the list without the amount", dict(ev, list=bare)),
        "U universality": ("a frame vector in the vacuum", framed(ev)),
        "M magnitude": ("the list without the amount", dict(ev, list=bare)),
        "C one cone": ("a frame vector in the vacuum", framed(ev)),
        "T tensor carrier": ("the amount with the label's coupling", dict(ev, amount=labelled)),
        "P no preferred location": ("the amount read as massive below the boundary",
                                    dict(ev, amount=dict(ev["amount"],
                                                         mass_response_one=ev["amount"]["mass_response_massive"]))),
    }
    return {name: (case, verdicts(e)[name]) for name, (case, e) in cases.items()}


def run(report=print):
    ev = evidence()
    v, vf, fail = verdicts(ev), verdicts(framed(ev)), failing(ev)
    for name in NAMES:
        case, verdict = fail.get(name, ("no computed rival reading", "-"))
        report(f"{v[name]:4s} {name:24s} {case}: {verdict}")
    return dict(evidence=ev, verdicts=v, framed=vf, failing=fail)


if __name__ == "__main__":
    run()
