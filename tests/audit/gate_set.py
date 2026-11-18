"""
gate_set.py: the closure gates themselves, each with its premise and the rival it removes.

THE QUESTION. The wiring scan (interference/wiring.py) and the coherence
scan (interference/searchlight/embedding_uniqueness.py) each close on one
survivor under a set of gates.  This module asks the next question: does
each gate earn its place?  A gate is not a free choice when it states a
premise of the chain and does work, meaning that removing it lets a rival
through.

THE STATEMENTS (computed below, over the same complete menus the scans use).

  G1 (what each gate removes).  For each gate the scan is run with that
     gate removed and the rivals that pass are printed.  A gate SELECTS
     when its removal lets a rival through.  It SHARES the work with
     another gate when removing the two together lets a rival through, so
     either one is enough.  It is a CONSEQUENCE when the survivor of the
     others satisfies it on this menu and it removes no rival.
  G2 (every gate names its premise).  Each gate carries its premise: R1,
     "the definition of what 'zero free parameters' means for a lattice
     vertex operator algebra" (the Octavian Coherence Gate paper, R1, whose
     wording the lattice premises quote), self-presence, or quantum numbers
     crystallising on breaking rather than being written on the vacuum.

The menus are the scans' own: the lattices of rank 8 and the even
unimodular lattices of rank 16 (Witt, Milnor), the classified conformal
pairs of E8(1), the levels k = 1..12, and the coherence table of every
rank-8 parent with every pair of distinct table algebras of total rank at
most 8.  No measured value enters.

Usage: python3 tests/audit/gate_set.py     (instant)
"""

import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "probes"))
import _bootstrap  # noqa: E402,F401
from interference import wiring as W  # noqa: E402
from interference.searchlight import embedding_uniqueness as EU  # noqa: E402


def _connected(G):
    """The Gram matrix's basis graph is connected (an irreducible lattice in this basis)."""
    n = len(G)
    seen, stack = {0}, [0]
    while stack:
        i = stack.pop()
        for j in range(n):
            if j not in seen and abs(G[i, j]) > 1e-12:
                seen.add(j)
                stack.append(j)
    return len(seen) == n


def _lattice_facts(name, G, rank):
    if G is None:                       # rank 16: even unimodular by classification
        return dict(integral=True, even=True, unimodular=True, positive=True,
                    irreducible=(name != "E8+E8"), minimal_rank=(rank == 8))
    return dict(integral=bool(np.allclose(G, np.round(G))),
                even=all(int(round(G[i, i])) % 2 == 0 for i in range(len(G))),
                unimodular=round(float(np.linalg.det(G))) == 1,
                positive=bool(np.all(np.linalg.eigvalsh(G) > 0)),
                irreducible=_connected(G),
                minimal_rank=(rank == 8))


def _coherence_rows():
    """Every (parent, G, F) of the coherence scan's menu with its six gate booleans."""
    alg = EU.ALGEBRAS
    rank8 = [p for p, d in alg.items() if d["rank"] == 8]
    subs = [a for a, d in alg.items() if d["rank"] <= 8 and a not in ("D3", "C2")]
    rows = []
    for p in rank8:
        pd, cP = alg[p], EU.central_charge(p)
        for g, f in itertools.combinations(subs, 2):
            gd, fd = alg[g], alg[f]
            if gd["rank"] + fd["rank"] > pd["rank"]:
                continue
            bridge = gd["fund_dim"] * fd["fund_dim"]
            rows.append(((p, g, f), {
                "(a) c_coset = 0": cP - EU.central_charge(g) - EU.central_charge(f) == 0,
                "(b) h_bridge = 1": EU.h_fund(g) + EU.h_fund(f) == 1,
                "(c) every mode assigned": gd["dim"] + fd["dim"] + bridge == pd["dim"],
                "(d) bridge > 0": bridge > 0,
                "(e) self adjoint parent": pd["fund_dim"] == pd["dim"],
                "(f) trivial centres": gd["centre_order"] == 1 and fd["centre_order"] == 1,
            }))
    return rows


def menus():
    """Each stage: its candidates with gate booleans, and each gate's stated premise (or None)."""
    lat = [((n, r), _lattice_facts(n, G, r)) for n, G, r in W.LATTICES]
    pairs = [((n,), {"emergent Z3": z, "non-pointed": np_}) for n, z, np_, _ in W.PAIRS]
    levels = [((k,), {"modular extension": k % 3 == 0, "h_J = 1": k == 3}) for k in range(1, 13)]
    premises = {
        "integral": "R1: integrality ensures well defined structure constants",
        "even": "R1: evenness ensures locality (integral operator dimensions)",
        "unimodular": "R1: unimodularity ensures every coherence channel has a partner",
        "positive": "R1: positive definiteness excludes negative norm states",
        "irreducible": "R1: a sum of two lattices carries a free relative normalisation",
        "minimal_rank": ("R1: minimality excludes reducible lattices whose relative OPE "
                         "normalizations would be free parameters"),
        "emergent Z3": "quantum numbers crystallise on breaking: an inherited Z3 is written on the vacuum",
        "non-pointed": "the quantum dimensions are read from a nonpointed sector",
        "modular extension": "the simple current extends the theory only at integer h_J",
        "h_J = 1": "R1: h_J >= 2 is a W generator with a free normalisation",
        "(a) c_coset = 0": "every mode is assigned: no coset sector is left over",
        "(b) h_bridge = 1": "R1: the bridge is a current, not a W generator",
        "(c) every mode assigned": "every mode is assigned: gauge + matter + bridge",
        "(d) bridge > 0": "the closure carries matter",
        "(e) self adjoint parent": "self-presence: the parent is its own deepest representation",
        "(f) trivial centres": "quantum numbers crystallise on breaking, not written on the vacuum",
    }
    return {"lattice": lat, "pair": pairs, "level": levels, "coherence": _coherence_rows()}, premises


def _survivors(rows, dropped=()):
    return [key for key, gates in rows if all(v for g, v in gates.items() if g not in dropped)]


def audit():
    stages, premises = menus()
    out = {}
    for stage, rows in stages.items():
        base = _survivors(rows)
        names = list(rows[0][1])
        for g in names:
            rivals = [k for k in _survivors(rows, (g,)) if k not in base]
            status, partner = "selects", None
            if not rivals:
                status = "consequence"
                for h in names:
                    if h != g and [k for k in _survivors(rows, (g, h)) if k not in _survivors(rows, (h,))]:
                        status, partner = "shares", h
                        break
            out[(stage, g)] = dict(survivors=base, rivals=rivals, status=status, partner=partner,
                                   premise=premises[g])
    return out


def _silent(*a, **k):
    pass


def run(report=print):
    res = audit()
    for (stage, g), r in res.items():
        tail = (f"opens {len(r['rivals'])}: {', '.join('/'.join(map(str, x)) for x in r['rivals'][:4])}"
                if r["rivals"] else (f"with {r['partner']}" if r["partner"] else "removes no rival"))
        report(f"  {stage:9s} {g:26s} {r['status']:11s} {tail}")
        report(f"  {'':9s} {'':26s} premise: {r['premise'] or 'none stated in its module'}")
    working = [k for k, r in res.items() if r["status"] != "consequence"]
    unstated = [k for k in working if res[k]["premise"] is None]
    assert not unstated, unstated
    report(f"RESULT: {len(res)} gates, {len(working)} of them select or share the selection, "
           "and every gate names its premise.")
    return dict(gates=res, working=working, unstated=unstated)


if __name__ == "__main__":
    run()
